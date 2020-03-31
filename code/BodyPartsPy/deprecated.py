"""Functions I don't want to use anymore."""
import networkx as nx


def mesh2graph(mesh,
               n_mesh_sample=5000, n_nodes=25,
               cut=0, sampling="surface"):
    """
    Get a graph from a mesh.

    Parameters
    ----------
    mesh : trimesh.base.Trimesh
        A mesh we want to convert into a graph.
    n_mesh_sample : int
        Number of samples taken from the mesh to fit.
    n_nodes : int
        Number of nodes of the resulting graph.
    zcut : float
        z-level at which to cut the mesh
    sampling : str
        What part of the mesh to sample.
        Must be one of "surface", "interior", "both"

    """
    # which sampling to use
    assert sampling in ["surface", "interior", "both"]
    get_sample = {"both": get_sample_mixed,
                  "surface": get_sample_surface,
                  "interior": get_sample_inside_naive,
                  }[sampling]

    # get main sample
    points = get_sample(mesh, n_mesh_sample)

    # case one one chunk
    if cut == 0:
        return _get_tree_onepiece(points, n_nodes=n_nodes)

    # split the data into two chunks
    df = pd.DataFrame(points)
    pca = PCA(n_components=3).fit(points)
    pca0 = pca.transform(points).T[0]
    points1 = df.loc[pca0 < 0].values
    points2 = df.loc[pca0 >= 0].values

    # get the split vertex
    plane_normal = pca.components_[0]
    plane_origin = pca.inverse_transform([0, 0, 0])
    section = mesh.section(plane_normal, plane_origin)
    G = nx.Graph()
    G.add_edges_from([line.end_points for line in section.entities])

    # make sure we only cutted one leg
    # if not, return no cut version
    if not len(list(nx.connected_components(G))) == 1:
        return _get_tree_onepiece(points, n_nodes=n_nodes)

    v = section.vertices.mean(axis=0)

    # get two trees
    nodes1, edges1 = _get_tree_onepiece(points=points1,
                                        n_nodes=n_nodes,
                                        )
    nodes2, edges2 = _get_tree_onepiece(points=points2,
                                        n_nodes=n_nodes,
                                        )
    # join the two trees
    nodes, edges = join_graphs(nodes1, edges1, nodes2, edges2, v)

    return nodes, edges


def check_chain(edges):
    """
    Check if the edges define a linear chain.

    Parameters
    ----------
    edges : np.array, (n_edges, 2) or nx.Graph
        Indices defining the edges of a graph,
        or networkx.Graph object

    """
    if isinstance(edges, nx.Graph):
        G = edges
    else:
        G = nx.Graph(([list(x) for x in edges]))

    # number of CC
    num_cc = len(list(nx.connected_components(G)))

    # true if only nodes with degree 1 or 2 exist
    is_chain = set(dict(G.degree).values()) == {1, 2}

    return ((num_cc == 1) and is_chain)
