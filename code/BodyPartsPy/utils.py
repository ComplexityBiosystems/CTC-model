"""
Various utility functions.

Francesc Font-Clos
May 2018
"""
import numpy as np
import pandas as pd
from scipy.spatial.distance import norm
import networkx as nx
from sklearn.decomposition import PCA

from trimesh.creation import icosphere
from trimesh.transformations import rotation_matrix
from trimesh.creation import cylinder


def remove_superflous_node(net, node):
    r"""
    Remove a superflous node from a network.

    A superflous node $s$ is defined simply as
    having degree two. This function removes it
    while preserving the topology of the network:

    + remove the edge (u, s)
    + remove the edge (s, v)
    + add a new edge (u, v)

    """
    (u, du), (v, dv) = net[node].items()
    l0 = du["length"]
    l1 = dv["length"]
    r0 = du["radius"]
    r1 = dv["radius"]
    l = l0 + l1
    # resistors in series
    r = ((l0 + l1) / ((l0 / r0**4) + (l1 / r1**4))) ** (1/4.)
    net.remove_edge(u, node)
    net.remove_edge(v, node)
    net.add_edge(u, v, radius=r, length=l)
    net.remove_node(node)


def infer_pressure(df, target_r, bin_s=0.1):
    """
    Infer the pressure at given radius.

    Computes the average of the pressure for nodes
    within certain tolerance of a target radius.

    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe containing the radius and pressure of all nodes.
    target_r: float
        Value of radius for which we want to infer the pressure.
    bin_s: float
        Relative radius tollerance when binning nodes. Nodes whose radius is
        within (1 +/- bin_s) * target_r are used to infer target_r.

    Returns
    -------
    target_p: float
        Value of pressure inferred from the given data.

    """
    # make sure we have values for the r target
    assert target_r >= df.radius.min()
    assert target_r <= df.radius.max()
    # slice the df within +-10% of the target
    c2 = (df.radius < (1 + bin_s) * target_r)
    c1 = (df.radius >= (1 - bin_s) * target_r)
    df_slice = df.loc[c1 & c2]
    inferred_p = df_slice.pressure.mean()
    return inferred_p


def infer_radius(net, u, v):
    """
    Inferr the radius of an edge that crosses bodyparts.
    """
    nn_edges = list(net.edges([u, v]))
    # remove edges that cross bodyparts
    nn_edges = [(uu, vv) for uu, vv in nn_edges
                if uu.split("_")[0] == vv.split("_")[0]]
    r = np.mean([net.get_edge_data(*ne)["radius"]
                 for ne in nn_edges])
    return r


def prune_graph(nodes, edges, max_antenna_remove=1):
    """
    Remove one-edge antennas from a graph.

    Deletes all edges that are just length-one
    antenna, without breaking the graph.

    Parameters
    ----------
    nodes : np.array, (n_nodes, 3)
        Positions of the nodes.
    edges : np.array, (n_edges, 2)
        Indices defining the edges.

    """
    # turn into proper graph
    G = nx.Graph([list(x) for x in edges])
    nG, _ = _prune_graph_nx(G, max_antenna_remove=max_antenna_remove)

    index_mapping = dict([[b, a] for a, b in list(enumerate(list(nG.nodes)))])
    nodes_prunned = nodes[list(nG.nodes)]
    edges_prunned = np.array([[index_mapping[a], index_mapping[b]]
                              for a, b in nG.edges])
    return nodes_prunned, edges_prunned


def _prune_graph_nx(G, max_antenna_remove=1):
    # hubs are nodes where branching occurs
    hubs = [k for k, v in dict(G.degree).items() if v > 2]
    # subgraphs after removing hubs
    H = G.copy()
    H.remove_nodes_from(hubs)
    subgraphs = list(nx.connected_component_subgraphs(H))

    # make prunned graph
    nG = G.copy()
    removed_nodes = []
    for subgraph in subgraphs:
        # remove nodes only if the resulting graph remains in one piece
        if subgraph.number_of_nodes() <= max_antenna_remove:
            tmpG = nG.copy()
            tmpG.remove_nodes_from(subgraph.nodes)
            if len(list(nx.connected_components(tmpG))) == 1:
                nG.remove_nodes_from(subgraph.nodes)
                removed_nodes.append(subgraph.nodes)
    return nG, removed_nodes


def get_radius(mesh, nodes, edge):
    """
    Compute the radius of the mesh at a given edge.

    Given the mesh and a set of nodes, edges obtained
    e.g. via get_graph(mesh, ...), for a given edge
    u, v we approximate the radius as follows:
    1. c = (u+v)/2
    2. n = v-u
    2. cut mesh throuch c with normal n
    3. compute length of intersection (only one CC)
    4. infer radius

    This is a rought approximation but seems the best we can
    do in this situation. Also, a rought approximation will
    be more than enough for our needs.

    """
    u, v = nodes[edge]
    section = mesh.section(v - u, (u + v) / 2)

    connected_comps = list(nx.connected_components(section.vertex_graph))
    cut_centers = [section.vertices[list(H)].mean(axis=0)
                   for H in connected_comps]
    mycut_idx = np.argmin(norm(cut_centers - v, axis=1))

    lengths = [norm(section.vertices[line.end_points[0]] -
                    section.vertices[line.end_points[1]])
               for line in section.entities
               if (line.end_points[0] in connected_comps[mycut_idx])]

    radius = np.sum(lengths) / (2 * np.math.pi)
    return radius


def join_graphs(nodes1, edges1, nodes2, edges2, vertex):
    """
    Join two graphs using an auxiliary vertex.

    At the moment simply look for the closest node
    of each tree to the aux vertex, and link.
    """
    closest_n1 = np.argmin(norm(nodes1 - vertex, axis=1))
    closest_n2 = np.argmin(norm(nodes2 - vertex, axis=1))
    nodes = np.concatenate([nodes1, nodes2, [vertex]])
    edges = np.concatenate([edges1, edges2 + len(nodes1),
                            [[len(nodes) - 1, closest_n1]],
                            [[len(nodes) - 1, len(nodes1) + closest_n2]]])
    return nodes, edges


def _can_split_data(mesh, points):
    """Check if data can be safely splitted."""
    pca = PCA(n_components=3).fit(points)
    plane_normal = pca.components_[0]
    plane_origin = pca.inverse_transform([0, 0, 0])
    section = mesh.section(plane_normal, plane_origin)
    G = nx.Graph()
    G.add_edges_from([line.end_points for line in section.entities])
    # make sure we only cutted one leg
    # if not, return no cut version
    if not len(list(nx.connected_components(G))) == 1:
        return False
    else:
        return True


def _split_data(mesh, points):
    """Split data using PCA0."""
    if not _can_split_data(mesh, points):
        raise ValueError("Please verify data can be split beforehand.")
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
    v = section.vertices.mean(axis=0)
    return points1, points2, v


def _sample_random_edge_point(nodes, edges):
    U = nodes[edges.T[0]]
    V = nodes[edges.T[1]]
    t = np.random.uniform(size=len(edges))
    return U + np.array([tt * uv for tt, uv in zip(t, V - U)])


def _get_node_sphere(node, radius=0.5):
    sp = icosphere(subdivisions=2)
    sp.apply_scale(radius)
    sp.apply_translation(node)
    return sp


def _get_edge_cylinder(node0, node1, radius=0.1):
    center = (node0 + node1) / 2
    normal = node1 - node0
    cyl = _get_cylinder(radius=radius, normal=normal,
                        center=center, height=norm(normal), sections=6)
    return cyl


def _angle_between(a, b):
    """Compute angle between two 3d vectors."""
    import math
    from numpy.linalg import norm as npnorm
    arccosInput = np.dot(a, b) / npnorm(a) / npnorm(b)
    arccosInput = 1.0 if arccosInput > 1.0 else arccosInput
    arccosInput = -1.0 if arccosInput < -1.0 else arccosInput
    return math.acos(arccosInput)


def _get_cylinder(radius, normal, center, height=0.1, sections=10):
    """Get a thin disk at a given location and orientation."""
    rotation = rotation_matrix(angle=_angle_between([0, 0, 1], normal),
                               direction=np.cross([0, 0, 1], normal),
                               )
    cyl = cylinder(radius=radius, height=height, transform=rotation,
                   sections=sections)
    cyl.vertices += center
    return cyl


def _get_edge_plate(node0, node1, radius, h=0.1):
    """
    Create a flat cylinder between two nodes.

    The cylinder is placed at the mid-point of the
    edge joining the nodes, with its axis in the
    edge direction. Given an edge and a fitted
    radius, this function can be used to test
    visually if the radius is ok.

    Parameters
    ----------
    node0 : np.array, 3
        First node of an edge.
    node1 : np.array, 3
        Second node of an edges.
    radius : float
        Radius of the cylinder.
    h : float
        Height of the cylinder.

    """
    center = (node0 + node1) / 2
    normal = node1 - node0
    cyl = _get_cylinder(radius=radius, normal=normal,
                        center=center, height=h, sections=32)
    return cyl
