"""
Functions to convert between meshes and graphs.

Mostly making use of Elastic Principal Trees.
ADD REFERENCES

Francesc Font-Clos
April 2018

"""
from .sampling import get_sample_mixed
from .sampling import get_sample_surface
from .sampling import get_sample_inside_naive
from .utils import join_graphs
from .utils import _get_node_sphere
from .utils import _get_edge_cylinder
from .utils import _can_split_data
from .utils import _split_data

import numpy as np

import os
import sys
# this is a stupid hack to get the level where
# I usually place al git repos
ROOTPATH = os.path.abspath(".").split("MetastasisModel")[0]
sys.path.append(os.path.join(ROOTPATH))


def mesh2graph(mesh,
               n_mesh_sample=5000, n_nodes=25,
               cut=0, sampling="surface", **kwargs):
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
    cut : int
        How many times to split the mesh using PCA.
        Defaults to zero, max is 2 (thus 4 pieces).
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
        return _get_tree_onepiece(points, n_nodes=n_nodes, **kwargs)

    # check that we can actually cut the data
    # before doing it
    if not _can_split_data(mesh, points):
        return mesh2graph(mesh=mesh, n_mesh_sample=n_mesh_sample,
                          n_nodes=n_nodes, cut=0, sampling=sampling, **kwargs)
    # split the data into two chunks
    points1, points2, v = _split_data(mesh, points)

    # if no more cuts asked, lets start joining graphs
    if cut == 1:
        # get two trees
        nodes1, edges1 = _get_tree_onepiece(
            points=points1, n_nodes=n_nodes, **kwargs)
        nodes2, edges2 = _get_tree_onepiece(
            points=points2, n_nodes=n_nodes, **kwargs)
        # join the two trees
        nodes, edges = join_graphs(nodes1, edges1, nodes2, edges2, v)
        return nodes, edges

    # again, check that we can split before doing it
    if (not _can_split_data(mesh, points1)) or\
       (not _can_split_data(mesh, points2)):
        return mesh2graph(mesh=mesh, n_mesh_sample=n_mesh_sample,
                          n_nodes=n_nodes, cut=1, sampling=sampling, **kwargs)

    # split the data further
    points11, points21, v1 = _split_data(mesh, points1)
    points12, points22, v2 = _split_data(mesh, points2)

    if cut == 2:
        # get four trees
        nodes11, edges11 = _get_tree_onepiece(
            points=points11, n_nodes=n_nodes, **kwargs)
        nodes21, edges21 = _get_tree_onepiece(
            points=points21, n_nodes=n_nodes, **kwargs)
        nodes12, edges12 = _get_tree_onepiece(
            points=points12, n_nodes=n_nodes, **kwargs)
        nodes22, edges22 = _get_tree_onepiece(
            points=points22, n_nodes=n_nodes, **kwargs)
        # join the four smallest trees into 2
        nodes1, edges1 = join_graphs(nodes11, edges11, nodes21, edges21, v1)
        nodes2, edges2 = join_graphs(nodes12, edges12, nodes22, edges22, v2)
        # join the two medium trees into 1
        nodes, edges = join_graphs(nodes1, edges1, nodes2, edges2, v)
        return nodes, edges

    else:
        raise ValueError("Cut must be 0, 1, or 2")


def graph2mesh(nodes, edges, node_radius=0.3, edge_radius=0.1):
    """
    Turn a graph into a mesh.

    This function transforms a graph into a mesh by drawing
    nodes as spheres and edges as cylinders.

    Parameters
    ----------
    nodes : np.array, (n_nodes, 3)
        Node positions.
    edges : np.array, (n_edges, 2)
        Indices defining the edges.
    node_radius : float
        Radius for spheres representing nodes.
    edge_radius : float
        Radius for cylinders representing edges.

    Returns
    -------
    scene : trimesh.base.Trimesh
        Mesh representing the input graph.

    """
    scene = _get_node_sphere(nodes[0], radius=node_radius)
    for node in nodes[1:]:
        scene += _get_node_sphere(node, radius=node_radius)
    for i, j in edges:
        scene += _get_edge_cylinder(nodes[i], nodes[j], radius=edge_radius)
    return scene


def _get_tree_onepiece(points, n_nodes=25, engine="python", **kwargs):
    """
    Get a tree from a set of points.

    By default uses the python implementation.
    """
    if "opoints" in kwargs.keys():
        _points = kwargs["opoints"]
    else:
        _points = points
    if engine == "python":
        return _get_tree_onepiece_python(points=_points, n_nodes=n_nodes,
                                         **kwargs)
    elif engine == "R":
        return _get_tree_onepiece_R(points=_points, n_nodes=n_nodes)
    else:
        return ValueError(f"I do not know the {engine} engine")


def _get_tree_onepiece_python(points, n_nodes=25, **kwargs):
    """
    Get a tree from a set of points.

    Wrapping around computeElasticPrincipalGraph
    """
    from ElPiGraph.computeElasticPrincipalGraph \
        import computeElasticPrincipalGraph
    nodes, matrix = computeElasticPrincipalGraph(
        points, NumNodes=n_nodes,
        drawPCAview=False, drawAccuracyComplexity=False,
        drawEnergy=False, MaxNumberOfIterations=10, eps=0.01,
        verbose=False)
    # get edges in u,v format
    edges = []
    for u, row in enumerate(matrix):
        for v, element in enumerate(row):
            if element > 0 and u < v:
                edges.append([u, v])
    edges = np.array(edges)
    return nodes, edges


def _get_tree_onepiece_R(points, n_nodes=25):
    """
    Get a tree from a set of points.

    Wrapping around ElPiGraph.R computeElasticPrincipalTree.
    """
    from rpy2.robjects.packages import importr
    from rpy2.robjects import numpy2ri
    elpi = importr("ElPiGraph.R")

    numpy2ri.activate()
    tmp = elpi.computeElasticPrincipalTree(X=points, NumNodes=n_nodes,
                                           drawAccuracyComplexity=False,
                                           drawEnergy=False,
                                           drawPCAView=False,
                                           verbose=False
                                           )
    numpy2ri.deactivate()
    nodes = np.array(tmp[0][0])
    edges = np.array(tmp[0][1][0]) - 1

    return nodes, edges
