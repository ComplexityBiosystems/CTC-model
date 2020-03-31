"""
Quality Control.

Functions to check quality of automatized graph extraction.

Francesc Font-Clos
May 2018

"""
from .utils import _sample_random_edge_point
from .utils import get_radius

import numpy as np
from scipy.spatial.distance import norm


def check_nodes_inside(mesh, nodes, edges):
    """
    Check if all nodes are inside the mesh.

    Parameters
    ----------
    mesh : trimesh.base.Trimesh
        The mesh used to create a graph.
    nodes : np.array, (n_nodes, 3)
        Positions of the graph nodes.
    edges : np.array, (n_edges, 2)
        Indices that define the edges.

    Returns
    -------
    nodes_inside : bool
        Whether all nodes are inside the mesh.

    """
    nodes_inside = all(mesh.contains(nodes))
    return nodes_inside


def check_edges_inside(mesh, nodes, edges):
    """
    Check if all edges are inside the mesh.

    Parameters
    ----------
    mesh : trimesh.base.Trimesh
        The mesh used to create a graph.
    nodes : np.array, (n_nodes, 3)
        Positions of the graph nodes.
    edges : np.array, (n_edges, 2)
        Indices that define the edges.

    Returns
    -------
    edges_inside : bool
        Whether all nodes are inside the mesh.

    """
    edge_points = np.concatenate([_sample_random_edge_point(nodes, edges)
                                  for _ in range(100)])
    edges_inside = all(mesh.contains(edge_points))
    return edges_inside


def check_emptyness(mesh, nodes, edges):
    """
    Try to guess if there are empty spaces.

    Computes the distance of each vertex to its closest
    graph node, rescaled by average mesh radius. Assumes
    something is wrong if there are vertex with dist > th

    Parameters
    ----------
    ...

    """
    radiuses = [get_radius(mesh, nodes, edge) for edge in edges]
    mean_rad = np.mean(radiuses)
    vertices_dist_to_closest_node = [min(norm(nodes - v, axis=1))
                                     for v in mesh.vertices]
    vertices_dist_to_closest_node /= mean_rad
    return max(vertices_dist_to_closest_node)
