"""
Uniformly sample meshes.

Francesc Font-Clos
May 2018
"""
import numpy as np
from trimesh.sample import volume_mesh


def get_sample_mixed(mesh, n_points):
    """
    Get a uniform sample of a mesh.

    Samples a mesh combining interior points
    with points in the surface.

    Parameters
    ----------
    mesh : trimesh.base.Trimesh
        The mesh we want to sample.
    n_points : int
        The number of samples

    Notes
    -----
    At the moment, when calling get_sample_inside()
    we don't get n_points, but much less....

    """
    # get_sample_inside = get_sample_inside_naive (DEPRECATED)
    points_inside = get_sample_inside(mesh, n_points)
    points_surface = get_sample_surface(mesh, len(points_inside))
    return np.concatenate((points_inside, points_surface))


def get_sample_inside_naive(mesh, n_points):
    """
    Get a uniform sample of points inside a mesh.

    This is the most naive, brute-force implementation
    possible.
    """
    m = mesh.vertices.min(axis=0)
    M = mesh.vertices.max(axis=0)

    points = m + (M - m) * np.random.uniform(size=(n_points, 3))
    points_inside = points[mesh.contains(points)]
    return points_inside


def get_sample_inside(mesh, n_points):
    """Sample with rejection method."""
    data = volume_mesh(mesh=mesh, count=n_points)
    print("# Interior sampling efficiency: %.2f" % (len(data) / n_points))
    return data


def get_sample_surface(mesh, n_points):
    """Get a uniform sample of points in the surface of a mesh."""
    return np.array(mesh.sample(n_points))
