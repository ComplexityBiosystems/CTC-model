"""
BodyPartPy Visualization.

Functions to visualize meshes, graphs, etc

"""
from .bodysystem import BodySystem

from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances
import requests

from trimesh.base import Trimesh
import numpy as np
from typing import Sequence
from typing import List
from typing import Union
# from scipy.spatial.distance import norm
from trimesh.ray.ray_pyembree import RayMeshIntersector

import networkx as nx


def random_hex_colors(n_colors=10):
    """Get some random hex colors."""
    res = requests.get(f'http://www.colr.org/json/colors/random/{n_colors}')
    data = res.json()['colors']
    colors_hex = [d["hex"] for d in data]
    return colors_hex


def view_pca(mesh, nodes, edges, ax=None, pca=None):
    """
    View a (mesh, nodes, edges) set in PCA space.

    Parameters
    ----------
    mesh : trimesh.base.Trimesh
        The original mesh.
    nodes : np.array, (n, 3)
        The nodes positions.
    edges : np.array, (n, 2)
        Pairse of indices (i, j) definiing the
        graphss' edges
    ax : matplotlib.axes.Axes
        Where to plot the pca projection.
    pca : sklearn.decomposition.PCA
        Fitted pca instance.

    """
    if pca is None:
        pca = PCA(n_components=2).fit(mesh.vertices)

    # plot mesh
    ax.scatter(*pca.transform(mesh.sample(10000)).T, color="0.7", s=1)

    # plot nodes
    ax.scatter(*pca.transform(nodes).T, color="red", s=10)

    # plot edges
    for i, j in edges:
        ax.plot(*pca.transform(nodes[[i, j]]).T, color="0.2")


class AttributeInterpolator():
    def __init__(
        self,
        mesh: Trimesh,
        att_points: np.ndarray,
    ) -> None:
        assert len(att_points) > 1
        self.mesh = mesh
        self.att_points = att_points
        try:
            distances = euclidean_distances(
                self.mesh.vertices,
                self.att_points
            )
            self.closest_idx = np.argmin(distances, axis=1)
        except MemoryError:
            # do it by chunks
            def chunker(seq, size):
                return (seq[pos:pos + size] for pos in range(0, len(seq), size))

            closest_idx = []
            for vertices_chunk in chunker(self.mesh.vertices, 5000):
                print("# DEBUG: ", vertices_chunk.shape, self.att_points.shape)
                distances = euclidean_distances(
                    vertices_chunk,
                    self.att_points
                )
                closest_idx.append(np.argmin(distances, axis=1))
                del distances
            self.closest_idx = np.concatenate(closest_idx)

    def interpolate(self, att_values: Sequence) -> np.ndarray:
        assert len(att_values) == len(self.att_points)
        return np.array(att_values)[self.closest_idx]


class BodySystemAttributeInterpolator(AttributeInterpolator):
    def __init__(
        self,
        bodysystem: BodySystem,
        FJcodes: Union[List[str], None] = None
    ) -> None:
        """
        Create a bodysystem-aware attribute interpolator given a bodysystem and a list of bodyparts.

        Parameters
        ----------
        bodysystem : BodySystem
            The bodysystem with solved flows. 
        FJcodes : List[str]
            List of bodypart labels.

        Returns
        -------
        AttributeInterpolator
        """
        if FJcodes is not None:
            bodyparts = bodysystem.get_bodyparts(FJcodes)
        else:
            FJcodes = FJcodes = list(set(
                "_".join(x.split("_")[:2])
                for x in bodysystem.network_solved.nodes
                if x[:7] != "source_"
            ))
            bodyparts = bodysystem.get_bodyparts(FJcodes)

        mesh = np.sum([x.mesh for x in bodyparts])
        self.edges_list = list(np.concatenate([
            list(x.network.edges)
            for x in bodyparts
        ]))
        self.network_solved = bodysystem.network_solved
        edges_positions = (np.concatenate(
            [(x.pmesh.vertices) for x in bodyparts]))

        self.attribute_interpolator = AttributeInterpolator(
            mesh=mesh,
            att_points=edges_positions,
        )

    def interpolate_edge_attribute(
        self,
        edge_attribute: Union[str, dict]
    )-> np.ndarray:
        if isinstance(edge_attribute, str):
            tmp = nx.get_edge_attributes(
                self.network_solved,
                edge_attribute
            )
        elif isinstance(edge_attribute, dict):
            tmp = edge_attribute
        else:
            raise RuntimeError(
                "Wrong edge attribute type (must be str or dict")
        raw_att_values = []
        for x in self.edges_list:
            try:
                raw_att_values.append(tmp[(x[0], x[1])])
            except KeyError:
                raw_att_values.append(tmp[(x[1], x[0])])

        raw_att_values = np.concatenate([[x]*64 for x in raw_att_values])
        att_values = raw_att_values
        interpolated_attribute = self.attribute_interpolator.interpolate(
            att_values=att_values
        )
        return interpolated_attribute
