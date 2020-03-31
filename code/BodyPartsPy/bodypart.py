"""
A class for a body part.

Francesc Font-Clos
May 2018
"""
from .io import load_mesh

from .conversions import mesh2graph
from .conversions import graph2mesh

from .utils import prune_graph
from .utils import get_radius
from .utils import _get_edge_plate

from .qc import check_emptyness
from .qc import check_edges_inside
from .qc import check_nodes_inside

import os
import pandas as pd
import numpy as np
import networkx as nx


class BodyPart(object):
    """
    A class for a body part.

    Lots of other info.
    """

    def __init__(self, FJcode, where, use_embree=False,
                 ignore_disconnected_objects=False):
        """
        Create BodyPart object.

        Parameters
        ----------
        a

        """
        mesh = load_mesh(FJcode, where=where, use_embree=use_embree)
        if len(mesh.split()) != 1 and ignore_disconnected_objects is False:
            raise ValueError(
                f"Mesh {FJcode} has more than one connected component.\n"
                f"Loaded from '{where}'.")

        self.mesh = mesh
        self.FJcode = FJcode
        self.where = where
        self.mark_as("not fitted")

    # def __getstate__(self):
    #     ### not used anymore, set use_embree=False ###
    #     """Delete properties that cannot be pickled."""
    #     odict = self.__dict__.copy()
    #     for mesh_kwd in ["mesh", "gmesh", "pmesh"]:
    #         del odict[mesh_kwd].visual
    #         del odict[mesh_kwd].ray
    #     return odict

    def fit(self,
            n_nodes=25,
            PCAsplit=0,
            n_mesh_sample=5000,
            sampling="surface",
            prune=True,
            edges=None,
            nodes=None,
            **kwargs
            ):
        """
        Fit a graph to the mesh.

        Uses Elastic Principal Graph python's implementation by default.

        Parameters
        ----------
        n_nodes :  int
            Number of nodes.
        PCAsplit : int
            Number of PCA splits to perform before fitting.
        n_mesh_sample : int
            Number of points to sample the mesh.
        sampling : str
            see other docstrings.

        """
        # fit the graph
        if edges is None or nodes is None:
            nodes, edges = mesh2graph(self.mesh,
                                      n_mesh_sample=n_mesh_sample,
                                      n_nodes=n_nodes,
                                      cut=PCAsplit,
                                      sampling=sampling,
                                      **kwargs
                                      )
        # prune the graph
        if prune:
            nodes, edges = prune_graph(nodes, edges)
        else:
            edges = np.array(edges)
        gmesh = graph2mesh(nodes, edges)

        # save data
        self.fit_params = {"n_nodes": n_nodes,
                           "cut": PCAsplit,
                           "n_mesh_sample": n_mesh_sample,
                           "sampling": sampling}
        self.nodes = nodes
        self.edges = edges
        self.gmesh = gmesh
        self.n_nodes = n_nodes
        self.n_edges = len(edges)

        # quality control
        self.max_emptyness_ratio = check_emptyness(self.mesh, self.nodes,
                                                   self.edges)
        self.nodes_inside = check_nodes_inside(self.mesh, self.nodes,
                                               self.edges)
        self.edges_inside = check_edges_inside(self.mesh, self.nodes,
                                               self.edges)
        # radiuses
        radiuses = [get_radius(self.mesh, self.nodes, edge)
                    for edge in self.edges]
        self.radiuses = radiuses

        # edge plates
        for i, ((u, v), r) in enumerate(zip(edges, radiuses)):
            if i == 0:
                scene = _get_edge_plate(nodes[u], nodes[v], r)
            else:
                scene += _get_edge_plate(nodes[u], nodes[v], r)
        self.pmesh = scene
        self.mark_as("not verified")

    def mark_as(self, condition):
        """
        Mark the object according to verifycation degree.

        Parameters
        ----------
        condition : str
            In which condition the BodyPart object is in.

            If the object has not been fitted, use "not fitted".

            If the object has been fitted but not inspected,
            use "not verified".

            If the fit is OK, use "accepted".

            If the fit is NOT OK, use one of the following:
            "increase nodes", "decrease nodes", "increase cut",
            "decrease cut" or a generic "not accepted".

        """
        assert condition in ["accepted", "not accepted",
                             "increase nodes", "decrease nodes",
                             "increase cut", "decrease cut",
                             "not fitted", "not verified"
                             ]
        self.condition = condition

    def to_pickle(self, where):
        """
        Save bodypary in pickled form.

        This function will save a pickled copy of self
        into a given path. The filename is constructed
        using self.FJcode.

        Parameters
        ----------
        path : str
            Where to save the object.

        """
        pd.to_pickle(self, os.path.join(where, self.FJcode + ".p"))

    def _create_network(self):
        assert self.condition == "accepted"
        nx_edges = [[self.FJcode + "_node" + str(u),
                     self.FJcode + "_node" + str(v)]
                    for u, v in self.edges]
        net = nx.from_edgelist(nx_edges)

        # set nodes attributes
        nodes_position_dict = {self.FJcode + "_node" + str(a): b
                               for a, b in enumerate(self.nodes)}
        nx.set_node_attributes(net, nodes_position_dict, name="position")

        # set edge attribute
        edges_radiuses_dict = dict(zip(net.edges, self.radiuses))
        nx.set_edge_attributes(net, edges_radiuses_dict, name="radius")

        self.network = net
