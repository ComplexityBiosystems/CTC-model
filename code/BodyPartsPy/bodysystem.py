"""
A class for a body system.

A body system is simply a set of body parts
that have already been fitted and verified.

Francesc Font-Clos
May 2018
"""
from .bodybridge import BodyBridge
from .bodypart import BodyPart

import os
from os.path import basename
from os.path import splitext
from os.path import join
from itertools import combinations
import glob

import numpy as np
import pandas as pd

import networkx as nx
from trimesh.ray import ray_pyembree

from scipy.sparse import csc_matrix
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import inv as sparse_inv

from .utils import infer_pressure

from typing import Union
from typing import List


class BodySystem:
    """
    A class for a body system.

    Lots of info.

    """

    def __init__(self, FJcodes=None, pickled_dir=None):
        """
        Create a BodySystem object.

        Using a list of pickled BodyPart objects,
        create a BodySystem object.

        """
        assert pickled_dir is not None
        if FJcodes is None:
            FJcodes = [splitext(basename(x))[0]
                       for x in glob.iglob(join(pickled_dir, "*.p"))]

        # load the bodyparts
        bps = []
        nobps = []
        for FJcode in FJcodes:
            path = os.path.join(pickled_dir, FJcode + ".p")
            bp = pd.read_pickle(path)
            if bp.condition == "accepted":
                bps.append(bp)
            else:
                nobps.append(bp)
        self.bodyparts = bps
        self.unconfirmed_bodyparts = nobps

        # load the pyembree ray intersectors
        for bp in bps:
            bp.mesh.ray = ray_pyembree.RayMeshIntersector(bp.mesh)

        # initiate the networkx nets for all bodyparts
        # initiate a higher-levle netowrkx network where
        # nodes are FJcodes and edges are accepted bridges
        self.supernet = nx.Graph()
        for bp in self.bodyparts:
            self.supernet.add_node(bp.FJcode)
            bp._create_network()

        # join all networks
        self.network = nx.algorithms.operators.compose_all(
            [bp.network for bp in self.bodyparts])

        # set status to not solved
        self.solved = False

    def get_bodypart(self, FJcode):
        """
        Retrieve a bodypart via FJcode.

        Parameters
        ----------
        FJcode: str
            The FJcode of the bodypart to be retrieved.

        Returns
        -------
        bp: bodypartspy.bodypart.BodyPart
            The requested bodypart such that bp.FJcode == FJcode

        Raises
        ------
        Error if none or more than one are found.

        """
        # append _submesh0 if necessary
        if FJcode.find("submesh") == -1:
            FJcode += "_submesh0"
        selected_bps = [x for x in self.bodyparts if x.FJcode == FJcode]
        if len(selected_bps) != 1:
            raise ValueError(f"The requested bodypart {FJcode} does not exist or is"
                             "not unique.")
        return selected_bps[0]

    def get_bodyparts(self, FJcodes: List[str]) -> List[BodyPart]:
        """
        Retrieve a list of bodyparts via FJcode.

        Parameters
        ----------
        FJcode: list[str]
            List of FJcodes of the bodyparts to be retrieved.

        Returns
        -------
        selected_bps: list[BodyPartsPy.bodypart.BodyPart]
            The requested bodyparts.

        Raises
        ------
        Error if none or more than one are found.

        """
        selected_bps = [self.get_bodypart(FJcode) for FJcode in FJcodes]
        return selected_bps

    def find_bridges(self):
        """Find all bp-to-bp connections."""
        bridges = []
        for bp0, bp1 in combinations(self.bodyparts, 2):
            if any(bp0.mesh.contains(bp1.mesh.vertices)) or \
               any(bp1.mesh.contains(bp0.mesh.vertices)):
                bridges.append(BodyBridge(bp0, bp1))
        self.bodybridges = bridges

    def accept_all_bridges(self):
        """Add edges corresponding to all bridges."""
        for bridge in self.bodybridges:
            # add edge to network
            u = f"{bridge.bp0.FJcode}_node{bridge.candidate0}"
            v = f"{bridge.bp1.FJcode}_node{bridge.candidate1}"
            self.network.add_edge(u, v)
            # add edge to the supernetwork (FJcode level)
            self.supernet.add_edge(bridge.bp0.FJcode, bridge.bp1.FJcode)
            bridge.mark_as("accepted")

    def fill_missing_radiuses(self):
        """
        Fill missing values in radius edge attribute.

        This function infers missing radius information
        from first neighbours.
        """
        net = self.network
        for bridge in self.bodybridges:
            if bridge.condition == "accepted":
                u = bridge.bp0.FJcode + "_node" + str(bridge.candidate0)
                v = bridge.bp1.FJcode + "_node" + str(bridge.candidate1)
                nn_edges = list(net.edges([u, v]))
                # remove edges that cross bodyparts
                nn_edges = [(uu, vv) for uu, vv in nn_edges
                            if uu.split("_")[0] == vv.split("_")[0]]
                r = np.mean([net.get_edge_data(*ne)["radius"]
                             for ne in nn_edges])
                nx.set_edge_attributes(net, {(u, v): {"radius": r}})

    def compute_edge_lengths(self):
        """
        Copmpute the length of all edges.

        This function is intended to be executed
        after all bridges have been added.
        """
        net = self.network
        node_pos = nx.get_node_attributes(net, "position")
        for u, v in net.edges:
            edge_length = np.sqrt(sum((node_pos[u] - node_pos[v])**2))
            nx.set_edge_attributes(net, {(u, v): {"length": edge_length}})

    def solve_flows(self,
                    root: Union[List[str], str] = None,
                    viscosity=1.,
                    iterative_min_rad=1,
                    bin_s=0.3,
                    cc_idx=0,
                    force_root: bool = False
                    ):
        """
        Solve the flows problem.

        Parameters
        ----------
        root: str
            Label of root node in the format "FJ????_submesh?_node?".
        viscosity: float
            The viscosity of blood, relative to water = 1.
        iterative_min_rad: float
            In iterative solving mode, sets the radius below wich nodes'
            preassure is *not* updated iteratively.
        bin_s: float
            In iterative solving mode, radius tolerance to infer pressures.
        cc_idx: int
            Which connected compoent (from large to small) to use
        force_root: bool
            If passed source node is not a viable root, force it to be by
            adding a phantom node.
        """
        # root must be a list of nodes
        assert root is not None
        if isinstance(root, str):
            root = [root]
        # convert viscosity to (mmHg s) units
        self.viscosity = viscosity * 1e-6 * 7.5006
        # get largest connected component
        cc = nx.connected_component_subgraphs(self.network)
        sorted_cc = sorted(cc, key=len, reverse=True)
        net = sorted_cc[cc_idx]

        # if root has not degree 1, force it
        if force_root:
            _root = []
            for roo in root:
                if net.degree[roo] > 1:
                    sroo = "source_"+roo
                    net.add_node(sroo)
                    r = np.mean([
                        net[roo][v]["radius"]
                        for v in net.neighbors(roo)
                    ])
                    net.add_edge(sroo, roo, radius=r, length=1)
                    _root.append(sroo)
                else:
                    _root.append(roo, )
            root = _root
            self.fill_missing_radiuses()

        # define internal and external nodes
        internal_nodes = [node for node, k in net.degree if k > 1]
        external_nodes = [node for node, k in net.degree if k == 1]
        for roo in root:
            if roo not in external_nodes:
                raise RuntimeError(
                    f"Node {roo} was passed as source node but is not a leaf")

        # index to label conversions for nodes
        nodes_idx2label = dict(enumerate(net.nodes))
        nodes_label2idx = {v: k for k, v in nodes_idx2label.items()}

        # sparse matrices to solve the linear algebra problem
        N = len(net.nodes)
        b = dok_matrix((N, 1))
        S = dok_matrix((N, N))

        # initial condition
        p = np.zeros(N) + 20
        # set pressure of root
        for roo in root:
            p[nodes_label2idx[roo]] = 100
        # if not first iteration, guess boundary conditions
        if self.solved:
            # auxiliar df to infer pressures
            p_dict = nx.get_node_attributes(self.network_solved, "pressure")
            r_dict = nx.get_node_attributes(self.network, "radius")
            df = pd.DataFrame(data=[r_dict, p_dict],
                              index=["radius", "pressure"]).T.dropna()
            # set boundary conditions of external nodes except root
            for node in external_nodes:
                if node not in root:
                    i = nodes_label2idx[node]
                    target_r = df.loc[node, "radius"]
                    if target_r > iterative_min_rad:
                        inferred_p = infer_pressure(df=df,
                                                    target_r=target_r,
                                                    bin_s=bin_s)
                        p[i] = inferred_p

        # fill in coefficients
        for u, v in net.edges:
            r = net[u][v]["radius"]
            length = net[u][v]["length"]
            cij = np.math.pi * (2 * r) ** 4 / (128 * self.viscosity * length)
            i = nodes_label2idx[u]
            j = nodes_label2idx[v]
            if u in internal_nodes:
                S[i, i] += cij
                S[j, i] += -cij
            else:
                b[j] += cij * p[i]
            if v in internal_nodes:
                S[j, j] += cij
                S[i, j] += -cij
            else:
                b[i] += cij * p[j]
        bb = [x in internal_nodes for x in net.nodes]
        SS = csc_matrix(S.todense()[bb].T[bb].T)
        # solve the system
        # i.e. find pressure at each node
        SS_inv = sparse_inv(SS)
        sol = SS_inv.dot(csc_matrix(b.todense()[bb]))
        pressure_internal = dict(zip(internal_nodes, (sol.toarray().T[0])))
        pressure_external = dict(zip(external_nodes,
                                     p[[x in external_nodes
                                        for x in net.nodes]]))
        pressure_all_nodes = {**pressure_internal, **pressure_external}
        # fill in network with pressure and flows
        dinet = nx.DiGraph()
        dinet.add_nodes_from(net.nodes)
        # copy over node positions
        dd = nx.get_node_attributes(net, "position")
        nx.set_node_attributes(dinet, dd, "position")
        # write pressures dict
        nx.set_node_attributes(dinet, pressure_all_nodes, "pressure")
        # mark the system as solved
        self.solved = True

        edge_attrs = {}
        pressure_dict = nx.get_node_attributes(dinet, "pressure")
        for i, (u, v) in enumerate(net.edges):
            # get its properties
            pu = pressure_dict[u]
            pv = pressure_dict[v]
            r = net[u][v]["radius"]
            length = net[u][v]["length"]
            cij = np.math.pi * (2 * r) ** 4 / (128 * self.viscosity * length)
            flow = cij * (pu - pv)
            if flow < 0:
                flow = -flow
                u, v = v, u
            speed = flow / (np.math.pi * r**2)
            time = length / speed
            # create the edge
            dinet.add_edge(u, v)
            # save properties in dict entry
            edge_attrs[(u, v)] = {"flow": flow,
                                  "speed": speed,
                                  "time": time,
                                  "radius": r,
                                  "length": length
                                  }
        nx.set_edge_attributes(dinet, edge_attrs)
        self.network_solved = dinet

    def compute_node_radius(self):
        """
        Compute the radius of nodes.

        The radius of a node is defined as the average of the radius of its
        input/output edges.

        """
        fradius = {}
        for bp in self.bodyparts:
            bpnet = bp.network
            for node in bpnet.nodes:
                radius = np.mean([bpnet[node][nn]["radius"]
                                  for nn in nx.neighbors(bpnet, node)])
                fradius[node] = radius
        nx.set_node_attributes(self.network, fradius, "radius")

    def compute_node_flows(self):
        """
        Compute the total flow going throught nodes.

        As the sum of the flow of its outgoing edges.

        """
        flow_dict = {}
        net_s = self.network_solved
        for node in net_s.nodes:
            flow_out = np.sum([net_s[e[0]][e[1]]["flow"]
                               for e in net_s.out_edges(node)])
            flow_in = np.sum([net_s[e[0]][e[1]]["flow"]
                              for e in net_s.in_edges(node)])
            if np.isclose(flow_in, flow_out):
                flow_dict[node] = flow_in
            elif flow_in == 0:
                flow_dict[node] = flow_out
            elif flow_out == 0:
                flow_dict[node] = flow_in

        nx.set_node_attributes(net_s, flow_dict, "flow")
