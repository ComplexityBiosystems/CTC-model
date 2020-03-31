"""
Class to read flows and flow ratios from a solved directed flow graph.

Francesc Font-Clos
Feb 2019
"""

import networkx as nx
import numpy as np
import pandas as pd


class FlowReader:
    def __init__(self, net: nx.DiGraph, source: str):
        assert source in net.nodes
        self.net = net
        self.source = source
        self.total_flow = self._node_out_flow(source)

    @staticmethod
    def from_pickle(path_to_pickled_net: str, source: str) -> "FlowReader":
        net = pd.read_pickle(path_to_pickled_net)
        return FlowReader(net=net, source=source)

    def mesh_flow_ratio(self, FJcode: str) -> float:
        if FJcode.find("_submesh") == -1:
            FJcode += "_submesh0"
        nodes = [
            x for x in self.net.nodes
            if x.find(FJcode) > -1
        ]
        max_in_flow = max((self.node_in_flow_ratio(node)) for node in nodes)
        max_out_flow = max((self.node_out_flow_ratio(node)) for node in nodes)
        return max(max_in_flow, max_out_flow)

    def node_in_flow_ratio(self, node: str) -> float:
        """Relative in flow of a node"""
        in_flow = self._node_in_flow(node)
        return in_flow / self.total_flow

    def node_out_flow_ratio(self, node: str) -> float:
        """Relative out flow of a node"""
        out_flow = self._node_out_flow(node)
        return out_flow / self.total_flow

    def _node_in_flow(self, node: str) -> float:
        """Absolute in flow of a node."""
        in_flow: float = np.sum([
            self.net.edges[edge]["flow"]
            for edge in self.net.in_edges(node)
        ])
        return in_flow

    def _node_out_flow(self, node: str) -> float:
        """Absolute out flow of a node."""
        out_flow: float = np.sum([
            self.net.edges[edge]["flow"]
            for edge in self.net.out_edges(node)
        ])
        return out_flow
