"""
Tracer() class.

For particles that travel around
a BodySystem following the flow, and leave a *trace*.

Francesc Font-Clos
Nov 2018
"""
import networkx as nx
import numpy as np
import pandas as pd


class Tracer():
    """
    Tracer particle.

    An object flowing around the bodysystem
    that leaves a trace.

    Attributes
    ----------
    net: nx.DiGraph()
        Must be taken from some Bodysystem.network_solved
    position: np.array(3)
        Current position of the tracer.
    time: float
        Current simulation time.
    node: string
        Current node where the tracer is.
    _flow: float
        Input flow of the current node.
    Notes
    -----
    _positions, _times and _nodes are lists holding the whole history.

    """

    def __init__(self,
                 net=None,
                 initial_time=0,
                 source_node=0,
                 flow_reversed=False
                 ):
        """
        Initialize a tracer particle.

        Parameters
        ----------

        """
        assert net is not None
        try:
            assert source_node in net.nodes
        except AssertionError:
            print(source_node)

        self.flow_reversed = flow_reversed
        if not flow_reversed:
            self.net = net
        else:
            self.net = net.reverse()
        self.time = initial_time
        self._times = [initial_time]
        self.node = source_node
        self._nodes = [source_node]
        self._reached_sink = False
        self._attached = False

    def travel(self):
        """Hop until you reach a sink."""
        # load temporal positiosn and flows dicts
        # for fast access during travel time
        self._flows_dict = nx.get_edge_attributes(self.net, "flow")
        self._pos_dict = nx.get_node_attributes(self.net, "position")
        while not (self._reached_sink or self._attached):
            self._hop()
        # fill in the 3d positions
        self._positions = [self._pos_dict[x] for x in self._nodes]
        self.position = self._positions[-1]
        # remove the dicts
        del self._pos_dict
        del self._flows_dict
        return self.node

    def _hop(self):
        """
        Make a hop.

        If only one outgoing edge, move to that node.
        If more than one, choose proportionally to flows

        """
        # find my neighbours
        neighbours = list(self.net.neighbors(self.node))
        # decide where to go next

        if not neighbours:
            self._reached_sink = True
            return None
        elif len(neighbours) == 1:
            neighbour = neighbours[0]
        else:
            flows = [f for u, v, f in self.net.out_edges(
                self.node, data="flow")]
            out_flow = np.sum(flows)
            p = np.array(flows) / out_flow
            # we cant assume a and neighbours are in the same order
            a = [v for u, v, f in self.net.out_edges(self.node, data="flow")]
            neighbour = np.random.choice(a, p=p)
        # decide if we attached to the current edge
        u, v = (self.node, neighbour)
        attach_prob = self.net[u][v]['attach_prob']
        if np.random.uniform() < attach_prob:
            self._attached = True
        # update current attributes
        self.time += self.net[self.node][neighbour]["time"]
        self.node = neighbour
        self._nodes.append(neighbour)
        self._times.append(self.time)


def compute_prob_attachment(
    prefactor_vessel: float,
    prefactor_tree: float,
    traverse_time: float,
    radius: float,
    capillary_radius: float,
    method: str="approximate"
) -> float:
    prob_enter_capillary_system = compute_prob_enter_capillary_system(
        prefactor=prefactor_vessel,
        traverse_time=traverse_time,
        radius=radius
    )
    prob_attach_capillary_tree = compute_prob_attach_capillary_tree(
        prefactor=prefactor_tree,
        initial_radius=min(0.5, radius),
        capillary_radius=capillary_radius,
        method=method
    )
    prob_attach = prob_enter_capillary_system * prob_attach_capillary_tree
    return prob_attach

def compute_prob_enter_capillary_system(
    prefactor: float,
    traverse_time: float,
    radius: float
) -> float:
    return prefactor * traverse_time / radius ** 2

# capillary computations
def compute_prob_attach_capillary_tree(
    prefactor: float,
    initial_radius: float,
    capillary_radius: float,
    method: str="exact"
) -> float:
    if method == "exact":
        return _compute_exact_capillary_attachment_prob(
            prefactor=prefactor,
            initial_radius=initial_radius,
            capillary_radius=capillary_radius,
        )
    elif method == "approximate":
        return _compute_approximate_capillary_attachment_prob(
            prefactor=prefactor,
            initial_radius=initial_radius,
            capillary_radius=capillary_radius,
        )
    else:
        raise ValueError(f"Method {method} not known. Available methods: 'exact', 'approximate'.")

        
def _compute_approximate_capillary_attachment_prob(
    prefactor: float,
    initial_radius: float,
    capillary_radius: float,
):  
    # make sure we are not using large initial radius
    assert initial_radius < 1
    # make sure final radius is less than initial radius
    assert initial_radius > capillary_radius
    # computations obtained with exp-log trick
    argument = 1 / capillary_radius ** 2 - 1 / initial_radius ** 2
    result = 1 - np.exp(-prefactor * argument)
    return result

def _compute_exact_capillary_attachment_prob(
    prefactor: float,
    initial_radius: float,
    capillary_radius: float,
):
    raise NotImplementedError
