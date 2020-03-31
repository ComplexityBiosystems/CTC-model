import networkx as nx
from typing import Union, List, Dict
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import inv as sparse_inv


def solve_flows(
    net: nx.Graph,
    boundary_conditions: Dict[str, float],
    root: Union[List[str], str] = None,
    viscosity=1.,
):
    # root must be a list of nodes
    assert root is not None
    if isinstance(root, str):
        root = [root]
    # convert viscosity to (mmHg s) units
    _viscosity = viscosity * 1e-6 * 7.5006

    # define internal and external nodes
    internal_nodes = [node for node, k in net.degree if k > 1]
    external_nodes = [node for node, k in net.degree if k == 1]

    # index to label conversions for nodes
    nodes_idx2label = dict(enumerate(net.nodes))
    nodes_label2idx = {v: k for k, v in nodes_idx2label.items()}

    # sparse matrices to solve the linear algebra problem
    N = len(net.nodes)
    b = dok_matrix((N, 1))
    S = dok_matrix((N, N))

    # initial condition
    p = np.zeros(N)
    # set pressure of root
    for node, press in boundary_conditions.items():
        p[nodes_label2idx[node]] = press

    # fill in coefficients
    for u, v in net.edges:
        r = net[u][v]["radius"]
        length = net[u][v]["length"]
        cij = np.math.pi * (2 * r) ** 4 / (128 * _viscosity * length)
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

    edge_attrs = {}
    pressure_dict = nx.get_node_attributes(dinet, "pressure")
    for i, (u, v) in enumerate(net.edges):
        # get its properties
        pu = pressure_dict[u]
        pv = pressure_dict[v]
        r = net[u][v]["radius"]
        length = net[u][v]["length"]
        cij = np.math.pi * (2 * r) ** 4 / (128 * _viscosity * length)
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
    return dinet
