"""
Functions to find coefficients that match given flow ratios.

Francesc Font Clos
April 2019
"""
from .flows import solve_flows
import numpy as np
from scipy.optimize import minimize
import networkx as nx

from typing import List


def find_organ_coefficients(
    initial_guess: np.ndarray,
    organs_list: List,
    organs_dict: dict,
    organs_target: dict,
    net: nx.Graph,
    boundary_conditions: dict,
    maxiter: int,
    verbose: bool,
    tol=1e-6,
    method="Powell"
) -> np.ndarray:
    # run minimizer
    fixed_args = organs_list, organs_dict, organs_target, net, boundary_conditions
    res = minimize(
        fun=opt_fun,
        x0=initial_guess,
        args=fixed_args,
        callback=lambda x: print(x),
        tol=tol,
        options={"maxiter": maxiter, "disp": verbose},
        method=method
    )

    if verbose:
        if not res.success:
            print("# [ WARNING ] no success in optimization")
        else:
            print("# [ OK ] optimization successful")

    # get modified network
    coefs = res.x
    return coefs


def get_rel_flows(organs_dict, dinet, boundary_conditions):
    # get organs flow
    organs_flow = {k: 0 for k in organs_dict.keys()}
    for organ, edges in organs_dict.items():
        for edge in edges:
            try:
                organs_flow[organ] += dinet[edge[0]][edge[1]]["flow"]
            except:
                organs_flow[organ] += dinet[edge[1]][edge[0]]["flow"]

    # get total flow
    roots = [k for k, v in boundary_conditions.items() if v == 100]
    assert len(roots) == 1
    root = roots[0]
    _, _, total_flow = list(dinet.out_edges(root, data="flow"))[0]

    # get relative organs flow
    rel_organs_flow = {
        organ: organs_flow[organ] / total_flow
        for organ in organs_dict
    }
    return rel_organs_flow


def get_error(dinet, organs_dict, organs_target, boundary_conditions):
    # make sure target makes sense
    assert sum(organs_target.values()) <= 1
    # get relative flows
    rel_organs_flow = get_rel_flows(organs_dict, dinet, boundary_conditions)
    # get total squarred error
    error = np.sum([
        (organs_target[organ] - rel_organs_flow[organ]) ** 2
        for organ in organs_dict.keys()
    ])
    return error


def modify_net(coefs, organs_list, organs_dict, net):
    nnet = net.copy()
    for coef, organ in zip(coefs, organs_list):
        new_radii = {e: coef**2 * net[e[0]][e[1]]
                     ["radius"] for e in organs_dict[organ]}
        nx.set_edge_attributes(nnet, new_radii, name="radius")
    return nnet


def opt_fun(
    coefs,  # parameters being optimized
    organs_list, organs_dict, organs_target, net, boundary_conditions  # fixed parameters
):
    nnet = modify_net(coefs, organs_list, organs_dict, net)
    dinet = solve_flows(
        net=nnet,
        boundary_conditions=boundary_conditions,
        root=list(boundary_conditions.keys())
    )
    error = get_error(dinet, organs_dict, organs_target, boundary_conditions)
    return error
