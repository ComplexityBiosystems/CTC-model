"""
Find the optimal organ connection radii.

Finds radii for organ connections given a set of target
blood flow fractions

Francesc Font-Clos
Oct 2018
"""
import pandas as pd
import numpy as np
import argparse

from BodyPartsPy.optimize import find_organ_coefficients
from BodyPartsPy.optimize import modify_net
from BodyPartsPy.optimize import get_rel_flows
from BodyPartsPy.utils import remove_superflous_node

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Find optimal organ connection radii",
        add_help=False)

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        "--input_network",
        help="Path to pickled input network.",
        required=True,
        type=str)
    required.add_argument(
        "--output_network",
        help="Path to pickled output, optimized network.",
        required=True,
        type=str)
    required.add_argument(
        "--output_coefs",
        help="Path to store optimal coefficients",
        required=True,
        type=str)

    required.add_argument(
        "--boundary_conditions",
        help="Path to csv boundary conditions",
        required=True,
        type=str)
    required.add_argument(
        "--connections",
        help="Path to pickled list of artery-to-vein connections",
        required=True,
        type=str)
    required.add_argument(
        "--organs_nodes_dict",
        help="Path to pickled dictionary from organs to network nodes",
        required=True,
        type=str)
    required.add_argument(
        "--organs_flow_fraction_target",
        help="Path to pickled dictionary from organs to target blood flow fractions",
        required=True,
        type=str)
    required.add_argument(
        "--maxiter",
        help="Maximum number of iterations in minimizer",
        required=True,
        type=int
    )
    # Optional params
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        '-v', "--verbose",
        help="This is a binary switch flag.",
        action='store_true'
    )
    optional.add_argument(
        "--method",
        help="optimization method for scipy.optimize.minime",
        required=False,
        default="Powell",
        type=str
    )
    # Help added manually so that it shows up as optional arg
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help='Show this help message and exit.')

    # parse the arguments
    args = parser.parse_args()
    verbose = args.verbose

    # load data
    net = pd.read_pickle(args.input_network)
    boundary_conditions = dict(
        pd.read_csv(args.boundary_conditions, header=None).values
    )
    connections = pd.read_pickle(args.connections)
    organs_dict = pd.read_pickle(args.organs_nodes_dict)
    # key order should not be random in latest python versions
    # but keeping the list of keys anyway
    organs_list = list(organs_dict.keys())
    organs_target = pd.read_pickle(args.organs_flow_fraction_target)

    # fix double edges in connections
    unique_connections = set([
        tuple(sorted((u, v)))
        for u, v, _, _ in connections
    ])
    net.remove_edges_from([
        (u, v)
        for u, v, _, _ in connections
    ])
    net.add_edges_from(unique_connections, radius=1, length=1)
    original_net = net.copy()
    # nodes that cannot be removed from the network
    # for simplification because they are connected to
    # the other bodysystem
    cannot_remove = []
    for u, v in unique_connections:
        cannot_remove.append(u)
        cannot_remove.append(v)
    cannot_remove = set(cannot_remove)

    if verbose:
        print(f"# Original network")
        print(f"#   nodes {net.number_of_nodes()}")
        print(f"#   edges {net.number_of_edges()}")

    # remove nodes that have degree two
    nodes_to_remove = [
        node for node, k in net.degree()
        if k == 2 and node not in cannot_remove]
    for ntr in nodes_to_remove:
        remove_superflous_node(net, ntr)

    if verbose:
        print(f"# Simplified network")
        print(f"#   nodes {net.number_of_nodes()}")
        print(f"#   edges {net.number_of_edges()}")

    # find the optimal coefficients
    coefs = find_organ_coefficients(
        initial_guess=[1] * len(organs_dict),
        organs_list=organs_list,
        organs_dict=organs_dict,
        organs_target=organs_target,
        net=net,
        boundary_conditions=boundary_conditions,
        maxiter=args.maxiter,
        verbose=True,
        method=args.method
    )

    # reload original network (not simplified)
    # and modify it with found coefficients
    net_modified = modify_net(
        coefs=coefs,
        organs_list=organs_list,
        organs_dict=organs_dict,
        net=original_net
    )

    # save output
    pd.to_pickle(net_modified, args.output_network)
    np.savetxt(args.output_coefs, coefs)
