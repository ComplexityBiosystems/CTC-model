"""
Solve the flows problem in a network:

Find pressures p_i for all vertices and flows s(i, j) for all edges, such that:
+ input flow equals output flow for all non-leaf nodes.
+ pressure at leaf nodes is given as boundary conditions


Francesc Font-Clos
Feb 2019
"""

import pandas as pd
import argparse
from BodyPartsPy.flows import solve_flows
import os
import pathlib

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Solve the flows problem in a network",
        add_help=False)

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        "--network_pickled",
        help="The network to be solved (pickled nx.Graph)",
        required=True,
        type=str)
    required.add_argument(
        "--boundary_conditions_csv",
        help="Pressure of the network's leaf nodes (csv: node_label,pressure)",
        required=True,
        type=str)
    required.add_argument(
        "--output_file",
        help="Path to output file (pickled nx.Digraph)",
        required=True,
        type=str)

    # Optional params
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        '-v', "--verbose",
        help="This is a binary switch flag.",
        action='store_true'
    )

    # Help added manually so that it shows up as optional arg
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help='Show this help message and exit.')

    # parse the arguments
    args = parser.parse_args()
    network_pickled = pathlib.Path(args.network_pickled)
    boundary_conditions_csv = pathlib.Path(args.boundary_conditions_csv)
    verbose = args.verbose
    output_file = pathlib.Path(args.output_file)

    # basic checks
    assert network_pickled.is_file()
    assert boundary_conditions_csv.is_file()

    # load data
    network = pd.read_pickle(network_pickled)
    boundary_conditions = dict(pd.read_csv(
        boundary_conditions_csv, header=None).values)

    solved_network = solve_flows(
        network,
        root=list(boundary_conditions.keys()),
        boundary_conditions=boundary_conditions
    )

    pd.to_pickle(solved_network, output_file)
