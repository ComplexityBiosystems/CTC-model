"""
Find artery/vein contact points for one organ.

Finds all network nodes that are *inside* a given organ.
Should be used with *solved* bodysystems as otherwise
there are no guarantees the contact points will get a flow later.

Francesc Font-Clos
Oct 2018

"""
import argparse
import os
from typing import List

import numpy as np
import pandas as pd
import networkx as nx
from collections import OrderedDict
import trimesh


def find_contact_points(
    dinet: nx.DiGraph,
    organ_path: str,
    contact_threshold: float = 5,
    allowed_nodes: List[str] = None,
    chunk_size: int = 500
) -> List[str]:
    """Find which nodes of a solved flows network are making contanct with an organ.

    Parameters      
    ----------
    dinet : networkx.DiGraph
        Directed network of solved flows.
    organ_path : str
        Path to organ's .obj file
    contact_threshold : float, optional
        Tolerance in mm to consider there is contact. (the default is 5)
    allowed_nodes : List[str]
        List of nodes allowed to be considered contact points
    chunk_size : int
        Iterate over points in chunks of size chunk_size.

    Returns
    -------
    List[str]
        List of contact points (node labels).
    """

    node_pos_dict = OrderedDict(
        nx.get_node_attributes(dinet, "position"))
    if allowed_nodes is None:
        nodes = np.array(list(node_pos_dict.keys()))
        positions = np.array(list(node_pos_dict.values()))
    else:
        nodes = []
        positions = []
        for node, pos in node_pos_dict.items():
            if node in allowed_nodes:
                nodes.append(node)
                positions.append(pos)
        nodes = np.array(nodes)
        positions = np.array(positions)

    organ = trimesh.load(organ_path)

    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    is_inside = []
    for pos in chunker(positions, chunk_size):
        print(f"organ: {organ_path}, chunk_size: {len(pos)}")
        is_inside.append(
            organ.nearest.signed_distance(pos) > -1 * contact_threshold
        )
    is_inside = np.concatenate(is_inside)
    return nodes[is_inside]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Find all artery/vein-organ contact points.",
        add_help=False)

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        "--solved_flows_network",
        help="Path to netowrkx solved flows net (pickled)",
        default=".",
        required=True,
        type=str)
    required.add_argument(
        "-g", "--organ",
        help="Path to organ's .obj file",
        default=".",
        required=True,
        type=str)
    required.add_argument(
        "-o", "--output_dir",
        help="Directory where to save output."
             " Filename is generated from input names.",
        default=".",
        required=True,
        type=str)
    required.add_argument(
        "-t", "--tolerance",
        help="Tolerance in mm to consider there is contact",
        required=True,
        type=float)

    # Optional params
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "--pickled_connections",
        help="Restrict contact points to AV connections",
        default=None,
        required=False,
        type=str)
    optional.add_argument(
        "-v", "--verbose",
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
    verbose = args.verbose
    organ_path = args.organ
    organ_name, file_ext_obj = os.path.splitext(os.path.basename(organ_path))
    network_path = args.solved_flows_network
    network_name, file_ext_p = os.path.splitext(
        os.path.basename(network_path))
    threshold = args.tolerance

    # make sure file formats and names make sense
    assert file_ext_obj == ".obj"
    assert file_ext_p == ".p"

    # load the network
    dinet = pd.read_pickle(network_path)

    allowed_nodes = args.pickled_connections
    if allowed_nodes is not None:
        connections = pd.read_pickle(allowed_nodes)
        tmp1 = [
            a
            for a, _, _, _ in connections
        ]
        tmp2 = [
            b
            for _, b, _, _ in connections
        ]
        allowed_nodes = tmp1 + tmp2

    # find nodes of bodysystem that are inside the organ
    nodes_inside = find_contact_points(
        dinet=dinet,
        organ_path=organ_path,
        contact_threshold=threshold,
        allowed_nodes=allowed_nodes
    )

    # save output
    out_name = "_".join(
        [organ_name, network_name, f"threshold{threshold}", "contacts"])
    out_file = os.path.join(args.output_dir, out_name) + ".txt"
    np.savetxt(out_file, nodes_inside, fmt="%s")

    # print info
    if verbose:
        print(f"# Input bodysystem: {network_path}")
        print(f"# Input organ: {organ_path}")
        print(f"#")
        print(f"# Number of contact points found: {len(nodes_inside)}")
        print(f"#")
        print(f"# Output file: {out_file}")
