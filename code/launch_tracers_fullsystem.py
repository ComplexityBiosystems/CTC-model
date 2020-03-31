"""
Launch a large number of tracers from a given organ into the blood stream.

This script uses a bodysystem object filled with edge attributes to determine
attachment probabilities, and a contact_points file to define possible starting
locations.

Francesc Font-Clos
Jan 2019

"""
import argparse
import os
from pathlib import Path
import pandas as pd
import networkx as nx
import glob
from BodyPartsPy.tracer import Tracer
import numpy as np
from tqdm import tqdm


def launch_tracers(
    net: nx.DiGraph,
    organ_name: str,
    path_to_contacts: str,
    num_tracers: int,
    verbose: bool = False
):
    """Launche many tracers for a given organ.

    Parameters
    ----------
    path_to_artery_contacts : str
        Path to list of artery contant points with organ
    path_to_vein_contacts : str
        Path to list of vein contant points with organ
    num_tracers : int
        Number of tracers to launch. They will be splitted
        among available contact points.

    """
    data = []
    # arteries network
    for system_name, path_to_contacts in [
        ("main", path_to_contacts),
    ]:

        # arteries contact points
        contact_points = np.loadtxt(path_to_contacts, dtype=str)
        if len(contact_points) > 0:
            for _ in tqdm(range(num_tracers), desc=f"Launching tracers ({organ_name}, {system_name}): ", ascii=True):
                source_node = np.random.choice(contact_points)
                # create tracer
                tracer = Tracer(
                    net=net,
                    source_node=source_node
                )
                final_node = tracer.travel()
                data.append([
                    organ_name,
                    system_name,
                    source_node,
                    tracer._reached_sink,
                    tracer._attached,
                    final_node,
                    tracer._positions,
                    tracer._nodes,
                    tracer._times
                ])
                if verbose:
                    print(f"""
 # system: {system_name}
 # source_node: {source_node}
 # _reached_sink: {tracer._reached_sink}                
 # _attached: {tracer._attached}
 # final_node: {final_node}
 # final_pos: {tracer.position}
 
                """)
    df = pd.DataFrame(
        data=data,
        columns=[
            "organ_name",
            "system_name",
            "source_node",
            "_reached_sink",
            "_attached",
            "final_node",
            "_positions",
            "_nodes",
            "_times"
        ]
    )
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Launch many tracers from an organ into the blood stream",
        add_help=False)

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        "--organ",
        help="Name of recognized organ. See Makefile for list of available organs",
        required=True,
        type=str
    )
    required.add_argument(
        "--network",
        help="Full body network with attachment probabilities",
        required=True,
        type=str
    )

    required.add_argument(
        "-n", "--num_tracers",
        help="Number of tracers (to be splitted among available starting points",
        required=True,
        type=int
    )
    required.add_argument(
        "--output_dir",
        help="Path to output dir",
        required=True,
        type=str
    )

    # Optional params
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        "--verbose",
        help="This is a binary switch flag.",
        action='store_true'
    )
    optional.add_argument(
        "--drop_trajectories",
        help="do not store full trajectories, just final position",
        action='store_true'
    )

    optional.add_argument(
        "--keep_times",
        help="additionaly store simulation time",
        action='store_true'
    )

    optional.add_argument(
        "--map_endpoints",
        help="map the endpoints to organs",
        action='store_true'
    )

    optional.add_argument(
        "--suffix",
        help="Suffix to add to output files",
        required=False,
        default="",
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
    organ = args.organ
    num_tracers = args.num_tracers
    verbose = args.verbose
    output_dir = args.output_dir
    path_to_network = args.network
    drop_trajectories = args.drop_trajectories
    keep_times = args.keep_times
    map_endpoints = args.map_endpoints

    # check that all needed files exist
    contacts_dir = Path(os.path.join(
        os.path.split(__file__)[0], "../output/data/organs_contact_points/"))
    path_to_contacts = contacts_dir / \
        f"{organ}_solved_main_network_threshold10.0_contacts.txt"
    assert os.path.isdir(contacts_dir)
    assert os.path.isfile(path_to_contacts)
    assert os.path.isdir(output_dir)
    assert os.path.isfile(path_to_network)
    # MAIN
    # generate output filename
    output_file_name_full = \
        f"{organ}_tracers_{num_tracers}{args.suffix}.p"
    path_to_output_full = Path(output_dir) / output_file_name_full
    # make sure we're not overwritting
    assert not os.path.isfile(path_to_output_full)
    # read the network data
    net = pd.read_pickle(path_to_network)
    # launch the tracers
    df = launch_tracers(
        net=net,
        organ_name=organ,
        path_to_contacts=path_to_contacts,
        num_tracers=num_tracers,
        verbose=verbose
    )

    # POSTPROCESSING
    # reduce file size dropping trajectories
    if drop_trajectories:
        df.drop(["_positions", "_nodes"], axis=1, inplace=True)
    # reduce file size dropping times
    if not keep_times:
        df.drop("_times", axis=1, inplace=True)
    # map endpoints to organs
    if map_endpoints:
        contact_points_dict = {
            (fpath.split("/")[-1].split("_solved_")[0]
             ): set(pd.read_csv(fpath, header=None).values.T[0])
            for fpath in glob.glob("../output/data/organs_contact_points/*_solved_main_network_threshold10.0_contacts.txt")
        }
        _x = {}
        for organ, nodes in contact_points_dict.items():
            for node in nodes:
                _x[node] = organ
        df["target_organ"] = df.final_node.apply(
            lambda x: _x[x] if x in _x else "skin"
        )
    df.to_pickle(path_to_output_full)
