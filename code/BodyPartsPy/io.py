"""
Input/Output.

Functions to hangle input/output operations
such as reading .obj files, etc

"""
import trimesh
import os
import glob
from numpy.random import choice
from os.path import basename


ISA_FJ_PATH = "../data/BodyParts3D/isa_BP3D_4.0_obj_99/"
PARTOF_FJ_PATH = "../data/BodyParts3D/partof_BP3D_4.0_obj_99/"


def load_mesh(FJcode, where="isa", use_embree=True):
    """
    Load a mesh by FJ code.

    By default it searches the 'partof' data dir.
    It also knows about the 'isa' data dir an can
    search in custom locations
    """
    if where == "partof":
        return _load_partof_mesh(FJcode)
    elif where == "isa":
        return _load_isa_mesh(FJcode)
    else:
        path = os.path.join(where, FJcode) + ".obj"
        mesh = trimesh.load(path, use_embree=use_embree)
        return mesh


def get_random_FJcode(where="isa"):
    """
    Get a random FJ code.

    Give one random FJcode from the .obj files
    found in 'where'.

    Parameters
    ----------
    where :  str
        Either "isa", "partof" or a custom directory.

    Returns
    -------
    A random FJ code.

    """
    if where == "isa":
        where = ISA_FJ_PATH
    elif where == "partof":
        where = PARTOF_FJ_PATH
    paths = glob.glob(os.path.join(where, "*.obj"))

    return choice([basename(x).replace(".obj", "") for x in paths])


def _load_partof_mesh(FJcode):
    path = os.path.join(PARTOF_FJ_PATH, FJcode) + ".obj"
    mesh = trimesh.load(path)
    return mesh


def _load_isa_mesh(FJcode):
    path = os.path.join(ISA_FJ_PATH, FJcode) + ".obj"
    mesh = trimesh.load(path)
    return mesh
