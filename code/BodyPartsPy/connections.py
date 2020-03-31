"""
Mesh connections.

Functions to stitch together graphs obtained
with the mesh2graph module.

"""
from .io import load_mesh
from .io import PARTOF_FJ_PATH
from .io import ISA_FJ_PATH

import trimesh
import os
import glob


def get_touching_objs(FJcode, where="partof"):
    """
    Get list of objects that touch a given object.

    Parameters
    ----------
    FJcode : str
        Object being tested, must be in .../where/FJcode.obj
    where : str
        Path to dir containing objects to test contact.

    """
    objs = []
    mesh0 = load_mesh(FJcode=FJcode, where=where)

    # where to look for rest of objects
    if where == "partof":
        path_to_objs_dir = PARTOF_FJ_PATH
    elif where == "isa":
        path_to_objs_dir = ISA_FJ_PATH
    else:
        path_to_objs_dir = where

    # not using io.load_mesh because of glob pattern
    for path in glob.glob(os.path.join(path_to_objs_dir, "*.obj")):
        mesh1 = trimesh.load(path)
        shared = mesh0.contains(mesh1.vertices)
        if any(shared):
            other = os.path.basename(path).replace(".obj", "")
            objs.append(other)
    objs.remove(FJcode)
    return objs
