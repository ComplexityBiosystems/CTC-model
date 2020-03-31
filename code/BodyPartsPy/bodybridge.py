"""Class to hod connections between bodyparts."""

import numpy as np
from scipy.spatial.distance import norm


class BodyBridge(object):
    """
    A bridge is an edge connecting to bodyparts.

    We put bridges connecting nodes close to intersections
    with bodyparts.
    """

    def __init__(self, bp0, bp1):
        """
        Create a bridge from bp0 to bp1.

        blah
        """
        self.verified = False
        self.bp0 = bp0
        self.bp1 = bp1
        # find intersection
        self.find_intersection_naive()
        # find the candidates
        self.find_candidates()
        # initial condition
        self.mark_as("not verified")

    def set_center_no_inter(self):
        """
        Set the center of the bridge.

        To be used only when there is not intersection between the bodyparts.

        """
        m0 = self.bp0.mesh.vertices
        m1 = self.bp1.mesh.vertices
        a = np.array([norm(m1 - v, axis=1) for v in m0])
        i0, i1 = np.unravel_index(np.argmin(a, axis=None), a.shape)
        self.center = 0.5 * (m0[i0] + m1[i1])

    def find_intersection_naive(self):
        """
        Find the intersection between the two bodyparts.

        This funciton does NOT depend on blender. It naively
        finds the set of points of bp0 in bp1 and
        the points of bp1 in bp0.
        """
        vertices0 = self.bp0.mesh.vertices
        vertices1 = self.bp1.mesh.vertices
        points0 = vertices0[self.bp1.mesh.contains(vertices0)]
        points1 = vertices1[self.bp0.mesh.contains(vertices1)]
        points = np.concatenate([points0, points1])
        if len(points) > 0:
            self.center = np.mean(points, axis=0)
        else:
            self.set_center_no_inter()

    def find_candidates(self):
        """Find candidate nodes to be joined."""
        cand0 = np.argmin(norm(self.bp0.nodes - self.center, axis=1))
        cand1 = np.argmin(norm(self.bp1.nodes - self.center, axis=1))
        self.candidate0 = cand0
        self.candidate1 = cand1

    def mark_as(self, condition):
        """
        Mark the object according to verifycation degree.

        Parameters
        ----------
        condition : str
            In which condition the BodyBridge object is in.

            If the bridge has not been inspected,
            use "not verified".

            If the bridge has been inspected an is OK,
            use "accepted"

            If the bridge has been insepcted and is not OK,
            use "not accepted"


        """
        assert condition in ["accepted", "not accepted",
                             "not verified"
                             ]
        self.condition = condition

    def _find_intersection_blender(self):
        """
        Find the intersection between the two bodyparts.

        This funciton uses the trimesh.base.Trimesh.intersection function
        which depends on blender.
        """
        # find intersection
        intersection = self.bp0.mesh.intersection(self.bp1.mesh)
        center = intersection.center_mass
        self.center = center
