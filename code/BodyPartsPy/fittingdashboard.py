"""A widgets interface to verify fits."""
import os
from os.path import basename
from os.path import splitext
import pandas as pd
import glob
import numpy as np

from BodyPartsPy import BodyPart

import k3d

import ipywidgets as widgets


class FittingDashboard(object):
    """A class to hold a fitting board."""

    def __init__(self, splitted_dir=None, pickled_dir=None, max_objects=10):
        """Initialize a fitting board."""
        self.splitted_dir = splitted_dir
        self.pickled_dir = pickled_dir
        self.max_objects = max_objects
        # find FJ codes
        paths = glob.glob(os.path.join(splitted_dir, "*.obj"))
        FJcodes = [splitext(basename(path))[0] for path in paths]
        self.FJcodes_all = FJcodes
        self.FJcodes_subset = np.random.choice(self.FJcodes_all, replace=False,
                                               size=self.max_objects)
        # load everything
        self._load_objects()
        self._load_k3d_objects()
        self._load_widgets()

    def _load_objects(self):
        """Load the data."""
        bps = [self._load_object(FJcode=FJcode)
               for FJcode in self.FJcodes_subset]
        bps = [bp for bp in bps
               if bp.condition not in ["accepted", "not accepted"]]
        self.bps = bps
        self.FJcodes = [bp.FJcode for bp in self.bps]

    def _load_object(self, FJcode=None):
        """Load an object that might have been already fitted, or refit it."""
        n_nodes = 15
        cut = 0
        # case we have already fitted the object
        if os.path.isfile(os.path.join(self.pickled_dir, FJcode + ".p")):
            bp = pd.read_pickle(os.path.join(self.pickled_dir, FJcode + ".p"))

            # case object was already verified
            if bp.condition in ["accepted", "not accepted"]:
                return bp

            # increase nodes by 30%
            elif bp.condition == "increase nodes":
                n_nodes = int(bp.fit_params["n_nodes"] * 1.3 + 0.5)
                cut = bp.fit_params["cut"]

            # decrease nodes by 30%
            elif bp.condition == "decrease nodes":
                n_nodes = int(bp.fit_params["n_nodes"] * 0.7 + 0.5)
                cut = bp.fit_params["cut"]

            # increase cut by 1
            elif bp.condition == "increase cut":
                n_nodes = bp.fit_params["n_nodes"]
                cut = bp.fit_params["cut"] + 1
                if cut >= 2:
                    cut = 2

            # decrease cut by 1
            elif bp.condition == "decrease cut":
                n_nodes = bp.fit_params["n_nodes"]
                cut = bp.fit_params["cut"] - 1
                if cut <= 0:
                    cut = 0
            bp.fit(n_nodes=n_nodes, PCAsplit=cut)
            return bp

        bp = BodyPart(FJcode=FJcode, where=self.splitted_dir)
        bp.fit(n_nodes=n_nodes, PCAsplit=cut)
        return bp

    # VIZ FUNCTIONS
    def _get_k3d_meshes(self, bp):
        self.m1 = k3d.mesh(bp.mesh.vertices, bp.mesh.faces,
                           wireframe=True, color=int("262626", 16))
        self.m2 = k3d.mesh(bp.gmesh.vertices, bp.gmesh.faces,
                           wireframe=True, color=int("ce7500", 16))
        self.m3 = k3d.mesh(bp.pmesh.vertices, bp.pmesh.faces,
                           wireframe=True, color=int("3a5587", 16))

    def _show_current_function(self, a):
        # erase previous if necessary
        if self.m1 is not None:
            self.plot -= self.m1
        if self.m2 is not None:
            self.plot -= self.m2
        if self.m3 is not None:
            self.plot -= self.m3
        # load current
        self._get_k3d_meshes(self.bps[self.FJ_selector.value])
        # plot it
        self.plot += self.m1
        self.plot += self.m2
        self.plot += self.m3

    def _load_widgets(self):
        # the object selector
        self.FJ_selector = widgets.Select(options=[(label, i) for i, label
                                                   in enumerate(self.FJcodes)])

        # accept and show next button
        self.accept_and_next_widget = widgets.Button(
            description="Accept & Show next",
            button_style="success")
        self.accept_and_next_widget.on_click(self._accept_and_next_function)

        self.decline_and_next_widget = widgets.Button(
            description="Decline & Show next", button_style="danger")
        self.decline_and_next_widget.on_click(self._decline_and_next_function)

        self.mark_as_more_nodes_widget = widgets.Button(
            description="MORE nodes & Show next", button_style="warning")
        self.mark_as_more_nodes_widget.on_click(
            self._mark_as_more_nodes_function)

        self.mark_as_less_nodes_widget = widgets.Button(
            description="LESS nodes & Show next", button_style="warning")
        self.mark_as_less_nodes_widget.on_click(
            self._mark_as_less_nodes_function)

        self.mark_as_more_cut_widget = widgets.Button(
            description="MORE PCAsplit & Show next", button_style="warning")
        self.mark_as_more_cut_widget.on_click(self._mark_as_more_cut_function)

        self.mark_as_less_cut_widget = widgets.Button(
            description="LESS PCAsplit & Show next", button_style="warning")
        self.mark_as_less_cut_widget.on_click(self._mark_as_less_cut_function)

        self.show_current_widget = widgets.Button(
            description="Show selected", button_style="info")
        self.show_current_widget.on_click(self._show_current_function)

        self.restart_widget = widgets.Button(
            description="Save & Restart", button_style="danger")
        self.restart_widget.on_click(self._restart_function)

        self.buttons = [
            self.accept_and_next_widget,
            self.decline_and_next_widget,
            self.show_current_widget,
            self.mark_as_less_cut_widget,
            self.mark_as_more_cut_widget,
            self.mark_as_less_nodes_widget,
            self.mark_as_more_nodes_widget,
            self.restart_widget
        ]

    def _load_k3d_objects(self):
        # text objects
        self.finish_text = k3d.text2d("FINISHED",
                                      position=(.5, .5),
                                      size=3, )
        self.processing_text = k3d.text2d("PROCESSING",
                                          position=(.5, .5),
                                          size=3, )

        # the plot area
        self.plot = k3d.plot()

        # gosht objects
        self.m1 = None
        self.m2 = None
        self.m3 = None

    # WIDGET FUNCTIONS
    def _accept_and_next_function(self, a):
        # mark current as accepted
        self.bps[self.FJ_selector.value].mark_as("accepted")
        self._show_next()

    def _decline_and_next_function(self, a):
        # mark current as accepted
        self.bps[self.FJ_selector.value].mark_as("not accepted")
        self._show_next()

    def _mark_as_less_cut_function(self, a):
        # mark current as accepted
        self.bps[self.FJ_selector.value].mark_as("decrease cut")
        self._show_next()

    def _mark_as_more_cut_function(self, a):
        # mark current as accepted
        self.bps[self.FJ_selector.value].mark_as("increase cut")
        self._show_next()

    def _mark_as_less_nodes_function(self, a):
        # mark current as accepted
        self.bps[self.FJ_selector.value].mark_as("decrease nodes")
        self._show_next()

    def _mark_as_more_nodes_function(self, a):
        # mark current as accepted
        self.bps[self.FJ_selector.value].mark_as("increase nodes")
        self._show_next()

    def _show_next(self):
        # case we have not finished
        if self.FJ_selector.value < len(self.FJcodes) - 1:
            self.FJ_selector.value += 1
            self._show_current_function(None)
        else:
            self.plot += self.finish_text

    def _restart_function(self, a):
        self.plot -= self.finish_text
        self.plot -= self.m1
        self.plot -= self.m2
        self.plot -= self.m3

        for bp in self.bps:
            bp.to_pickle(where=self.pickled_dir)

        self.plot += self.processing_text
        self._load_objects()
        self.plot -= self.processing_text
        self.FJ_selector.value = 0
        self.FJ_selector.options = [(label, i)
                                    for i, label in enumerate(self.FJcodes)]
