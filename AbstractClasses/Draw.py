#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TGraphErrors, gStyle
from os.path import join


class Draw:

    def __init__(self, config, prog_dir):

        self.Config = config
        self.Dir = prog_dir
        self.Objects = []
        self.ResultsDir = None

    def set_titles(self, on=None):
        gStyle.SetOptTitle(self.Config.getboolean('SAVE', 'activate_title') if on is None else on)

    def set_save_directory(self, name):
        self.ResultsDir = join(self.Dir, name)
