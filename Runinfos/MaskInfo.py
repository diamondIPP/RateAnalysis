#!/usr/bin/env python

"""
Class MaskInfo for 2014 September PSI Testbeam Analysis.
"""


# ##############################
# Imports
###############################

import json
import os.path

from Helper.Initializer import initializer
from DataTypes import data_types


###############################
# MaskInfo
###############################

class MaskInfo:
    # Dictionary of all available masks
    #  -newly created objects register automatically
    #  -keys are of the format "diamond-data_type-mask_time", ie:
    #    S30-DATA-2030
    masks = {}


    # Helper function: create_name
    # Create name (used as key for masks dict) from diamond, data_type and mask_time
    @staticmethod
    def create_name(diamond,
                    bias_sign,
                    data_type,
                    mask_time):
        if data_type != 0 and data_type != 1:
            data_type = 0
        key = "{0}-{1}-{2}-{3:04d}".format(diamond, bias_sign, data_types[data_type], mask_time)
        return key

    # End of create_name


    # initializer - a member variable is automatically created
    #  for each argument of the constructor
    @initializer
    def __init__(self,
                 data_type,  # [int] (key from data_types dictionary)
                 bias_sign,
                 diamond,  # [string] (example: "S30")
                 mask_time,  # [int] (example: 2030)
                 min_row_roc0,  # [int]
                 max_row_roc0,  # [int]
                 min_col_roc0,  # [int]
                 max_col_roc0,  # [int]
                 min_row_roc4,  # [int]
                 max_row_roc4,  # [int]
                 min_col_roc4,  # [int]
                 max_col_roc4,  # [int]
                 min_x=-2,  # [float]
                 max_x=-2,  # [float]
                 min_y=-2,  # [float]
                 max_y=-2,
                 areas=[]):  # [float]

        # Add to masks dictionary
        name = self.create_name(self.diamond,
                                self.bias_sign,
                                self.data_type,
                                self.mask_time)
        MaskInfo.masks[name] = self

    # End of __init__

    # Dump all MaskInfos (the content of the masks dictionary)
    #  to a file using json
    @classmethod
    def dump(cls, filename):
        print 'save MaskInfo:',filename
        f = open(filename, "w")
        f.write(json.dumps(cls.masks,
                           default=lambda o: o.__dict__,
                           sort_keys=True,
                           indent=4))
        f.close()

    # End of to_JSON


    # Read all MaskInfos from a file and use to intialize objects
    @classmethod
    def load(cls, filename):
        # first get the dictionary from the file..
        f = open(filename, "r")
        data = json.load(f)
        f.close()
        # ..then intialize the individual MaskInfo objects from it
        for k, v in data.iteritems():
            # print k, v
            MaskInfo(**v)
        # print '\tLoad: ',k, v

        # End of to_JSON


# End of class MaskInfo  


# Run a simple test when called from command line
if __name__ == "__main__":
    fname = "masks.json"
    if os.path.isfile(fname):
        print 'Load MaskInfo: ', fname
        MaskInfo.load(fname)

    # MaskInfo(1, "IIa-2", 2518, -1, 65, 79, 1, 50, 1, 50, 1, 50, -0.2, 0.2, 0., 0.4)
    MaskInfo.dump(fname)
    print MaskInfo.masks
