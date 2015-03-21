from Helper.Initializer import initializer


class Pad2DHistConfig(object):
    @initializer
    def __init__(self,bins_x,min_x,max_x,bins_y=None,min_y=None,max_y=None):
        self.binsxy = bins_x
        if bins_y == None:
            self.bins_y = bins_x
        if min_y == None:
            self.min_y = min_x
        if max_y == None:
            self.max_y = max_x

class MeanSignalHistConfig(object):
    @initializer
    def __init__(self,bins,min_x,max_x):
        pass
