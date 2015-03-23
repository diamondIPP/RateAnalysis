from Helper.Initializer import initializer


class Pad2DHistConfig(object):
    @initializer
    def __init__(self,bins_x,bins_y=None):
        self.binsxy = bins_x
        if bins_y == None:
            self.bins_y = bins_x

class MeanSignalHistConfig(object):
    @initializer
    def __init__(self,bins,min_x,max_x):
        pass
