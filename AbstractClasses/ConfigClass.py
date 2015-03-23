#from Helper.Initializer import initializer


class Pad2DHistConfig(object):
    #@initializer
    def __init__(self,bins_x,bins_y=None):
        self.bins_x = bins_x
        self.binsxy = bins_x
        if bins_y == None:
            self.bins_y = bins_x
        else:
            self.bins_y = bins_y

class MeanSignalHistConfig(object):
    #@initializer
    def __init__(self,bins,min_x,max_x):
        pass
