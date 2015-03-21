import ROOT


def initialize_ROOT():
    print "initializing ROOT.."
    # keep a pointer to the original TCanvas constructor
    oldinit = ROOT.TCanvas.__init__

    # define a new TCanvas class (inheriting from the original one),
    # setting the memory ownership in the constructor
    class GarbageCollectionResistentCanvas(ROOT.TCanvas):
        def __init__(self, *args):
            self.ANewDefinedVariable = 5
            oldinit(self, *args)
            ROOT.SetOwnership(self, False)

    # replace the old TCanvas class by the new one
    ROOT.TCanvas = GarbageCollectionResistentCanvas
