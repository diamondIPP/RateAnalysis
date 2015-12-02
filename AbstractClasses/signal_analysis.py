__author__ = 'micha'

# ==============================================
# IMPORTS
# ==============================================
from AbstractClasses.Elementary import Elementary
from ROOT import gROOT, TGraphErrors, TCanvas, TH2D, gStyle, TF1, TH1F
from newAnalysis import Analysis
from collections import OrderedDict
from array import array


# ==============================================
# MAIN CLASS
# ==============================================
class SignalAnalysis(Elementary):

    def __init__(self, analysis, channel, canvas=None, binning=5000):
        Elementary.__init__(self)

        # main
        self.analysis = analysis
        self.analysis = Analysis(392)
        self.channel = channel
        self.run = self.analysis.run.run_number
        self.tree = self.analysis.run.tree

        # stuff
        self.draw_option = 'COLZ'
        self.cut = self.analysis.GetNEventsCut(channel)
        self.binning = binning
        self.n_bins = int(self.analysis.run.endEvent) / int(self.binning)
        self.polarity = self.get_polarity()

        # names
        self.integral_names = self.get_integral_names()
        self.signal_num = self.get_signal_number('b', 3)
        self.signal_name = '{pol}*IntegralValues[{num}]'.format(pol=self.polarity, num=self.signal_num)


        # self.show_signal_histo()
        # self.make_histos()


    def get_integral_names(self):
        names = OrderedDict()
        self.tree.GetEntry(0)
        for i, name in enumerate(self.tree.IntegralNames):
            names[name] = i
        return names

    def get_polarity(self):
        self.tree.GetEntry(0)
        return self.tree.polarities[self.channel]

    def show_integral_names(self):
        for key, value in self.integral_names.iteritems():
            print str(value).zfill(3), key
        return

    def get_signal_number(self, region, integral):
        assert region in 'abcdef', 'wrong region'
        assert integral in [1, 2, 3], 'wrong integral'
        name = 'ch{ch}_signal_{reg}_PeakIntegral{int}'.format(ch=self.channel, reg=region, int=integral)
        return self.integral_names[name]

    def show_signal_histo(self):
        self.histo = TH1F('signalx', 'signal', 400, -100, 300)
        self.tree.Draw("{name}>>signalx".format(name=self.signal_name), 'n_tracks == 1')

    def get_binning(self):
        jumps = self.analysis.cut[self.channel].jump_ranges
        n_jumps = len(jumps['start'])
        bins = [0, self.analysis.GetMinEventCut()]
        ind = 0
        for start, stop in zip(jumps['start'], jumps['stop']):
            gap = stop - start
            # continue if first start and stop outside minevent
            if stop < bins[-1]:
                ind += 1
                continue
            # if there is a jump from the start
            if start < bins[-1] < stop:
                bins.append(stop)
                ind += 1
                continue
            # add bins until hit interrupt
            while bins[-1] + self.binning < start:
                bins.append(bins[-1] + self.binning)
            # two jumps shortly after one another
            if ind < n_jumps - 2:
                next_start = jumps['start'][ind + 1]
                next_stop = jumps['stop'][ind + 1]
                if bins[-1] + self.binning + gap > next_start:
                    print 'double gap', bins[-1]
                    gap2 = next_stop - next_start
                    bins.append(bins[-1] + self.binning + gap + gap2)
                else:
                    bins.append(bins[-1] + self.binning + gap)
            else:
                bins.append(bins[-1] + self.binning + gap)
            ind += 1
        return bins

    def get_time_binning(self):
        bins = self.get_binning()
        time_bins = []
        for event in bins:
            time_bins.append(self.analysis.run.GetTimeAtEvent(event))
        return time_bins

    def make_histos(self):
        # 2D Histos
        name = "signaltime_" + str(self.run)
        xbins = array('d', self.get_time_binning())
        # self.signaltime = TH2D(name, "signaltime", self.n_bins, 0, self.analysis.run.totalTime, 1000, -50, 300)
        self.signaltime = TH2D(name, "signaltime", len(xbins) - 1, xbins, 1000, -50, 300)
        self.analysis.run.tree.Draw("{name}:time>>{histo}".format(histo=name, name=self.signal_name),
                                    self.analysis.GetCut(self.channel),
                                    self.draw_option)


if __name__ == "__main__":
    z = SignalAnalysis(2, 0)