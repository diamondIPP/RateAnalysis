from math import ceil, sqrt

from ROOT import gStyle
from numpy import arange, array, delete, append, concatenate

from utils import increased_range, get_root_vec


class Plots:
    def __init__(self, run=None, cut=None):

        self.Run = run
        self.Cut = cut
        self.Config = run.MainConfig

        self.NDevices = run.NPlanes
        self.NEntries = run.NEntries
        self.NCols = run.NPixels[0]
        self.NRows = run.NPixels[1]

        # Binning
        self.BinSize = self.Config.getint('PLOTS', 'bin size')
        self.RawBinning = self.get_raw_binning(self.BinSize)
        if self.Cut is not None:
            self.Binning = self.get_binning()
            self.TimeBinning = self.get_time_binning()
            self.NBins = len(self.Binning) - 1

        self.Settings = {'ph1Dmin': -5000,
                         'ph1Dmax': 100000,
                         'ph1Dbins': 105000 / 500,
                         'maxPadPh': 500,
                         'minPadPh': -50,
                         'ph1DbinsSi': 160,
                         'ph1DminSi': 0,
                         'ph1DmaxSi': 90000,
                         'nEventsAv': 20000,
                         'event_bins': max(int(ceil(float(self.NEntries) / 10)), 200),
                         'event_min': 0,
                         'event_max': self.NEntries,
                         'maxphplots': int(ceil(8 * self.NEntries / 100)),
                         'nBinsX': 52,
                         'nBinsY': 80,
                         'xmin': -3900,
                         'xmax': 3900,
                         'ymin': -4000,
                         'ymax': 4000,
                         'nBinCol': 51,
                         'nCols': 52,
                         'nRows': 80,
                         'minCol': 0,
                         'maxCol': 51,
                         'nBinRow': 79,
                         'minRow': 0,
                         'maxRow': 79,
                         'num_diff_cluster_sizes': 4,
                         'vcalBins': [int((1350 * 47) / 500), -100, 1250],
                         'globalCoods': [-.5025, .5175, -.505, .515],
                         'xpix': .015,
                         'ypix': .01}
        self.Settings['phBins'] = [self.Settings['ph1Dbins'], self.Settings['ph1Dmin'], self.Settings['ph1Dmax']]
        self.Settings['2DBins'] = [self.Settings['nCols'], - .5, self.Settings['nCols'] - .5, self.Settings['nRows'], - .5, self.Settings['nRows'] - .5]
        self.Settings['2DBinsX'] = [self.Settings['nCols'], - .5, self.Settings['nCols'] - .5] * 2
        self.Settings['2DBinsY'] = [self.Settings['nRows'], - .5, self.Settings['nRows'] - .5] * 2
        self.Settings['event_bins'] = int(ceil(float(self.NEntries) / 5000)) if self.NEntries <= 100000 else \
            int(ceil(float(self.NEntries) / 100)) if self.NEntries <= 500000 else int(ceil(float(self.NEntries) / self.Settings['nEventsAv']))
        self.Settings['deltaX'] = float(self.Settings['xmax'] - self.Settings['xmin']) / self.Settings['nBinsX']
        self.Settings['deltaY'] = float(self.Settings['ymax'] - self.Settings['ymin']) / self.Settings['nBinsY']
        self.FluxRange = [1, 40000]

    def root_setup(self):
        gStyle.SetPalette(self.Config.get('PLOTS', 'palette'))
        gStyle.SetNumberContours(self.Config.get('PLOTS', 'contours'))

    @staticmethod
    def get_arrays(lst):
        arrays = []
        for i in xrange(len(lst) / 3):
            step = (lst[3 * i + 2] - lst[3 * i + 1]) / float(lst[3 * i])
            arrays.append(arange(lst[3 * i + 1], lst[3 * i + 2] + step, step, 'd'))
        ret_val = []
        for arr in arrays:
            ret_val.append(len(arr) - 1)
            ret_val.append(arr)
        return ret_val

    # ----------------------------------------
    # region BINNING
    def get_raw_binning(self, evts_per_bin, start_event=0, end_event=None, time_bins=False):
        evts_per_bin = self.BinSize if evts_per_bin is None else evts_per_bin
        binning = arange(start_event, self.NEntries if end_event is None else end_event, self.NEntries / (self.NEntries / evts_per_bin), 'd')
        return [len(binning) - 1, array([self.Run.get_time_at_event(int(ev)) for ev in binning], 'd') if time_bins else binning]

    def get_raw_time_binning(self, bin_width, rel_time=False, start_time=0, end_time=None):
        """ returns bins with fixed time width. bin_width in secons """
        end_time = self.Run.EndTime if end_time is None else end_time
        bins = range(self.Run.StartTime + start_time, end_time, bin_width)
        bins += [end_time] if bins[-1] != end_time else []
        return [len(bins) - 1, array(bins, 'd') - (self.Run.StartTime if rel_time else 0)]

    def get_binning(self, evts_per_bin=None):
        evts_per_bin = self.BinSize if evts_per_bin is None else evts_per_bin
        if self.Cut is None:
            return self.get_raw_binning(evts_per_bin)
        jumps = filter(lambda x: x[1] > self.Cut.get_min_event(), self.Cut.Interruptions)  # filter out interruptions outside event range
        first_jump_end = jumps[0][1] if jumps else 0
        events = arange(min(self.Cut.get_min_event(), first_jump_end), self.Cut.get_max_event())
        for start, stop in jumps:  # remove beam interruptions from events
            events = delete(events, arange(start, stop + 1))
        events = events[::evts_per_bin]  # slice out every nth entry
        return append(events, self.Cut.get_max_event()) if self.Cut.get_max_event() != events[-1] else events

    def set_bin_size(self, value):
        if value is None:
            return
        self.BinSize = value
        self.Binning = self.get_binning()
        self.TimeBinning = self.get_time_binning()
        self.NBins = len(self.Binning)
        return value

    def get_time_binning(self):
        time_bins = []
        for event in self.Binning:
            time_bins.append(self.Run.get_time_at_event(event))
        return array(time_bins, 'd')

    def get_time_bins(self, evts_per_bin=None):
        self.set_bin_size(evts_per_bin)
        return [self.NBins - 1, self.TimeBinning]

    def get_event_bins(self, binning=None):
        self.set_bin_size(binning)
        return [self.NBins - 1, array(self.Binning, 'd')]

    def get_ph_bins(self, bin_width=1):
        return [int(ceil((self.Settings['maxPadPh'] - self.Settings['minPadPh']) / float(bin_width))), self.Settings['minPadPh'], self.Settings['maxPadPh']]

    def get_tcal_bins(self):
        return [int(round(sqrt(len(self.Run.TCal))))] + increased_range([min(self.Run.TCal), max(self.Run.TCal)], .1, .1)

    def get_global_bins(self, res, mode=None, arrays=False, mm=False):
        x, y = array(self.Settings['globalCoods'][:2]), array(self.Settings['globalCoods'][2:])
        size_x = self.Settings['xpix'] if mode in ['x', None] else self.Settings['ypix']
        size_y = self.Settings['ypix'] if mode in ['y', None] else self.Settings['xpix']
        x, y, size_x, size_y = [val * (10 if mm else 1) for val in [x, y, size_x, size_y]]
        bins = [int(ceil((x[1] - x[0]) / size_x * sqrt(12) / res)), x[0], x[1], int(ceil((y[1] - y[0]) / size_y * sqrt(12) / res)), y[0], y[1]]
        x_arr, y_arr = arange(x[0], x[1], size_x, 'd'), arange(y[0], y[1], size_y, 'd')
        return bins if not arrays else [len(x_arr) - 1, x_arr, len(y_arr) - 1, y_arr]

    def get_pulser_bins(self, n_pulser):
        n = self.Run.Tree.Draw('Entry$', 'pulser', 'goff')
        pulser_events = get_root_vec(self.Run.Tree, n)
        bins = concatenate([[0], pulser_events[arange(n_pulser, n, n_pulser)], [pulser_events[-1]]])
        bins = bins[:-1] if bins[-1] == bins[-2] else bins
        return bins.size - 1, bins
    # endregion BINNING
    # ----------------------------------------
