#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the binning for the analysis
# created on Oct 28th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gStyle
from numpy import arange, array, delete, append, concatenate, sqrt

from utils import get_root_vec


class Bins:
    def __init__(self, run=None, cut=None):

        self.Run = run
        self.Cut = cut
        self.Config = run.MainConfig

        self.NDevices = run.NPlanes
        self.NEntries = run.NEntries
        self.NCols, self.NRows = run.NPixels
        self.PX, self.PY = run.PixelSize

        # Binning
        self.BinSize = self.Config.getint('PLOTS', 'bin size')
        self.RawBinning = self.get_raw(self.BinSize)
        if self.Cut is not None:
            self.Binning = self.create()
            self.TimeBinning = self.load_time_binning()
            self.NBins = len(self.Binning) - 1
        self.NPulser = None
        self.Pulser = None

        # Pixel
        self.MaxADC = 2**8 - 1
        self.MinPH = -5000
        self.MaxPH = 100000
        self.PHBinWidth = 200
        self.MinVcal = -100
        self.MaxVcal = 1250
        self.AdcToVcal = 47.

        # Pad Pulse Height
        self.MinPadPH = -50
        self.MaxPadPH = 500
        self.PadPHBinWidth = 1

        # Miscellaneous
        self.GlobalCoods = [-.5025, .5175, -.505, .515]  # range in x and y in telescope coordinates [cm]
        self.FluxRange = [1, 40000]

        self.root_setup()

    # ----------------------------------------
    # region INIT
    def root_setup(self):
        gStyle.SetPalette(self.Config.getint('PLOTS', 'palette'))
        gStyle.SetNumberContours(self.Config.getint('PLOTS', 'contours'))

    def load_time_binning(self):
        return array([self.Run.get_time_at_event(event) for event in self.Binning], 'd')

    def create(self, evts_per_bin=None):
        evts_per_bin = self.BinSize if evts_per_bin is None else evts_per_bin
        if self.Cut is None:
            return self.get_raw(evts_per_bin)
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
        self.Binning = self.get()
        self.TimeBinning = self.load_time_binning()
        self.NBins = len(self.Binning)
        return value
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GENERAL
    def get_raw(self, evts_per_bin, start_event=0, end_event=None, time_bins=False):
        evts_per_bin = self.BinSize if evts_per_bin is None else evts_per_bin
        binning = arange(start_event, self.NEntries if end_event is None else end_event, self.NEntries / (self.NEntries / evts_per_bin), 'd')
        return [len(binning) - 1, array([self.Run.get_time_at_event(int(ev)) for ev in binning], 'd') if time_bins else binning]

    def get_raw_time(self, bin_width, rel_time=False, start_time=0, end_time=None):
        """ returns bins with fixed time width. bin_width in secons """
        end_time = self.Run.EndTime if end_time is None else end_time
        bins = arange(self.Run.StartTime + start_time, end_time, bin_width)
        bins = concatenate((bins, [end_time] if bins[-1] != end_time else []))
        return [bins.size - 1, bins - (self.Run.StartTime if rel_time else 0)]

    def get(self, bin_width=None, vs_time=False):
        return self.get_time(bin_width) if vs_time else self.get_event(bin_width)

    def get_event(self, bin_width=None):
        self.set_bin_size(bin_width)
        return [self.NBins - 1, array(self.Binning, 'd')]

    def get_time(self, bin_width=None):
        self.set_bin_size(bin_width)
        return [self.NBins - 1, self.TimeBinning]

    def get_pixel(self):
        x, y = (arange(-.5, n) for n in self.Run.NPixels)
        return [x.size - 1, x, y.size - 1, y]

    def get_pixel_x(self):
        return [self.Run.NPixels[0], arange(-.5, self.Run.NPixels[0])]

    def get_pixel_y(self):
        return [self.Run.NPixels[1], arange(-.5, self.Run.NPixels[1])]

    def get_tcal(self):
        return make_bins(0, self.Run.NSamples, 1)

    def get_native_global(self, mm=False):
        return self.get_global(res_fac=1, mm=mm)

    def get_global(self, res_fac, mm=False):
        return self.get_global_x(res_fac, mm) + self.get_global_y(res_fac, mm)

    def get_global_x(self, res_fac=None, mm=False):
        """ calculates the global telescope bins
        :param res_fac: telescope resolution in um
        :param mm: use mm as unit of return [default cm]
        :return: [nbins, bin_array] """
        res = self.PX * (1 / sqrt(12) if res_fac is None else res_fac)
        bins = arange(self.GlobalCoods[0], self.GlobalCoods[1] + res / 1000., res)
        return [bins.size, bins * (10 if mm else 1)]

    def get_global_y(self, res_fac=None, mm=False):
        res = self.PY * (1 / sqrt(12) if res_fac is None else res_fac)
        bins = arange(self.GlobalCoods[2], self.GlobalCoods[3] + res / 1000., res)
        return [bins.size, bins * (10 if mm else 1)]

    def get_pulser(self, n_pulser):
        if self.NPulser != n_pulser or self.Pulser is None:
            n = self.Run.Tree.Draw('Entry$', 'pulser', 'goff')
            pulser_events = get_root_vec(self.Run.Tree, n, dtype='u4')
            bins = concatenate([[0], pulser_events[arange(n_pulser, n, n_pulser)], [pulser_events[-1]]])
            bins = bins[:-1] if bins[-1] == bins[-2] else bins
            self.Pulser = bins
            self.NPulser = n_pulser
        return [self.Pulser.size - 1, array(self.Pulser, 'd')]
    # endregion GENERAL
    # ----------------------------------------

    # ----------------------------------------
    # region PIXEL
    def get_adc_bins(self):
        return make_bins(0, self.MaxADC, bin_width=1)

    def get_vcal_bins(self):
        return make_bins(self.MinVcal, self.MaxVcal, int(self.PHBinWidth / self.AdcToVcal))

    def get_ph(self, bin_width=None):
        return make_bins(self.MinPH, self.MaxPH, self.PHBinWidth if bin_width is None else bin_width)
    # endregion PIXEL
    # ----------------------------------------

    def get_pad_ph(self, bin_width=None):
        return make_bins(self.MinPadPH, self.MaxPadPH, self.PadPHBinWidth if bin_width is None else bin_width)


def make_bins(min_val, max_val, bin_width):
    bins = arange(min_val, max_val + bin_width / 100., bin_width, dtype='d')
    return [bins.size - 1, bins]
