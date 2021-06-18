#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the bin_width for the analysis
# created on Oct 28th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from helpers.draw import make_bins, choose, array, append, load_main_config, diff, save_pickle
from numpy import quantile
from src.dut import Plane
from src.sub_analysis import SubAnalysis


class Bins(SubAnalysis):

    Config = load_main_config()
    Size = Config.get_value('PLOTS', 'bin size', int)
    FluxRange = Config.get_list('PLOTS', 'flux range')

    # Pixel
    ADCBits = 2 ** 8
    PHRange = array([-5000, 100000])
    PHWidth = 200
    VcalToEl = 47.

    # Pad
    PadPHRange = [-50, 500]
    PadPHBinWidth = 1

    def __init__(self, analysis=None):

        if analysis is not None:
            super().__init__(analysis, pickle_dir='Bins')

            # Run Info
            self.NEvents = self.Run.NEvents

    # ----------------------------------------
    # region GET
    @staticmethod
    def w(bin_size):
        return Bins.get_size(bin_size)

    @staticmethod
    def get_size(bin_size):
        return choose(bin_size, Bins.Size)

    def get(self, time=True, bin_size=None, cut=None):
        return (self.get_time if time else self.get_events)(bin_size, cut)

    def get_time(self, bin_size=None, cut=None):
        i, events = self.get_events(bin_size, cut)
        return [i, self.Run.Time[events.astype('i4')] / 1000]

    def get_end_time(self, bin_size=None, cut=None):
        return self.Run.Time[int(self.get_events(bin_size, cut)[1][-1])] / 1000

    @save_pickle(suf_args='[0, 1]', print_dur=True)
    def get_events(self, bin_size=None, cut=None, _redo=False):
        events = self.get_tree_vec('Entry$', self.Cut(cut), 'i4')
        b = events[::self.get_size(bin_size)]
        bins = append(b, events[-1] if events[events > b[-1]].size > self.get_size(bin_size) / 4 else [])
        return [bins.size - 1, array(bins, 'd')]

    def get_raw(self, bin_width=None, start_event=0, end_event=None, vs_time=False, t_from_event=False):
        return self.get_raw_time(bin_width, start_event, end_event, t_from_event) if vs_time else self.get_raw_event(bin_width, start_event, end_event)

    def get_raw_event(self, bin_width=None, start_event=0, end_event=None):
        end_event, bin_width = choose(end_event, self.NEvents), Bins.get_size(bin_width)
        return Bins.make(start_event, end_event, bin_width, last=end_event - start_event % bin_width > bin_width / 4)

    def get_raw_time(self, bin_width=None, start_time=0, end_time=None, t_from_event=False):
        """ returns bins with fixed time width. bin_width in seconds """
        if t_from_event:
            i, events = self.get_raw_event(bin_width, start_time, end_time)
            return [i, self.Run.Time[events.astype('i4')]]
        else:
            start, end = self.Run.StartTime + start_time, self.Run.EndTime if end_time is None else self.Run.StartTime + end_time
            return Bins.make(start, end, bin_width, last=end - start % bin_width > bin_width / 4)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region TELESCOPE
    @staticmethod
    def get_pixel_x(bin_width=1, aspect_ratio=False):
        return Bins.make(-Plane.get_xpix(aspect_ratio) - .5, Plane.NCols + Plane.get_xpix(aspect_ratio) - .5, bin_width)

    @staticmethod
    def get_pixel_y(bin_width=1, aspect_ratio=False):
        return Bins.make(-Plane.get_ypix(aspect_ratio) - .5, Plane.NRows + Plane.get_ypix(aspect_ratio) - .5, bin_width)

    @staticmethod
    def get_pixel(bin_width=1, aspect_ratio=False):
        return Bins.get_pixel_x(bin_width, aspect_ratio) + Bins.get_pixel_y(bin_width, aspect_ratio)

    @staticmethod
    def get_global(res=None, square=False):
        return Bins.get_global_x(res) + (Bins.get_global_x(res) if square else Bins.get_global_y(res))

    @staticmethod
    def get_native_global():
        return Bins.get_global()

    @staticmethod
    def get_global_cood(mode, res=None):
        return Bins.get_global_x(res) if mode == 'x' else Bins.get_global_y(res)

    @staticmethod
    def get_global_x(res=None):
        wmax = Plane.WMax * .6  # to keep the aspect ratio
        return Bins.make(-wmax, wmax, choose(res, Plane.R0) * Plane.PX)

    @staticmethod
    def get_global_y(res=None):
        wmax = Plane.WMax * .6  # to keep the aspect ratio
        return Bins.make(-wmax, wmax, choose(res, Plane.R0) * Plane.PY)

    @staticmethod
    def get_angle(bin_size=.1, max_angle=4):
        return Bins.make(-max_angle, max_angle, bin_size)

    @staticmethod
    def get_chi2(bin_size=None, chi_max=100):
        return Bins.make(0, chi_max, choose(bin_size, .1))
    # endregion TELESCOPE
    # ----------------------------------------

    # ----------------------------------------
    # region PIXEL
    @staticmethod
    def get_adc():
        return Bins.make(0, Bins.ADCBits)

    @staticmethod
    def get_vcal():
        return Bins.make(*Bins.PHRange / Bins.VcalToEl)

    @staticmethod
    def get_electrons(bin_width=None):
        return Bins.make(*Bins.PHRange, choose(bin_width, Bins.PHWidth))

    @staticmethod
    def get_ph(vcal=False, adc=False, bin_width=None):
        return Bins.get_vcal() if vcal else Bins.get_adc() if adc else Bins.get_electrons(bin_width)
    # endregion PIXEL
    # ----------------------------------------

    # ----------------------------------------
    # region PAD
    @staticmethod
    def get_pad_ph(bin_width=None):
        return Bins.make(*Bins.PadPHRange, choose(bin_width, Bins.PadPHBinWidth))

    def get_wf(self, bin_width=1):
        return Bins.make(0, self.Run.NSamples, bin_width)
    # endregion PAD
    # ----------------------------------------

    @staticmethod
    def make(min_val, max_val=None, bin_width=1, last=False, n=None, off=0):
        return make_bins(min_val, max_val, bin_width, last, n, off)

    @staticmethod
    def make2d(x, y, bs=None, off=0):
        bs = max(diff(sorted(x))) if bs is None else bs
        return Bins.make(min(x), max(x) + bs, bs, last=True, off=off) + Bins.make(min(y), max(y) + bs, bs, last=True, off=off)

    @staticmethod
    def find_width(x):
        return Bins.freedman_diaconis(x)

    @staticmethod
    def freedman_diaconis(x):
        return 2 * (quantile(x, .75) - quantile(x, .25)) / x.size ** (1 / 3)
