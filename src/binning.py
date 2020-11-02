#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the bin_width for the analysis
# created on Oct 28th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from json import loads
from numpy import arange, array, delete, append, concatenate, sqrt, zeros, ndarray, linspace
from helpers.utils import choose


class Bins:

    Config = None
    Size = 10000
    FluxRange = [5, 20000]

    # Pixel
    MaxADC = 2 ** 8 - 1
    MinPH = -5000
    MaxPH = 100000
    PHBinWidth = 200
    MinVcal = -100
    MaxVcal = 1250
    VcalToEl = 47.

    # Pad Pulse Height
    MinPadPH = -50
    MaxPadPH = 500
    PadPHBinWidth = 1

    # Miscellaneous
    GlobalCoods = [-.5025, .5175, -.505, .515]  # range in x and y in telescope coordinates [cm]

    def __init__(self, run=None, cut=None):

        self.Run = run
        self.Cut = cut
        Bins.Config = self.Run.MainConfig

        self.NCols, self.NRows = self.Run.NPixels
        self.PX, self.PY = self.Run.PixelSize
        Bins.Size = self.Config.getint('PLOTS', 'bin size')

        if run is not None:
            # Run Info
            self.NDevices = run.NPlanes
            self.NEvents = run.NEvents

            # Binning
            self.RawBinning = self.get_raw()
            if self.Cut is not None:
                self.Binning = self.create()
                self.TimeBinning = self.load_time_binning()
                self.NBins = len(self.Binning)

        Bins.FluxRange = loads(Bins.Config.get('PLOTS', 'flux range'))

    def __call__(self, bin_width=None):
        self.set_bin_size(bin_width)
        return Bins.Size

    # ----------------------------------------
    # region INIT

    def load_time_binning(self):
        return array([self.Run.get_time_at_event(event) for event in self.Binning], 'd')

    def create(self, evts_per_bin=None):
        evts_per_bin = choose(evts_per_bin, Bins.Size)
        if self.Cut is None:
            return self.get_raw(evts_per_bin)
        jumps = [r for r in self.Cut.get_interruptions_ranges() if r[1] > self.Cut.get_min_event()]  # filter out interruptions outside event range
        delete_events = concatenate([arange(self.Cut.get_min_event(), dtype='i8'), zeros(0, 'i8') if not len(jumps) else concatenate([arange(start, stop + 1) for start, stop in jumps])])
        events = delete(arange(self.Cut.get_max_event()), delete_events) if delete_events.size else arange(self.Cut.get_max_event())
        events = events[::evts_per_bin]  # slice out every nth entry
        return append(events, self.Cut.get_max_event()) if self.Cut.get_max_event() != events[-1] else events

    def set_bin_size(self, value):
        if value is None or value == Bins.Size:
            return
        Bins.Size = value
        self.Binning = self.create()
        self.TimeBinning = self.load_time_binning()
        self.NBins = len(self.Binning)
        return value
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GENERAL
    def get_raw(self, bin_width=None, start_event=0, end_event=None, vs_time=False, rel_time=False, t_from_event=False):
        return self.get_raw_time(bin_width, rel_time, start_event, end_event, t_from_event) if vs_time else self.get_raw_event(bin_width, start_event, end_event)

    def get_raw_event(self, bin_width=None, start_event=0, end_event=None):
        end_event, bin_width = choose(end_event, self.NEvents), choose(bin_width, Bins.Size)
        return Bins.make(start_event, end_event, bin_width, last=end_event - start_event % bin_width > bin_width / 4)

    def get_raw_time(self, bin_width, rel_time=False, start_time=0, end_time=None, t_from_event=False):
        """ returns bins with fixed time width. bin_width in seconds """
        offset = self.Run.StartTime if rel_time else 0
        if t_from_event:
            bins = array([self.Run.get_time_at_event(int(ev)) for ev in self.get_raw_event(bin_width, start_time, end_time)[1]], 'd')
            return [bins.size - 1, bins - offset]
        else:
            start, end = self.Run.StartTime + start_time - offset, choose(end_time, self.Run.EndTime) - offset
            return Bins.make(start, end, bin_width, last=end - start % bin_width > bin_width / 4)

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
        return Bins.make(0, self.Run.NSamples, 1)

    def get_native_global(self, mm=False):
        return self.get_global(res_fac=1, mm=mm)

    def get_global(self, res_fac=None, mm=False):
        return self.get_global_x(res_fac, mm) + self.get_global_y(res_fac, mm)

    def get_global_cood(self, mode, res_fac=None, mm=False):
        return self.get_global_x(res_fac, mm) if mode == 'x' else self.get_global_y(res_fac, mm)

    def get_global_x(self, res_fac=None, mm=False):
        """ calculates the global telescope bins
        :param res_fac: telescope resolution in um
        :param mm: use mm as unit of return [default cm]
        :return: [nbins, bin_array] """
        res = self.PX * (2 / sqrt(12) if res_fac is None else res_fac)
        bins = arange(self.GlobalCoods[0], self.GlobalCoods[1] + res / 1000., res)
        return [bins.size - 1, bins * (10 if mm else 1)]

    def get_global_y(self, res_fac=None, mm=False):
        res = self.PY * (2 / sqrt(12) if res_fac is None else res_fac)
        bins = arange(self.GlobalCoods[2], self.GlobalCoods[3] + res / 1000., res)
        return [bins.size - 1, bins * (10 if mm else 1)]

    @staticmethod
    def get_angle(bin_size=.1, max_angle=4):
        return Bins.make(-max_angle, max_angle, bin_size)
    # endregion GENERAL
    # ----------------------------------------

    # ----------------------------------------
    # region PIXEL
    def get_adc(self):
        return Bins.make(0, self.MaxADC, bin_width=1)

    def get_vcal(self):
        return Bins.make(self.MinVcal, self.MaxVcal, int(self.PHBinWidth / self.VcalToEl))

    def get_electrons(self, bin_width=None):
        return Bins.make(self.MinPH, self.MaxPH, self.PHBinWidth if bin_width is None else bin_width)

    def get_ph(self, vcal=False, adc=False, bin_width=None):
        return self.get_vcal() if vcal else self.get_adc() if adc else self.get_electrons(bin_width)
    # endregion PIXEL
    # ----------------------------------------

    # ----------------------------------------
    # region PAD
    def get_pad_ph(self, bin_width=None, mean_ph=None):
        return Bins.make(self.MinPadPH, self.MaxPadPH, choose(bin_width, self.PadPHBinWidth if mean_ph is None else mean_ph / 40))
    # endregion PAD
    # ----------------------------------------

    @staticmethod
    def make(min_val, max_val=None, bin_width=1, last=False, n=None):
        bins = array(min_val, 'd')
        if type(min_val) not in [ndarray, list]:
            min_val, max_val = choose(min_val, 0, decider=max_val), choose(max_val, min_val)
            bins = arange(min_val, max_val + (bin_width if last else 0), bin_width, dtype='d') if n is None else linspace(min_val, max_val, int(n) + 1, endpoint=True)
        return [bins.size - 1, bins]

    @staticmethod
    def make2d(x, y):
        return Bins.make(x[0], x[-1], x[1] - x[0], last=True) + Bins.make(y[0], y[-1], y[1] - y[0], last=True)
