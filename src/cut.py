# --------------------------------------------------------
#       general class to handle all the cut strings for the analysis
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from __future__ import print_function
from draw import *
from json import loads
from ROOT import TCut, TPie, TProfile
from InfoLegend import InfoLegend
from binning import Bins
from ConfigParser import NoOptionError
from numpy import histogram, histogram2d, split, where, linspace


class Cut:
    """ Contains methods to generate the cut strings for the TelescopeAnalysis and holds the dictionaries for the settings and all cut strings. """

    def __init__(self, parent):

        self.Analysis = parent
        self.Run = self.Analysis.Run
        self.TCString = self.Analysis.TCString
        self.InfoLegend = InfoLegend(parent)
        self.Bins = Bins(self.Run)

        # Configuration
        self.Config = self.Analysis.Config
        self.CutConfig = self.load_config()
        self.LowRateRun = None
        self.HighRateRun = None

        # Cut Strings
        self.CutStrings = CutStrings()

        # generate cut strings
        self.generate()

    def __call__(self, cut=None):
        return self.CutStrings() if cut is None else TCut(cut)

    def has(self, name):
        return bool(self.get(name).GetTitle())

    # ----------------------------------------
    # region CONFIG
    def load_config(self):
        return {'jump_range': loads(self.Config.get('CUT', 'exclude around jump')),
                'event_range': self.load_event_range(loads(self.Config.get('CUT', 'event range'))),
                'chi2_x': self.Config.getint('CUT', 'chi2X'),
                'chi2_y': self.Config.getint('CUT', 'chi2Y'),
                'track angle': self.Config.getint('CUT', 'track angle')}

    def update_config(self):
        pass

    def load_event_range(self, event_range=None):
        """ Generates the event range. Negative values are interpreted as minutes. Example: [-10, 700k] => 10 min < events < 700k. """
        event_range = [0, 0] if event_range is None else [self.Analysis.get_event_at_time(seconds=abs(value * 60)) if value < 0 else value for value in event_range]
        return [event_range[0], self.Run.NEntries if not event_range[1] else event_range[1]]

    def set_config(self, key, value):
        self.CutConfig[key] = value

    def get_config(self, option, section='CUT', typ=None):
        return None if not self.Config.has_option(section, option) else self.Config.get(section, option) if typ is None else typ(self.Config.get(section, option))

    def set_event_range(self, event_range):
        self.set_config('event_range', self.load_event_range(event_range))

    def load_fiducial(self, name='fiducial', warn=True):
        splits = (loads(self.Config.get('SPLIT', 'fiducial')) if self.Config.has_option('SPLIT', 'fiducial') else []) + [int(1e10)]
        n = next(i + 1 for i in xrange(len(splits)) if self.Run.Number <= splits[i])
        option = name if self.Config.has_option('CUT', name) and n == 1 else '{} {}'.format(name, n)
        return array(self.load_dut_config(option, warn=warn)) if self.load_dut_config(option, warn=warn) is not None else None

    def load_dut_config(self, option, store_true=False, warn=True):
        try:
            conf = loads(self.Config.get('CUT', option))
            dia = self.Analysis.DUT.Name
            return dia in conf if store_true else conf[dia]
        except (KeyError, NoOptionError):
            if warn:
                log_warning('No option {} in the analysis config for {}!'.format(option, make_tc_str(self.TCString)))
    # endregion CONFIG
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, name):
        return self.CutStrings.get(name)

    def get_event_range(self):
        """ :return: event range [fist event, last event], type [ndarray] """
        return array(self.CutConfig['event_range'])

    def get_min_event(self):
        """ :return: number of the first event, type [int] """
        return self.get_event_range()[0]

    def get_max_event(self):
        """ :return: number of the last event, type [int] """
        return self.get_event_range()[1]

    def get_n_events(self):
        """ :return: number of events in event range, type [int] """
        return self.get_max_event() - self.get_min_event()

    @staticmethod
    def get_track_var(num, mode, mm=False):
        return 'dia_track_{m}_local[{n}]{s}'.format(m=mode, n=num, s='*10' if mm else '')

    def get_track_vars(self, num, mm=False):
        return (self.get_track_var(num, v, mm) for v in ['y', 'x'])

    def get_beam_interruptions(self):
        """ :returns: list of raw interruptions, type [list[tup]]"""
        pickle_path = self.Analysis.make_pickle_path('BeamInterruptions', run=self.Run.Number, suf='_'.join(str(i) for i in self.CutConfig['jump_range']))
        return do_pickle(pickle_path, self.find_beam_interruptions)

    def get_interruptions_ranges(self):
        """ :returns: list of interruptions including safety margin from the AnalysisConfig. """
        range_pickle = self.Analysis.make_pickle_path('BeamInterruptions', 'Ranges', run=self.Run.Number, suf='_'.join(str(i) for i in self.CutConfig['jump_range']))
        return do_pickle(range_pickle, self.create_interruption_ranges, interruptions=self.get_beam_interruptions())

    def get_fiducial_size(self):
        xy = self.CutConfig['fiducial'] * 10  # in mm
        dx, dy = xy[1] - xy[0], xy[3] - xy[2]
        return concatenate([xy, [dx, dy, dx * dy]])

    def get_fiducial_area(self):
        return None if self.CutConfig['fiducial'] is None else self.get_fiducial_size()[-1] / 100.
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SET
    def set(self, name, value):
        self.CutStrings.set(name, value)

    def set_high_low_rate_run(self, high_run, low_run):
        self.LowRateRun = str(low_run)
        self.HighRateRun = str(high_run)

    def reset(self, name):
        self.CutStrings.reset(name)

    def update(self, name, value=None):
        self.CutStrings.set(name, value)

    def set_chi2(self, value):
        self.CutConfig['chi2_x'] = value
        self.CutConfig['chi2_y'] = value
        self.update('chi2_x', self.generate_chi2('x').Value)
        self.update('chi2_y', self.generate_chi2('y').Value)

    def reload(self):
        self.CutConfig = self.load_config()
        self.update_config()
        self.CutStrings.reset_all()
        self.generate()
        self.generate_dut()
    # endregion SET
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE
    def generate(self):
        """ Creates all cut strings. """

        # -- EVENT RANGE --
        self.CutStrings.register(self.generate_event_range(), level=10)
        self.CutStrings.register(self.generate_beam_interruptions(), 11)

        # -- EVENT ALIGNMENT --
        self.CutStrings.register(self.generate_aligned(), 12)

        # --TRACKS --
        self.CutStrings.register(self.generate_tracks(), 22)
        self.CutStrings.register(self.generate_chi2('x'), 72)
        self.CutStrings.register(self.generate_chi2('y'), 73)
        self.CutStrings.register(self.generate_track_angle('x'), 74)
        self.CutStrings.register(self.generate_track_angle('y'), 75)

    def generate_dut(self):
        pass

    @staticmethod
    def generate_tracks():
        return CutString('tracks', 'n_tracks == 1', 'only 1 track per event')

    def generate_event_range(self, min_event=None, max_event=None):
        event_range = [cfg if arg is None else arg for cfg, arg in zip(self.CutConfig['event_range'], [min_event, max_event])]
        description = '{:1.0f}k - {:1.0f}k'.format(*self.get_event_range() / 1000.)
        return CutString('event_range', 'event_number>={} && event_number<={}'.format(*event_range), description)

    def generate_chi2(self, mode='x', value=None):
        cut_value = self.calc_chi2(mode) if value is None else value
        description = 'chi2 in {} < {:1.1f} ({:d}% quantile)'.format(mode, cut_value, self.CutConfig['chi2_{}'.format(mode)])
        return CutString('chi2_{}'.format(mode), 'chi2_{}>=0'.format(mode) + ' && chi2_{mod}<{val}'.format(val=cut_value, mod=mode) if cut_value is not None else '', description)

    def generate_track_angle(self, mode='x'):
        amin, amax = -self.CutConfig['track angle'], self.CutConfig['track angle']
        string = '{v}>{} && {v}<{}'.format(amin, amax, v='angle_{}'.format(mode))
        description = '{:1.1f} < tracking angle in {} < {:1.1f} [degrees]'.format(amin, mode, amax)
        return CutString('track_angle_{}'.format(mode), string if self.CutConfig['track angle'] > 0 else '', description)

    def generate_beam_interruptions(self):
        """ This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied. """
        interruptions = self.get_interruptions_ranges()
        cut_string = TCut('')
        for interr in interruptions:
            cut_string += TCut('event_number<{low}||event_number>{high}'.format(low=interr[0], high=interr[1]))
        description = '{} ({:.1f}% of the events excluded)'.format(len(interruptions), 100. * sum(j - i for i, j in interruptions) / self.Run.NEntries)
        return CutString('beam_interruptions', cut_string, description)

    def generate_aligned(self):
        """ Cut to exclude events with a wrong event alignment. """
        description = '{:.1f}% of the events excluded'.format(100. * self.find_n_misaligned() / self.Run.NEntries) if self.find_n_misaligned() else ''
        return CutString('aligned', 'aligned[0]' if self.find_n_misaligned() else '', description)

    def generate_fiducial(self, center=False, n_planes=0):
        if self.CutConfig['fiducial'] is None:
            return CutString('fiducial', '', '')
        xy = self.CutConfig['fiducial'] + (([self.Bins.PX / 2] * 2 + [self.Bins.PY / 2] * 2) if center else 0)
        cut = self.Analysis.draw_box(xy[0], xy[2], xy[1], xy[3], line_color=2, width=3, name='fid{}'.format(self.Run.Number), show=False)
        cut.SetVarX(self.get_track_var(self.Analysis.DUT.Number - 1 - n_planes, 'x'))
        cut.SetVarY(self.get_track_var(self.Analysis.DUT.Number - 1 - n_planes, 'y'))
        self.Analysis.add(cut)
        description = 'x: [{},{}], y: [{},{}], area: {:.1f}mm x {:.1f}mm = {:.1f}mm2'.format(*self.get_fiducial_size())
        return CutString('fiducial', TCut(cut.GetName()) if cut is not None else '', description)

    def generate_jump_cut(self):
        cut_string = ''
        start_event = self.CutConfig['event_range'][0]
        for tup in self.get_beam_interruptions():
            if tup[1] > start_event:
                low = start_event if tup[0] < start_event else tup[0]
                cut_string += '&&' if cut_string else ''
                cut_string += '!(event_number<={up}&&event_number>={low})'.format(up=tup[1], low=low)
        return TCut(cut_string)

    def generate_flux_cut(self):
        return self.generate_custom(include=['beam_interruptions', 'event_range'], name='flux', prnt=False)

    def generate_custom(self, exclude=None, include=None, name='custom', prnt=True):
        self.Analysis.info('generated {name} cut with {num} cuts'.format(name=name, num=self.CutStrings.get_n_custom(exclude, include)), prnt=prnt)
        return self.CutStrings.generate_custom(exclude, include, name)

    def generate_consecutive(self):
        return self.CutStrings.consecutive()
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    def calc_chi2(self, mode='x', quantile=None):
        picklepath = self.Analysis.make_pickle_path('Cuts', 'Chi2', run=self.Run.Number, suf=mode.title())

        def f():
            t = self.Analysis.info('calculating chi2 cut in {mod} for run {run}...'.format(run=self.Run.Number, mod=mode), next_line=False)
            values = get_root_vec(self.Analysis.Tree, var='chi2_{}'.format(mode))
            chi2s = get_quantiles(values[values > -500], linspace(0, 100, 501))
            self.Analysis.add_to_info(t)
            return chi2s

        chi2 = do_pickle(picklepath, f)
        q = self.CutConfig['chi2_{mod}'.format(mod=mode.lower())] if quantile is None else quantile
        return chi2[q] if q != 100 else None

    def get_raw_pulse_height(self):
        n = self.Analysis.Tree.Draw(self.Analysis.generate_signal_name(), self.CutStrings(), 'goff')
        return make_ufloat(mean_sigma(self.Run.get_root_vec(n)))

    def find_zero_ph_event(self, redo=False):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'EventMax', self.Run.Number, self.Analysis.DUT.Number)

        def f():
            t = self.Analysis.info('Looking for signal drops of run {} ...'.format(self.Run.Number), next_line=False)
            signal = self.Analysis.generate_signal_name()
            p = TProfile('pphc', 'Pulse Height Evolution', *self.Analysis.Bins.get_raw_time(30))
            self.Analysis.Tree.Draw('{}:{}>>pphc'.format(signal, self.Analysis.get_t_var()), self.CutStrings(), 'goff')
            values = array([p.GetBinContent(i) for i in xrange(1, p.GetNbinsX() + 1)])
            i_start = next(i for i, v in enumerate(values) if v) + 1  # find the index of the first bin that is not zero
            ph = mean(values[i_start:(values.size + 9 * i_start) / 10])  # take the mean of the first 10% of the bins
            i_break = next((i + i_start for i, v in enumerate(values[i_start:]) if v < .2 * ph and v), None)
            self.Analysis.add_to_info(t)
            return None if ph < 10 or i_break is None else self.Analysis.get_event_at_time(p.GetBinCenter(i_break - 2), rel=True)

        return do_pickle(pickle_path, f, redo=redo)

    def find_beam_interruptions(self):
        return self.find_pad_beam_interruptions() if self.Run.Type == 'pad' else self.find_pixel_beam_interruptions()

    def find_pad_beam_interruptions(self, bin_width=100, max_thresh=.6):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        t = self.Analysis.info('Searching for beam interruptions of run {r} ...'.format(r=self.Run.Number), next_line=False)
        n = self.Analysis.Tree.Draw('Entry$:pulser', '', 'goff')
        x, y = get_root_vecs(self.Analysis.Tree, n, 2, dtype='i4')
        rates, x_bins, y_bins = histogram2d(x, y, bins=[arange(0, n, bin_width, dtype=int), 2])
        rates = rates[:, 1] / bin_width
        thresh = min(max_thresh, mean(rates) + .2)
        events = x_bins[:-1][rates > thresh] + bin_width / 2
        not_connected = where(concatenate([[False], events[:-1] != events[1:] - bin_width]))[0]  # find the events where the previous event is not related to the event (more than a bin width away)
        events = split(events, not_connected)  # events grouped into connecting events
        interruptions = [(ev[0], ev[0]) if ev.size == 1 else (ev[0], ev[-1]) for ev in events] if events[0].size else []
        self.Analysis.add_to_info(t)
        return array(interruptions, 'i4')

    def find_pixel_beam_interruptions(self, bin_width=10, threshold=.4):
        """ Finding beam interruptions by incestigation the event rate. """
        t_start = self.Analysis.info('Searching for beam interruptions of run {r} ...'.format(r=self.Run.Number), next_line=False)
        bin_values, time_bins = histogram(self.Run.Time / 1000, bins=self.Bins.get_raw_time(bin_width)[1])
        m = mean(bin_values[bin_values.argsort()][-20:-10])  # take the mean of the 20th to the 10th highest bin to get an estimate of the plateau
        deviating_bins = where(abs(1 - bin_values / m) > threshold)[0]
        times = time_bins[deviating_bins] + bin_width / 2 - self.Run.Time[0] / 1000  # shift to the center of the bin
        not_connected = where(concatenate([[False], deviating_bins[:-1] != deviating_bins[1:] - 1]))[0]  # find the bins that are not consecutive
        times = split(times, not_connected)
        interruptions = [[self.Analysis.get_event_at_time(v) for v in [t[0], t[0] if t.size == 1 else t[-1]]] for t in times] if len(times[0]) else []
        self.Analysis.add_to_info(t_start)
        return interruptions

    def create_interruption_ranges(self, interruptions):
        ranges = []
        for i, tup in enumerate(interruptions):
            t_start = max(0, self.Run.get_time_at_event(tup[0]) - self.Run.StartTime - self.CutConfig['jump_range'][0])
            t_stop = self.Run.get_time_at_event(tup[1]) - self.Run.StartTime + self.CutConfig['jump_range'][1]
            # if interruptions overlay just set the last stop to the current stop
            if i and t_start <= (ranges[-1][1]) + 10:
                ranges[-1][1] = t_stop
                continue
            ranges.append([t_start, t_stop])
        return [[self.Run.get_event_at_time(t) for t in tup] for tup in ranges]

    def find_n_misaligned(self):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'align', self.Run.Number)

        def f():
            return where(get_root_vec(self.Analysis.Tree, var='aligned[0]', dtype=bool) == 0)[0].size
        return do_pickle(pickle_path, f)

    # endregion COMPUTE
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW & ANALYSE
    def show_cuts(self, raw=False):
        rows = [[cut.Name, '{:5d}'.format(cut.Level), cut.Value if raw else cut.Description] for cut in self.CutStrings.get_strings()]
        print_table([row for row in rows if row[2]], ['Cut Name', 'Level', 'Description'])

    def draw_contributions(self, flat=False, short=False, show=True):
        set_root_output(show)
        contr = OrderedDict()
        n_events = self.Run.NEntries
        cut_events = 0
        for i, (key, cut) in enumerate(self.generate_consecutive().iteritems()):
            if key == 'raw':
                continue
            events = n_events - int(self.Analysis.Tree.Draw('1', cut, 'goff'))
            print(key.rjust(18), '{0:5d} {1:04.1f}%'.format(events - cut_events, (1. - float(events) / n_events) * 100.))
            contr[key.title().replace('_', ' ')] = (events - cut_events, self.Analysis.get_color())
            cut_events = events
        contr['Good Events'] = (n_events - cut_events, self.Analysis.get_color())
        sorted_contr = OrderedDict(sorted(OrderedDict(item for item in contr.iteritems() if item[1][0] >= (.03 * n_events if short else 0)).iteritems(), key=lambda x: x[1]))  # sort by size
        sorted_contr.update({'Other': (n_events - sum(v[0] for v in sorted_contr.values()), self.Analysis.get_color())} if short else {})
        sorted_contr = OrderedDict(sorted_contr.popitem(not i % 2) for i in xrange(len(sorted_contr)))  # sort by largest->smallest->next largest...
        pie = TPie('pie', 'Cut Contributions', len(sorted_contr), array([v[0] for v in sorted_contr.values()], 'f'), array([v[1] for v in sorted_contr.values()], 'i'))
        for i, label in enumerate(sorted_contr.iterkeys()):
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label)
        format_pie(pie, h=.04, r=.2, text_size=.025, angle3d=70, label_format='%txt (%perc)', angle_off=250)
        self.Analysis.save_histo(pie, 'CutContributions', draw_opt='{0}rsc'.format('3d' if not flat else ''), show=show)
        self.Analysis.reset_colors()
        return sorted_contr

    def draw_fid_cut(self, scale=10):
        cut = get_object('fid{}'.format(self.Run.Number))
        if cut:
            cut = deepcopy(cut)
            cut.SetName('fid{}'.format(scale))
            for i in xrange(cut.GetN()):
                cut.SetPoint(i, scale * cut.GetX()[i], scale * cut.GetY()[i])
            cut.Draw()
            self.Analysis.Objects.append(cut)
    # endregion SHOW & ANALYSE
    # ----------------------------------------


class CutString:

    def __init__(self, name, value, description='', level=1):
        self.Name = name
        self.Value = str(value)
        self.Level = level
        self.Description = description

    def __call__(self, cut=None):
        return TCut(self.Name, self.Value) if cut is None else TCut('')

    def __str__(self):
        return '{:2d}: {} cut'.format(self.Level, self.Name.replace('_', ' '))

    def __repr__(self):
        return self.__str__()

    def reset(self):
        self.Value = ''

    def set(self, value):
        self.Value = value

    def set_level(self, level):
        self.Level = level
        return self


class CutStrings(object):

    def __init__(self):
        self.Strings = OrderedDict()

    def __call__(self):
        cut_string = TCut('AllCuts', '')
        for cut in self.get_strings():
            cut_string += cut()
        return cut_string

    def __getitem__(self, item):
        return self.Strings.values()[item]

    def register(self, cut, level):
        self.Strings[cut.Name] = cut.set_level(level)
        self.sort()

    def sort(self):
        self.Strings = OrderedDict(sorted(self.Strings.iteritems(), key=lambda x: x[1].Level))

    def names(self):
        return self.Strings.keys()

    def get(self, name):
        return self.Strings[name]() if self.has_cut(name) else warning('There is no cut with the name "{name}"!'.format(name=name))

    def get_strings(self):
        return [cut for cut in self.Strings.values() if cut.Value]

    def get_n(self):
        return sum(cut.Value != '' for cut in self.get_strings())

    def get_n_custom(self, exclude, include):
        return sum(cut.Value != '' for cut in self.get_strings() if cut.Name not in make_list(exclude) and (include is None or cut.Name in make_list(include)))

    def consecutive(self):
        cuts = OrderedDict([('raw', TCut('0', ''))])
        for i, cut in enumerate(self.get_strings(), 1):
            new_cut = cuts.values()[i - 1] + cut()
            cut.Name.replace('interruptions', 'stops')
            cuts[cut.Name] = TCut('{n}'.format(n=i), str(new_cut))
        return cuts

    def has_cut(self, name):
        return name in self.Strings

    def reset(self, name):
        self.Strings[name].reset() if self.has_cut(name) else warning('There is no cut with the name "{name}"!'.format(name=name))

    def reset_all(self):
        for cut in self.get_strings():
            cut.reset()

    def set(self, name, value):
        self.Strings[name].set(value) if self.has_cut(name) else warning('There is no cut with the name "{name}"!'.format(name=name))

    def set_description(self, name, txt):
        self.Strings[name].set_description(txt) if self.has_cut(name) else warning('There is no cut with the name "{name}"!'.format(name=name))

    def generate_custom(self, exclude=None, include=None, name='custom'):
        cut_string = TCut(name, '')
        for cut in [cut for cut in self.get_strings() if cut.Name not in make_list(exclude)]:
            if include is not None and cut.Name not in include or not cut.Value:
                continue
            cut_string += cut()
        return cut_string


def invert(cut):
    return TCut('!({})'.format(str(cut)))
