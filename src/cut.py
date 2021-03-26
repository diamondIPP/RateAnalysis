# --------------------------------------------------------
#       general class to handle all the cut strings for the analysis
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import quantile
from ROOT import TCut, TPie
from helpers.draw import *
from src.binning import Bins
from src.sub_analysis import SubAnalysis
from src.dut import Plane


class Cut(SubAnalysis):
    """ Contains methods to generate the cut strings for the TelescopeAnalysis and holds the dictionaries for the settings and all cut strings. """

    def __init__(self, analysis):

        super().__init__(analysis, pickle_dir='Cuts')

        # Configuration
        self.LowRateRun = None

        # Cut Strings
        self.CutStrings = CutStrings()
        self.HasFid = False

        # generate cut strings
        self.generate()

    def __call__(self, cut=None):
        return self.CutStrings() if cut is None else TCut('' if cut == Ellipsis else cut)

    def __getitem__(self, name):
        return self.get(name)

    def has(self, name):
        return bool(self.get(name).GetTitle())

    def get_name(self, cut=None):
        return 'NoCut' if cut == '' else self(cut).GetName() if not self(cut).GetName().startswith('All') else ''

    # ----------------------------------------
    # region CONFIG
    def update_config(self):
        pass

    def get_config(self, option, section='CUT', dtype=str, default=None):
        return self.Config.get_value(section, option, dtype, default)

    def load_fiducial(self, name='fiducial'):
        splits = array(self.Config.get_list('SPLIT', 'fiducial'))
        n = next(iter(where(splits > self.Run.Number)[0]), splits.size) + 1
        option = name if self.Config.has_option('CUT', name) and n == 1 else '{} {}'.format(name, n)
        v = self.load_dut_config(option, warn=False)
        if v is None:
            warning('fiducial cut not defined!')
            return make_box_args(-1, -1, 1, 1)
        self.HasFid = True
        return make_box_args(*array(v)[[0, 2, 1, 3]]) if len(v) == 4 else array(v)

    def load_dut_config(self, option, warn=True):
        if option not in self.Config.options('CUT') or self.Ana.DUT.Name not in self.Config.get_list('CUT', option):
            return warning('No option {} in the analysis config for {}!'.format(option, make_tc_str(self.TCString)), prnt=warn)
        return self.Config.get_list('CUT', option)[self.Ana.DUT.Name]
    # endregion CONFIG
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, name, invert=False, warn=True):
        return self.CutStrings.get(name, invert, warn)

    def get_strings(self, with_raw=False):
        return self.CutStrings.get_strings(with_raw)

    def get_all(self, with_raw=False):
        return [cut() for cut in self.get_strings(with_raw)]

    def get_names(self, with_raw=False):
        return self.CutStrings.get_names(with_raw)

    def get_consecutive(self, short=False):
        cuts = self.CutStrings.consecutive()
        return {key: cut for key, cut in cuts.items() if key in self.get_short()} if short else cuts

    def get_size(self, name, excluded=True):
        n = self.Ana.get_n_entries(self.get(name) if name in self.get_names() else name)
        return self.Run.NEvents - n if excluded else n

    def get_sizes(self, consecutive=True, redo=False):
        def f():
            return array([self.get_size(cut) for cut in (self.get_consecutive().values() if consecutive else self.get_all())])
        return do_pickle(self.Ana.make_simple_pickle_path('Sizes', int(consecutive), 'Cuts'), f, redo=redo)

    def get_short(self, n=6, redo=False):
        """get a list of names of the n biggest cuts"""
        return ['raw'] + list(array(self.get_names())[diff(self.get_sizes(redo=redo)).argsort()][-n:])

    def get_event_range(self):
        """ :return: event range. Negative values are interpreted as minutes. Example: [-10, 700k] => 10 min < events < 700k, type [ndarray]"""
        low, high = [self.Ana.Run.get_event_at_time(seconds=abs(v * 60)) if v < 0 else v for v in self.Config.get_list('CUT', 'event range', default=[0, 0])]
        return array([low, self.Run.NEvents if high == 0 else high])

    def get_track_angle(self):
        return self.get_config('track angle', dtype=int)

    def get_chi2(self, mode='x', value=None):
        return choose(value, self.get_config('chi2{}'.format(mode.title()), dtype=int))

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
        return 'dia_track_{m}_local[{n}]{s}'.format(m=mode, n=num, s=' * 10' if mm else '')

    def get_track_vars(self, num, mm=False):
        return [self.get_track_var(num, v, mm) for v in ['x', 'y']]

    def get_beam_interruptions(self):
        """ :returns: list of raw interruptions, type [list[tup]]"""
        return do_pickle(self.make_simple_pickle_path('', sub_dir='BeamStops'), self.find_beam_interruptions)

    def get_interruptions_ranges(self):
        """ :returns: list of interruptions including safety margin from the AnalysisConfig. """
        pickle_path = self.make_simple_pickle_path('BeamStops', self.Config.get_value('CUT', 'exclude around jump', dtype=list))
        return do_pickle(pickle_path, self.create_interruption_ranges, interruptions=self.get_beam_interruptions())

    def get_fiducial_area(self):
        return poly_area(*self.load_fiducial())

    def get_fid(self, scale=10):
        cut = get_object('fid{}'.format(self.Run.Number))
        if cut:
            cut = deepcopy(cut)
            cut.SetName('fid{}{}'.format(self.Run.Number, scale))
            for i in range(cut.GetN()):
                cut.SetPoint(i, scale * cut.GetX()[i], scale * cut.GetY()[i])
            return Draw.add(cut)

    def get_raw_pulse_height(self):
        return mean_sigma(self.Run.get_tree_vec(var=self.Ana.get_signal_var(), cut=self()), err=False)[0]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SET
    def set(self, name, value):
        self.CutStrings.set(name, value)

    def set_fiducial(self, name):
        cut = self.generate_fiducial(fid_name=name)
        self.set('fiducial', cut.Value)
        self.CutStrings.set_description('fiducial', cut.Description)

    def set_low_rate_run(self, low_run):
        self.LowRateRun = str(low_run)

    def reset(self, name):
        self.CutStrings.reset(name)

    def update(self, name, value=None):
        self.CutStrings.set(name, value)

    def set_chi2(self, value):
        self.update('chi2_x', self.generate_chi2('x', value).Value)
        self.update('chi2_y', self.generate_chi2('y', value).Value)

    def reload(self):
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
        event_range = [choose(arg, cfg) for cfg, arg in zip(self.get_event_range(), [min_event, max_event])]
        description = '{:1.0f}k - {:1.0f}k'.format(*self.get_event_range() / 1000.)
        return CutString('event range', 'event_number>={} && event_number<={}'.format(*event_range), description)

    def generate_chi2s(self, q=None):
        string = (self.generate_chi2('x', q) + self.generate_chi2('y', q)).Value
        return CutString('chi2', string, 'chi2 x&y < {}'.format(self.get_chi2(value=q)))

    def generate_chi2(self, mode='x', q=None):
        cut_value = self.calc_chi2(mode, q)
        if cut_value is None:
            return CutString('chi2_{}'.format(mode), '')
        description = 'chi2 in {} < {:1.1f} ({:d}% quantile)'.format(mode, cut_value, self.get_chi2(mode, q))
        return CutString('chi2_{}'.format(mode), 'chi2_{}>=0'.format(mode) + ' && chi2_{mod}<{val}'.format(val=cut_value, mod=mode), description)

    def generate_track_angle(self, mode='x', amin=None, amax=None):
        amin, amax = (array([-1, 1]) * self.get_track_angle()) if amin is None else (amin, amax)
        string = '{v}>{} && {v}<{}'.format(amin, amax, v='angle_{}'.format(mode))
        description = '{:1.1f} < tracking angle in {} < {:1.1f} [degrees]'.format(amin, mode, amax)
        return CutString('track_angle_{}'.format(mode), string if amax > 0 else '', description)

    def generate_beam_interruptions(self):
        """ This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied. """
        interruptions = self.get_interruptions_ranges()
        cut_string = TCut('')
        for interr in interruptions:
            cut_string += TCut('event_number<{low}||event_number>{high}'.format(low=interr[0], high=interr[1]))
        description = '{} ({:.1f}% of the events excluded)'.format(len(interruptions), 100. * sum(j - i for i, j in interruptions) / self.Run.NEvents)
        return CutString('beam stops', cut_string, description)

    def generate_aligned(self):
        """ Cut to exclude events with a wrong event alignment. """
        description = '{:.1f}% of the events excluded'.format(100. * self.find_n_misaligned() / self.Run.NEvents) if self.find_n_misaligned() else ''
        return CutString('aligned', 'aligned[0]' if self.find_n_misaligned() else '', description)

    def generate_fiducial(self, center=False, n_planes=0, x=None, y=None, name=None, fid_name='fiducial', xvar=None, yvar=None):
        x = choose(x, self.load_fiducial(fid_name)[0]) + (Plane.PX / 20 if center else 0)
        y = choose(y, self.load_fiducial(fid_name)[1]) + (Plane.PY / 20 if center else 0)
        cut = Draw.polygon(x, y, line_color=2, width=3, name=choose(name, 'fid{}'.format(self.Run.Number)), show=False)
        cut.SetVarX(choose(xvar, self.get_track_var(self.Ana.DUT.Number - 1 - n_planes, 'x')))
        cut.SetVarY(choose(yvar, self.get_track_var(self.Ana.DUT.Number - 1 - n_planes, 'y')))
        description = 'x: {}, y: {}, area: {:.1f} mm2'.format(x * 10, y * 10, poly_area(x, y) * 100)
        return CutString(choose(name, 'fiducial'), TCut(cut.GetName()) if cut is not None else '', description)

    def generate_sub_fid(self, name, x1, x2, y1, y2):
        x, y = make_box_args(x1, y1, x2, y2)
        return self.generate_fiducial(x=x, y=y, name=name)() + self()

    def generate_jump_cut(self):
        cut_string = ''
        start_event = self.get_min_event()
        for tup in self.get_beam_interruptions():
            if tup[1] > start_event:
                low = start_event if tup[0] < start_event else tup[0]
                cut_string += '&&' if cut_string else ''
                cut_string += '!(event_number<={up}&&event_number>={low})'.format(up=tup[1], low=low)
        return TCut(cut_string)

    def generate_flux_cut(self):
        return self.generate_custom(include=['beam stops', 'event range'], name='flux', prnt=False)

    def generate_custom(self, exclude=None, include=None, invert=None, name='custom', add=None, prnt=True):
        self.Ana.info('generated {name} cut with {num} cuts'.format(name=name, num=self.CutStrings.get_n_custom(exclude, include)), prnt=prnt)
        return self.CutStrings.generate_custom(exclude, include, invert, name) + TCut(Cut.to_string(add))
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    def calc_chi2(self, mode='x', q=None):
        def f():
            t = self.Ana.info('calculating chi2 cut in {mod} for run {run} ...'.format(run=self.Run.Number, mod=mode), endl=False)
            values = self.Ana.get_tree_vec('chi2_{}'.format(mode))
            chi2s = quantile(values[values > -500], linspace(0, 1, 101))
            self.Ana.add_to_info(t)
            return chi2s
        chi2 = do_hdf5(self.Ana.make_simple_hdf5_path('Chi2{}'.format(mode.title()), sub_dir='Cuts'), f)
        q = choose(q, self.get_chi2(mode))
        return chi2[q] if q != 100 else None

    def find_zero_ph_event(self, redo=False):
        pickle_path = self.Ana.make_pickle_path('Cuts', 'EventMax', self.Run.Number, self.Ana.DUT.Number)

        def f():
            t0 = self.Ana.info('Looking for signal drops of run {} ...'.format(self.Run.Number), endl=False)
            ph, t = self.Ana.get_tree_vec(var=[self.Ana.get_signal_name(), self.Ana.get_t_var()], cut=self())
            time_bins, values = get_hist_vecs(self.Ana.Draw.profile(t, ph, Bins(self.Ana).get_raw_time(30), show=False), err=False)
            i_start = next(i for i, v in enumerate(values) if v) + 1  # find the index of the first bin that is not zero
            ph = abs(mean(values[i_start:(values.size + 9 * i_start) // 10]))  # take the mean of the first 10% of the bins
            i_break = next((i + i_start for i, v in enumerate(values[i_start:]) if abs(v) < .2 * ph and v), None)
            self.Ana.add_to_info(t0)
            return None if ph < 10 or i_break is None else self.Ana.get_event_at_time(time_bins[i_break - 1])

        return do_pickle(pickle_path, f, redo=redo)

    def find_beam_interruptions(self):
        return self.find_pixel_beam_interruptions()

    def find_pixel_beam_interruptions(self, bin_width=10, threshold=.4):
        """ Finding beam interruptions by incestigation the event rate. """
        t_start = self.Ana.info('Searching for beam interruptions of run {r} ...'.format(r=self.Run.Number), endl=False)
        bin_values, time_bins = histogram(self.Run.Time / 1000, bins=self.Bins.get_raw_time(bin_width)[1])
        m = mean(bin_values[bin_values.argsort()][-20:-10])  # take the mean of the 20th to the 10th highest bin to get an estimate of the plateau
        deviating_bins = where(abs(1 - bin_values / m) > threshold)[0]
        times = time_bins[deviating_bins] + bin_width / 2 - self.Run.Time[0] / 1000  # shift to the center of the bin
        not_connected = where(concatenate([[False], deviating_bins[:-1] != deviating_bins[1:] - 1]))[0]  # find the bins that are not consecutive
        times = split(times, not_connected)
        interruptions = [[self.Ana.get_event_at_time(v) for v in [t[0], t[0] if t.size == 1 else t[-1]]] for t in times] if len(times[0]) else []
        self.Ana.add_to_info(t_start)
        return interruptions

    def create_interruption_ranges(self, interruptions):
        ranges = []
        low, high = self.Config.get_value('CUT', 'exclude around jump', dtype=list)
        for i, tup in enumerate(interruptions):
            t_start = max(0, self.Run.get_time_at_event(tup[0]) - self.Run.StartTime - low)
            t_stop = self.Run.get_time_at_event(tup[1]) - self.Run.StartTime + high
            # if interruptions overlay just set the last stop to the current stop
            if i and t_start <= (ranges[-1][1]) + 10:
                ranges[-1][1] = t_stop
                continue
            ranges.append([t_start, t_stop])
        return [[self.Run.get_event_at_time(t) for t in tup] for tup in ranges]

    def find_n_misaligned(self):
        pickle_path = self.Ana.make_pickle_path('Cuts', 'align', self.Run.Number)

        def f():
            return where(get_tree_vec(self.Ana.Tree, var='aligned[0]', dtype=bool) == 0)[0].size
        return do_pickle(pickle_path, f)

    # endregion COMPUTE
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW & ANALYSE
    def show(self, raw=False):
        rows = [[cut.Name, '{:5d}'.format(cut.Level), cut.Value if raw else cut.Description] for cut in self.CutStrings.get_strings()]
        print_table([row for row in rows if row[2]], ['Cut Name', 'Level', 'Description'])

    def draw_contributions(self, flat=False, short=False, show=True):
        set_root_output(show)
        contr = OrderedDict()
        n, n_events = len(self.get_consecutive()), self.Run.NEvents
        cut_events = 0
        for i, (key, cut) in enumerate(self.get_consecutive().items()):
            if key == 'raw':
                continue
            events = n_events - int(self.Ana.Tree.Draw('1', cut, 'goff'))
            print(key.rjust(18), '{0:5d} {1:04.1f}%'.format(events - cut_events, (1. - float(events) / n_events) * 100.))
            contr[key.title().replace('_', ' ')] = (events - cut_events, self.Draw.get_color(n))
            cut_events = events
        contr['Good Events'] = (n_events - cut_events, self.Draw.get_color(n))
        sorted_contr = OrderedDict(sorted(OrderedDict(item for item in contr.items() if item[1][0] >= (.03 * n_events if short else 0)).items(), key=lambda x: x[1]))  # sort by size
        sorted_contr.update({'Other': (n_events - sum(v[0] for v in sorted_contr.values()), self.Draw.get_color(n))} if short else {})
        sorted_contr = OrderedDict(sorted_contr.popitem(not i % 2) for i in range(len(sorted_contr)))  # sort by largest->smallest->next largest...
        pie = TPie('pie', 'Cut Contributions', len(sorted_contr), array([v[0] for v in sorted_contr.values()], 'f'), array([v[1] for v in sorted_contr.values()], 'i'))
        for i, label in enumerate(sorted_contr.keys()):
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label)
        format_pie(pie, h=.04, r=.2, text_size=.025, angle3d=70, label_format='%txt (%perc)', angle_off=250)
        self.Draw(pie, 'CutContributions', draw_opt='{0}rsc'.format('3d' if not flat else ''), show=show)
        return sorted_contr

    def draw_fid(self, scale=10):
        self.get_fid(scale).Draw()
    # endregion SHOW & ANALYSE
    # ----------------------------------------

    @staticmethod
    def invert(cut):
        return TCut('!({})'.format(Cut.to_string(cut)))

    @staticmethod
    def to_string(cut):
        return cut.GetTitle() if type(cut) is TCut else choose(cut, '')


class CutString:

    def __init__(self, name, value, description='', level=1):
        self.Name = name
        self.Value = Cut.to_string(value)
        self.Level = level
        self.Description = description

    def __call__(self, cut=None):
        return TCut(self.Name, self.Value) if cut is None else TCut(cut)

    def __str__(self):
        return '{:2d}: {} cut'.format(self.Level, self.Name.replace('_', ' '))

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if other is not None:
            self.Value = (self() + (other() if type(other) is CutString else TCut(other))).GetTitle()
        return self

    def reset(self):
        self.Value = ''

    def set(self, value):
        self.Value = value

    def set_description(self, txt):
        self.Description = txt

    def set_level(self, level):
        self.Level = level
        return self

    def invert(self):
        return TCut('!{}'.format(self.Name), '!({})'.format(self.Value))


class CutStrings(object):

    def __init__(self):
        self.Strings = OrderedDict()

    def __call__(self):
        cut_string = TCut('AllCuts', '')
        for cut in self.get_strings():
            cut_string += cut()
        return cut_string

    def __getitem__(self, item):
        return self.get(item)

    def register(self, cut, level):
        if cut.Value:
            self.Strings[cut.Name] = cut.set_level(level)
            self.sort()

    def sort(self):
        self.Strings = OrderedDict(sorted(self.Strings.items(), key=lambda x: x[1].Level))

    def get_names(self, with_raw=False):
        return (['raw'] if with_raw else []) + list(self.Strings)

    def get(self, name, invert=False, warn=True):
        if self.has_cut(name):
            return self.Strings[name].invert() if invert else self.Strings[name]()
        warning('There is no cut with the name "{name}"!'.format(name=name), prnt=warn)
        return ''

    def get_strings(self, with_raw=False):
        return ([CutString('raw', '', 'raw values')] if with_raw else []) + list(self.Strings.values())

    def get_n(self):
        return sum(cut.Value != '' for cut in self.get_strings())

    def get_n_custom(self, exclude, include):
        return sum(cut.Value != '' for cut in self.get_strings() if cut.Name not in make_list(exclude) and (include is None or cut.Name in make_list(include)))

    def consecutive(self):
        cuts = OrderedDict([('raw', TCut('0', ''))])
        for i, cut in enumerate(self.get_strings(), 1):
            new_cut = list(cuts.values())[i - 1] + cut()
            cut.Name.replace('interruptions', 'stops')
            cuts[cut.Name] = TCut('{n}'.format(n=i), new_cut.GetTitle())
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

    def generate_custom(self, exclude=None, include=None, invert=None, name='custom'):
        cut_string = TCut(name, '')
        for cut in [cut for cut in self.get_strings() if cut.Name not in make_list(exclude)]:
            if include is not None and cut.Name not in include or not cut.Value:
                continue
            cut_string += cut() if cut.Name not in make_list(invert) else cut.invert()
        return cut_string
