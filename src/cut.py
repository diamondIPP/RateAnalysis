# --------------------------------------------------------
#       general class to handle all the cut strings for the analysis
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TCut
from numpy import delete, insert, vectorize

from plotting.draw import *
from helpers.utils import *
from src.binning import Bins
from src.dut import Plane
from src.sub_analysis import SubAnalysis
import plotting.latex as latex


class Cut(SubAnalysis):
    """ Contains methods to generate the cut strings for the TelescopeAnalysis and holds the dictionaries for the settings and all cut strings. """

    M = ['x', 'y']

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

    def __repr__(self):
        return f'{self.__class__.__name__} of {self.Run}'

    def has(self, name):
        return bool(self.get(name, warn=False).GetTitle())

    def get_name(self, cut=None):
        return 'NoCut' if cut == '' else self(cut).GetName() if not self(cut).GetName().startswith('All') else ''

    # ----------------------------------------
    # region CONFIG
    def update_config(self):
        pass

    def get_config(self, option, section='CUT', dtype: type = str, default=None, required=False):
        return critical(f'{option} cut not defined!') if required and not self.Config.has_option(section, option) else self.Config.get_value(section, option, dtype, default)

    def load_fiducial(self, name='fiducial'):
        splits = array(self.Config.get_list('SPLIT', 'fiducial'))
        n = next(iter(where(splits > self.Run.Number)[0]), splits.size) + 1
        option = name if self.Config.has_option('CUT', name) and n == 1 else '{} {}'.format(name, n)
        v = array(self.load_dut_config(option, warn=False))
        if v.ndim == 0:
            warning('fiducial cut not defined!')
            return make_box_args(-1, -1, 1, 1)
        self.HasFid = True
        return make_box_args(*v[[0, 2, 1, 3]]) if v.ndim == 1 else v

    def load_pixel_fid(self):
        fid = self.load_fiducial('pixel fiducial')
        return fid if fid[0][0] != -1 else make_box_args(0, 0, Plane.NCols - 1, Plane.NRows - 1)

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

    def get_consecutive(self, short=False, raw=True):
        cuts = self.CutStrings.consecutive()
        return {key: cut for key, cut in cuts.items() if key in self.get_short()} if short else {k: v for k, v in cuts.items() if raw or 'raw' not in k}

    def get_size(self, name, excluded=True):
        n = self.Ana.get_n_entries(self.get(name) if name in self.get_names() else name)
        return self.Run.NEvents - n if excluded else n

    @save_pickle('Sizes', suf_args='all')
    def get_sizes(self, consecutive=True, _redo=False):
        return array([self.get_size(cut) for cut in (self.get_consecutive().values() if consecutive else self.get_all())])

    def get_short(self, n=6, redo=False):
        """:returns a list of names of the <n> biggest cuts"""
        return ['raw'] + list(array(self.get_names())[diff(self.get_sizes(_redo=redo)).argsort()][-n:])

    def get_event_range(self):
        """ :return: event range. Negative values are interpreted as minutes. Example: [-10, 700k] => 10 min < events < 700k, type [ndarray]"""
        low, high = [self.Ana.Run.get_event_at_time(seconds=abs(v * 60)) if v < 0 else v for v in self.Config.get_list('CUT', 'event range', default=[0, 0])]
        return array([low, self.Run.NEvents if high == 0 else high])

    def get_track_angle(self, sigma=None):
        return sum([self.generate_track_angle(m, sigma)() for m in range(2)], start=TCut('track angle', ''))

    def get_chi2(self, m=0, value=None):
        return choose(value, self.get_config(f'chi2{Cut.M[m]}', dtype=int, required=True))

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
    def get_track_var(n, mode, mm=False, local=True, pixel=False):
        if pixel:
            return f'dia_track_{"col" if mode == "x" else "row"}[{n}]'
        return f'dia_track_{mode}{"_local" if local and not pixel else ""}[{n}]{" * 10" if mm else ""}'

    def get_track_vars(self, num, mm=False, local=True, pixel=False):
        return [self.get_track_var(num, v, mm, local, pixel) for v in ['x', 'y']]

    def get_fiducial_area(self):
        return poly_area(*self.load_fiducial())

    def get_fid(self, scale=10, name=None):
        cut = get_object(choose(name, f'fid{self.Run.Number}'))
        if cut:
            cut = deepcopy(cut)
            cut.SetName('fid{}{}'.format(self.Run.Number, scale))
            cut.SetFillStyle(1)
            for i in range(cut.GetN()):
                cut.SetPoint(i, scale * cut.GetX()[i], scale * cut.GetY()[i])
            return Draw.add(cut)

    def no_fid(self, fid=False, cut=None):
        return self(cut) if cut is not None or fid else self.exclude('fiducial')

    def no_chi2(self, cut=None):
        return self.exclude('chi2 x', 'chi2 y') if cut is None else self(cut)

    def ph(self, high, low=None):
        low, high = (high, low) if low is not None and low > high else (low, high)
        return Cut.make('', f'{self.Ana.get_signal_var()} < {high}') + (f'{self.Ana.get_signal_var()} > {low}' if low is not None else '')

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

    def remove(self, *name, ret=False):
        for cut in name:
            self.CutStrings.remove(cut)
        return self if ret else self.Ana.Bins.remove_pickle()

    def set_chi2(self, value):
        self.update('chi2 x', self.generate_chi2(0, value).Value)
        self.update('chi2 y', self.generate_chi2(1, value).Value)

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
        self.CutStrings.register(self.generate_chi2(0), 72)
        self.CutStrings.register(self.generate_chi2(1), 73)
        self.CutStrings.register(self.generate_track_angle(0), 74)
        self.CutStrings.register(self.generate_track_angle(1), 75)

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
        string = sum([self.generate_chi2(i, q).Value for i in range(2)], start=TCut())
        return CutString('chi2', string, f'chi2 x & y < {self.get_chi2(value=q)}')

    def generate_chi2(self, m=0, q=None):
        v, q, m = self.calc_chi2(m, q), self.get_chi2(m, q), Cut.M[m]
        if v is None:
            return CutString('chi2', '')
        return CutString(f'chi2 {m}', Cut.sum(f'chi2_{m} >= 0', f'chi2_{m} < {v}'), description=f'chi2 in {m} < {v:1.1f} ({q:d}% quantile)')

    def generate_track_angle(self, mode=0, sigma=None):
        n = choose(sigma, self.get_config('track angle sigma', dtype=float, required=True))
        m, s = self.calc_angle(mode)
        amin, amax = m - n * s, m + n * s
        string = f'angle_{Cut.M[mode]} > {amin} && angle_{Cut.M[mode]} < {amax}'
        description = f'{amin:1.1f} < tracking angle in {Cut.M[mode]} < {amax:1.1f} [degrees] ({n} sigma)'
        return CutString(f'track angle {Cut.M[mode]}', string, description)

    def generate_beam_interruptions(self):
        """ This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied. """
        interruptions = self.create_interruption_ranges()
        cut_string = TCut('')
        for interr in interruptions:
            cut_string += TCut('event_number<{} || event_number>{}'.format(*interr))
        description = '{} ({:.1f}% of the events excluded)'.format(len(interruptions), 100. * sum(j - i for i, j in interruptions) / self.Run.NEvents)
        return CutString('beam stops', cut_string, description)

    def generate_aligned(self):
        """ Cut to exclude events with a wrong event alignment. """
        if not self.has_branch('aligned') or self.Tree.GetBranch('aligned').ClassName() == 'TBranchElement':
            return CutString()
        description = '{:.1f}% of the events excluded'.format(100. * self.find_n_misaligned() / self.Run.NEvents) if self.find_n_misaligned() else ''
        return CutString('aligned', 'aligned' if self.find_n_misaligned() else '', description)

    def generate_fiducial(self, center=False, x=None, y=None, name=None, fid_name='fiducial', xvar=None, yvar=None):
        x = choose(x, self.load_fiducial(fid_name)[0]) + (Plane.PX / 20 if center else 0)
        y = choose(y, self.load_fiducial(fid_name)[1]) + (Plane.PY / 20 if center else 0)
        cut = Draw.polygon(x, y, line_color=2, width=3, name=choose(name, 'fid{}'.format(self.Run.Number)), show=False)
        cut.SetVarX(choose(xvar, self.get_track_var(self.Ana.DUT.Number - 1, 'x', mm=False)))
        cut.SetVarY(choose(yvar, self.get_track_var(self.Ana.DUT.Number - 1, 'y', mm=False)))
        description = 'x: {}, y: {}, area: {:.1f} mm2'.format(x * 10, y * 10, poly_area(x, y) * 100)
        return CutString(choose(name, 'fiducial'), TCut(cut.GetName()) if cut is not None else '', description)

    def generate_sub_fid(self, name, x1, x2, y1, y2):
        x, y = make_box_args(x1, y1, x2, y2)
        return self.generate_fiducial(x=x, y=y, name=name)() + self()

    def generate_jump_cut(self):
        cut_string = ''
        start_event = self.get_min_event()
        for tup in self.find_beam_interruptions():
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

    def exclude(self, *exclude, name='custom'):
        return self.CutStrings.generate_custom(exclude, name=name)

    def include(self, *include, name='incl'):
        return self.CutStrings.generate_custom(include=include, name=name)

    def add(self, *add, name='custom'):
        return TCut(name, Cut.to_string(self() + sum([Cut.make('', cut) for cut in add], start=TCut())))
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    @save_pickle('Chi2', print_dur=True, suf_args=0)
    def calc_chi2_(self, m=0, _redo=False):
        x = self.Ana.get_tree_vec(f'chi2_{Cut.M[m]}')
        return quantile(x[(x > -998) & (x < 100)], linspace(0, 1, 101))

    def calc_chi2(self, m=0, q=None, redo=False):
        q = choose(q, self.get_chi2(m))
        return self.calc_chi2_(m, _redo=redo)[q] if q != 100 else None

    @save_pickle('Angle', print_dur=True, suf_args='all')
    def calc_angle(self, m=0, _redo=False):
        x = self.get_tree_vec(f'angle_{Cut.M[m]}', cut=self())
        return mean_sigma(x[x > -900], err=False)

    @save_pickle('EventMax', print_dur=True)
    def find_signal_drops(self, thresh=.6, pol=1, _redo=False):
        ph, t = self.Ana.get_tree_vec(var=[self.Ana.get_signal_name(), self.Ana.get_t_var()], cut=self())
        p = self.Draw.profile(t, ph, Bins(self.Ana).get_raw_time(30), show=False)
        (x, y), n = hist_xy(p, err=False), bins.entries(p)
        x, y = x[(y != 0) & (n > 5)], pol * y[(y != 0) & (n > 5)]
        averages = cumsum(y) / (arange(y.size) + 1)
        j = next((i for i, is_bad in enumerate(y[2:] < thresh * averages[1:-1], 2) if is_bad), None)  # find next entry that is below the average of the previous
        sleep(.05)  # otherwise j is not initialised yet...
        return None if averages[min(1, y.size // 10)] < 10 or j is None else self.get_event_at_time(x[j - 2], rel=False)

    def find_beam_interruptions(self):
        return [[]]

    @save_pickle('BeamStopRange', print_dur=True)
    def create_interruption_ranges(self, _redo=False):
        low, high = self.Config.get_list('CUT', 'exclude around jump')
        x = array([[self.Run.get_time_at_event(ev) for ev in tup] for tup in self.find_beam_interruptions()]) - self.Run.StartTime
        x = (x + [-low, high] if x.size else x).flatten()  # extend range and flatten array
        x[x < 0] = 0  # time cannot be smaller than 0
        i = where(diff(x) < 10)[0]  # find all indices where the start of the next beam interruption is earlier than the end of the previous
        x = array([self.Run.get_event_at_time(t) for t in delete(x, concatenate([i, i + 1]))])  # remove the last stop and the next start if overlap
        return x.reshape((x.size // 2, 2))

    def interruption_times(self):
        return vectorize(self.Run.ev2t)(self.create_interruption_ranges())

    @save_pickle('Align')
    def find_n_misaligned(self, _redo=False):
        return where(self.get_tree_vec('aligned', dtype=bool) == 0)[0].size
    # endregion COMPUTE
    # ----------------------------------------

    # ----------------------------------------
    # region CONTRIBUTIONS
    def show(self, raw=False, latex_=False):
        rows, header = [[cut.Name, '{:5d}'.format(cut.Level), cut.Value if raw else cut.Description] for cut in self.CutStrings.get_strings()], ['Cut Name', 'Level', 'Description']
        rows = [row for row in rows if row[2]]
        print(latex.table(latex.bold(*header), rows)) if latex_ else print_table(rows, header)

    @update_pbar
    def get_contribution(self, cut, n_previous=0):
        return self.Run.NEvents - self.Tree.GetEntries(self(cut).GetTitle()) - n_previous

    @save_pickle('Contribution')
    def _get_contributions(self, _redo=False):
        cuts = self.get_consecutive(raw=False)
        PBAR.start(len(self.get_consecutive()) - 1, counter=True)
        n = [self.get_contribution(cut) for cut in self.get_consecutive().values()]
        return {name: i for name, i in zip(cuts, diff(n))}

    def get_contributions(self, threshold=None, redo=False):
        return {name: n for name, n in self._get_contributions(_redo=redo).items() if threshold is None or n > threshold * self.Run.NEvents}

    def show_contributions(self, redo=False, latex_=False, ret=False):
        abs_vals = 100 * (1 - (cumsum(list(self.get_contributions(redo=redo).values()))) / self.Run.NEvents)
        rows = [[name, f'{value:>6}', f'{value / self.Run.NEvents * 100: 10.2f}', f'{abs_: 3.2f}'] for (name, value), abs_ in zip(self.get_contributions().items(), abs_vals)]
        header = ['Cut', 'Events', 'Contr. [%]', 'Abs [%]']
        return rows if ret else print(latex.table(latex.bold(*header), rows)) if latex_ else None if print_table(rows, header=header) else None

    def show_2contributions(self, run, redo=False):
        from analyse import analysis_selector
        r1, r2 = array(self.show_contributions(redo, ret=True)), array(analysis_selector(run, self.DUT.Number, self.TCString, tree=True, prnt=False).Cut.show_contributions(redo, ret=True))
        r1 = insert(r1, where(r2[:, 0] == 'aligned')[0][0], ['aligned', '0', '      0.00', '0'], axis=0) if len(r1) != len(r2) and 'aligned' in r2 else r1
        r1, r2 = [append(r, [['good events', '0', f'{r[-1][-1]:>10}', '1']], axis=0) for r in [r1, r2]]
        print(latex.table(latex.bold('Cut Name', 'Low Rate', 'High Rate'), column_stack((array(r1)[:, [0, 2]], array(r2)[:, 2]))))

    def draw_contributions(self, flat=False, threshold=None, redo=False, **kwargs):
        cuts = self.get_contributions(threshold=threshold / 100, redo=redo)
        cuts['good events'] = self.Run.NEvents - sum(self.get_contributions().values())
        cuts = {**({'other': self.Run.NEvents - sum(cuts.values())} if threshold is not None else {}), **cuts}
        cuts = OrderedDict(sorted({name: (n, self.Draw.get_color(len(cuts))) for name, n in cuts.items()}.items(), key=lambda x: x[1]))  # add colors and sort by number of events
        sorted_cuts = OrderedDict(cuts.popitem(i % 2 == 1) for i in range(len(cuts)))  # sort by largest->smallest->next largest...
        self.Draw.pie(sorted_cuts, **prep_kw(kwargs, offset=.05, flat=flat, h=.04, r=.2, text_size=.025, angle3d=70, label_format='%txt (%perc)', angle_off=250))
        return sorted_cuts

    def draw_fid(self, scale=10):
        f = self.get_fid(scale)
        f.Draw()
        return f

    def get_p(self, cut, redo=False):
        cut = Cut.to_string(cut)
        n, n0 = [self.Ana.get_n_entries(c, _redo=redo) for c in [self(), self.exclude(cut, name=f'p-{cut}')]]
        return (n0 - n) / n

    @save_pickle('Purity')
    def purity(self, _redo=False):
        n = self.Ana.get_n_entries(self())
        return ufloat(n, sqrt(n - 1)) / self.Run.NEvents
    # endregion CONTRIBUTIONS
    # ----------------------------------------

    # ----------------------------------------
    # region ANA
    def compare_single_cuts(self, scale=True, redo=False):
        histos = [self.Ana.draw_signal_distribution(cut(), show=False, redo=redo, save=False) for cut in self.get_strings()]
        self.Draw.stack(histos, 'Single Cuts', self.get_names(), scale=scale)

    def compare_consecutive_cuts(self, scale=False, short=False, x_range=None, redo=False):
        cuts = self.get_consecutive(short)
        histos = [self.Ana.draw_signal_distribution(cut=cut, show=False, redo=redo, x_range=x_range, save=False) for cut in cuts.values()]
        self.Draw.stack(histos, 'Signal Distribution with Consecutive Cuts', cuts.keys(), scale)

    @quiet
    def draw_means(self, short=False, cuts=None, names=None, normalise=True, redo=False, **dkw):
        cuts, labels = choose(cuts, list(self.get_consecutive(short).values())), choose(names, self.get_consecutive(short).keys())
        PBAR.start(len(cuts), counter=True) if redo or not file_exists(self.make_simple_pickle_path('Fit', f'{cuts[-1].GetName()}_1', 'PH')) else do_nothing()
        x, y = arange(len(cuts)), array([self.Ana.get_pulse_height(cut=cut, redo=redo) for cut in cuts])
        y /= y[-1] if normalise else 1
        ytit = {'y_tit': 'Relative Pulse Height' if normalise else self.Ana.PhTit}
        return self.Draw.graph(x, y, title='PH Consec. Cuts', **ytit, **prep_kw(dkw, draw_opt='ap', gridy=True, x_range=[-1, len(y)], bin_labels=labels, xls=.07, bm=.27, file_name='CutMeans'))
    # endregion ANA
    # ----------------------------------------

    @staticmethod
    def invert(cut):
        cut = Cut.to_string(cut)
        return TCut(f'!({cut})' if cut else '')

    @staticmethod
    def to_string(cut: Any):
        return cut.GetTitle() if type(cut) is TCut else cut.Value if type(cut) is CutString else choose(cut, '')

    @staticmethod
    def make(name, cut, invert=False):
        return TCut(name, Cut.invert(cut) if invert else Cut.to_string(cut))

    @staticmethod
    def sum(*c, name='sum'):
        return sum([TCut(i) for i in c], start=TCut(name, ''))


class CutString:

    def __init__(self, name='', value='', description='', level=1):
        self.Name = name
        self.Value = Cut.to_string(value)
        self.Level = level
        self.Description = description

    def __call__(self, cut=None):
        return TCut(self.Name, self.Value) if cut is None else TCut(cut)

    def __str__(self):
        return self.Value

    def __repr__(self):
        return f'{self.Level:2d}: Cut {self.Description}'

    def __add__(self, other):
        if other is not None:
            self.Value = (self() + (other() if type(other) is CutString else TCut(other))).GetTitle()
        return self

    def reset(self):
        self.Value = ''

    def set(self, value):
        self.Value = Cut.to_string(value)

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

    def remove(self, name):
        if name not in self.Strings:
            return warning(f'there is no cut with name {name}', prnt=name is not None)
        self.Strings.pop(name)

    def sort(self):
        self.Strings = OrderedDict(sorted(self.Strings.items(), key=lambda x: x[1].Level))

    def get_names(self, with_raw=False):
        return (['raw'] if with_raw else []) + list(self.Strings)

    def get(self, name, invert=False, warn=True):
        if self.has_cut(name):
            return self.Strings[name].invert() if invert else self.Strings[name]()
        warning('There is no cut with the name "{name}"!'.format(name=name), prnt=warn)
        return TCut('')

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


def low_rate(func, high=False):
    @wraps(func)
    def wrapper(*args, **kwargs):
        ana, run = args[0], args[0].Run.get_high_rate_run(high)
        if run != ana.Run.Number:
            from pad.analysis import PadAnalysis
            return getattr(PadAnalysis(run, ana.DUT.Number, ana.TCString, prnt=False).Cut, func.__name__)(*args[1:], **kwargs)
        return func(*args, **kwargs)
    return wrapper


def high_rate(func):
    return low_rate(func, high=True)
