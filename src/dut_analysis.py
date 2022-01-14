#!/usr/bin/env python
# --------------------------------------------------------
#       parent class for the analysis of a single device under test
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import vectorize, meshgrid, digitize, histogram2d, lexsort, invert, any, unravel_index, argmax
from numpy.random import rand
from uncertainties.umath import sqrt as usqrt  # noqa

from src.currents import Currents
from src.cut import Cut
from src.run import Run
from src.tracks import Tracks
from src.telescope import Telescope
from plotting.save import *
from helpers.utils import *
from src.analysis import Analysis
from plotting.fit import Langau
from src.binning import Bins
from functools import wraps


def reload_tree(f):
    @wraps(f)
    def inner(ana, *args, **kwargs):
        ana.reload_tree_()
        return f(ana, *args, **kwargs)
    return inner


class DUTAnalysis(Analysis):
    def __init__(self, run_number, diamond_nr=1, testcampaign=None, load_tree=None, verbose=False, prnt=True):

        self.Run = self.init_run(run_number, testcampaign, load_tree, verbose)
        self.DUT = self.Run.DUTs[diamond_nr - 1]
        super(DUTAnalysis, self).__init__(testcampaign, sub_dir=join(self.DUT.Name, str(self.Run.Number)), verbose=verbose)

        self.Tree = self.Run.Tree
        self.StartTime = self.Run.StartTime if self.Tree.Hash() else time_stamp(self.Run.LogStart)

        self.print_start(run_number, prnt, dut=self.DUT.Name)
        self.update_config()

        # Sub-Analyses
        self.Currents = Currents(self)

        if self.Tree.Hash():
            self.NRocs = self.Run.NPlanes
            self.Cut = self.init_cut()
            self.StartEvent = self.Cut.get_min_event()
            self.EndEvent = self.Cut.get_max_event()
            self.Bins = Bins(self)
            self.Tracks = Tracks(self)
            self.Tel = Telescope(self)

            # alignment
            self.Alignment = None
            self.IsAligned = self.check_alignment()

    def __str__(self):
        return f'{self.__class__.__name__} of {self.Run} and {self.DUT}'

    def __repr__(self):
        return self.__str__()

    # ----------------------------------------
    # region INIT
    @staticmethod
    def init_run(run_number, testcampaign, load_tree, verbose):
        return Run(run_number, testcampaign, load_tree, verbose)

    def init_cut(self):
        return Cut(self)

    def update_config(self):
        pass

    @staticmethod
    def get_info_header():
        return ['Run', 'Type', 'Diamond', 'Flux [kHz/cm2]', 'HV [V]']

    def show_information(self, header=True, prnt=True, ret_row=False):
        rows = [[self.Run.Number, self.Run.Info['runtype'], self.DUT.Name, '{:14.1f}'.format(self.Run.Flux.n), '{:+6.0f}'.format(self.DUT.Bias)]]
        if ret_row:
            return rows[0]
        print_table(rows, self.get_info_header() if header else None, prnt=prnt)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region SAVE
    def save_all(self, print_link=True):
        self.save_data()
        self.save_plots(print_link)

    def save_data(self, data=None):
        if self.Draw.MountExists:
            data = choose(data, self.get_data())
            with h5py.File(join(self.Draw.ServerMountDir, 'data', 'data.hdf5'), 'a') as f:
                if self.TCString not in f:
                    f.create_group(self.TCString)
                if str(self.DUT.Number) not in f[self.TCString]:
                    f[self.TCString].create_dataset(str(self.DUT.Number), data=zeros((self.Run.get_max_run() + 1, len(data), 2), 'd'))
                f[self.TCString][str(self.DUT.Number)][self.Run.Number] = array([[v.n, v.s] for v in data], 'd')

    @reload_tree
    @quiet
    def save_plots(self, print_link=True):
        self.draw_hitmap(show=False)
        self.draw_signal_distribution(show=False)
        self.draw_signal_map(show=False)
        self.Currents.draw(show=False, fname='Current')
        self.draw_flux(save=self.has_branch('rate'), show=False)
        self.draw_pulse_height(show=False)
        self.Draw.print_http('plots.html', force_print=print_link)

    def save_tree(self, cut=None):
        f = TFile('test.root', 'RECREATE')
        t = self.Tree.CloneTree(0)
        n = self.Tree.Draw('Entry$', self.Cut(cut), 'goff')
        good_events = self.Run.get_tree_vec(n, dtype='i4')
        self.PBar.start(n)
        for i, ev in enumerate(good_events):
            self.Tree.GetEntry(ev)
            t.Fill()
            self.PBar.update(i)
        f.cd()
        t.Write()
        macro = self.Run.RootFile.Get('region_information')
        if macro:
            macro.Write()
        f.Write()
        f.Close()
        self.info('successfully saved tree with only cut events.')
    # endregion SAVE
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_data(self):
        return []

    def has_branch(self, branch):
        return self.Run.has_branch(branch)

    def get_time(self):
        return self.Run.get_time()

    def get_track_vars(self, mm=True, local=True, pixel=False):
        return self.Cut.get_track_vars(self.DUT.Number - 1, mm, local, pixel)

    def get_t_var(self):
        return 'time / 1000.' if self.Run.TimeOffset is None else '(time - {}) / 1000.'.format(self.Run.TimeOffset)

    def get_eff_var(self, *args, **kwargs):
        return ''

    def get_t_off(self, rel_time):
        return self.Run.StartTime if rel_time else 0

    def get_t_args(self, rel_time=False):
        return {'x_tit': 'Time [hh::mm]', 't_ax_off': self.get_t_off(rel_time)}

    def get_tree_vec(self, var, cut='', dtype=None, nentries=None, firstentry=0):
        return self.Run.get_tree_vec(var, cut, dtype, nentries, firstentry)

    def get_events(self, cut=None, redo=False):
        if type(cut) == str:
            return self.get_tree_vec('Entry$', cut, dtype='i4')
        cut = self.Cut(cut)
        return do_hdf5(self.make_simple_hdf5_path('', cut.GetName(), 'Events'), self.get_tree_vec, redo, dtype='i4', var='Entry$', cut=cut) if cut.GetTitle() else arange(self.Run.NEvents)

    def get_event_cut(self, cut=None, redo=False):
        return self.make_event_cut(self.get_events(cut, redo))

    def set_event_cut(self, events):
        """ selects the [events] in the TTree, undo with self.reset_entries()"""
        from ROOT import TEntryList
        elist = TEntryList()
        for ev in array(events, 'i'):
            elist.Enter(ev)
        self.Tree.SetEntryList(elist)

    def reset_entries(self):
        """ reset the TEntryList from the TTree """
        self.Tree.SetEntryList(0)

    def make_event_cut(self, events):
        c = zeros(self.Run.NEvents, dtype=bool)
        c[events] = True
        return c

    def get_pulser_cut(self, inv=False):
        cut = self.make_event_cut(self.get_events(self.Cut['pulser']))
        return invert(cut) if inv else cut

    def get_sub_events(self, cut):
        e = array(self.get_events())
        s = array(self.get_events(cut))
        cut = zeros(self.Run.NEvents + 1, bool)
        cut[s] = True
        return cut[e]

    def get_event_at_time(self, seconds, rel=False):
        return self.Run.get_event_at_time(seconds, rel)

    @save_pickle(sub_dir='Entries', suf_args='all')
    def get_n_entries(self, cut=None, _redo=False):
        return self.Tree.GetEntries(self.Cut(cut).GetTitle())

    def get_current(self):
        return self.Currents.get()

    def get_irradiation(self):
        return self.DUT.get_irradiation(self.TCString)

    def get_attenuator(self):
        return False

    def get_ph_var(self, ped=False):
        return ''

    def get_pulse_height(self, *args, **kwargs):
        return ufloat(0, 0)

    def get_sm_data(self, cut=None, fid=False):
        """ :return: signal map data as numpy array [[x], [y], [ph]] with units [[mm], [mm], [mV]]
            :param cut: applies all cuts if None is provided.
            :param fid: return only values within the fiducial region set in the AnalysisConfig.ini"""
        cut = self.Cut.generate_custom(exclude=['fiducial'], prnt=False) if not fid and cut is None else self.Cut(cut)
        return self.get_tree_vec(self.get_track_vars() + [self.get_ph_var()], cut)

    @save_pickle('Uniformity', sub_dir='Signal', suf_args='all')
    def get_uniformity(self, use_fwc=True, _redo=False):
        return self.draw_uniformity(use_fwc=use_fwc, redo=_redo, show=False)

    def get_track_length_var(self):
        dx2, dy2 = ['TMath::Power(TMath::Tan(TMath::DegToRad() * {}_{}), 2)'.format('slope' if self.Run.has_branch('slope_x') else 'angle', direction) for direction in ['x', 'y']]
        return '{} * TMath::Sqrt({} + {} + 1)'.format(self.DUT.Thickness, dx2, dy2)

    def get_flux(self, plane=None, corr=True, full_size=False, redo=False):
        return self.Tel.get_flux(plane, corr, True, full_size, _redo=redo)

    def get_ph_values(self, *args, **kwargs):
        """ :returns: all pulse height values for a given cut. [numpy.ndarray] """

    def get_signal_name(self, *args, **kwargs):
        """ :returns: the pulse height variable in the tree. [str] """

    def get_signal_var(self, *args, **kwargs):
        """ :returns: the pulse height variable in the tree + corrections. [str] """

    def get_split_ph(self, m=2):
        return get_2d_hist_vec(self.split_signal_map(m, show=0)[0])

    def get_next_dut(self, dut_nr=None):
        return next((dut for dut in self.Run.DUTs if dut.Number != self.DUT.Number), 1) if dut_nr is None else self.Run.DUTs[dut_nr - 1]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region ALIGNMENT
    def get_alignment(self):
        from src.event_alignment import EventAligment
        return EventAligment

    def init_alignment(self):
        if self.Alignment is None:
            self.Alignment = self.get_alignment()(self.Run.Converter)

    def get_aligned(self, *args, **kwargs):
        return self.get_alignment()(self.Run.Converter).get_aligned(*args, **kwargs)

    @save_pickle('Events', sub_dir='Alignment')
    def check_alignment_(self, _redo=False):
        return calc_eff(values=self.get_aligned())[0] > 99.7

    def check_alignment(self, redo=False):
        """ check if the events are aligned"""
        is_aligned = self.check_alignment_(_redo=redo)
        warning('Run {r} is misaligned :-('.format(r=self.Run.Number), prnt=not is_aligned)
        return is_aligned

    def draw_alignment(self, bin_size=1000, **kwargs):
        """ draw the aligment of telescope and DUT events """
        self.init_alignment()
        bins, y = self.Alignment.get_time_bins(bin_size=bin_size), self.Alignment.get_aligned(bin_size=bin_size)
        x, y = (bins[1][:-1] + diff(bins[1]) / 2)[:y.size].repeat(y + 1), full(sum(y + 1), 1)
        h = self.Draw.histo_2d(x, y, bins + [3, 0, 3], 'Event Alignment', **prep_kw(kwargs, x_tit='Time hh:mm', y_tit='Alignment', stats=False, l_off_y=99, center_y=True,
                               draw_opt='col', **Draw.mode(2, lm=.05, y_off=.3), pal=(3, array([1, 633, 418], 'i')), t_ax_off=0, rm=.03, z_range=[0, 2]))
        Draw.legend([Draw.box(0, 0, 0, 0, line_color=c, fillcolor=c) for c in [418, 633]], ['aligned', 'misaligned'], 'f')
        self.Draw.save_plots('EventAlignment')
        return h

    # endregion ALIGNMENT
    # ----------------------------------------

    # ----------------------------------------
    # region ALIASES
    def draw_chi2(self, *args, **kwargs):
        return self.Tracks.draw_chi2(*args, **kwargs)

    def draw_chi2s(self, *args, **kwargs):
        return self.Tracks.draw_chi2s(*args, **kwargs)

    def draw_angle(self, *args, **kwargs):
        return self.Tracks.draw_angle(*args, **kwargs)

    def draw_angles(self, *args, **kwargs):
        return self.Tracks.draw_angles(*args, **kwargs)

    def draw_occupancies(self, *args, **kwargs):
        return self.Tel.draw_occupancies(*args, **kwargs)

    def draw_cluster_occupancies(self, *args, **kwargs):
        return self.Tel.draw_cluster_occupancies(*args, **kwargs)

    def get_mean_angle(self, mode, redo):
        return self.Tracks.get_mean_angle(mode, redo)

    def draw_current(self, *args, **kwargs):
        return self.Currents.draw(*args, **kwargs)

    def draw_flux(self, *args, **kwargs):
        return self.Tel.draw_flux(*args, **kwargs)

    def get_pulse_height_trend(self, *args):
        return TProfile()

    def draw_pulse_height(self, *args, **kwargs):
        return TProfile()

    def draw_signal_distribution(self, *args, **kwargs):
        pass

    def draw_raw_signal_distribution(self, *args, **kwargs):
        return self.draw_signal_distribution(*args, **kwargs)

    def draw_fid_cut(self, scale=10):
        return self.Cut.draw_fid(scale)

    def get_efficiency(self, *args, **kwargs):
        pass

    def draw_beam_profile(self, *args, **kwargs):
        return self.Tracks.draw_beam_profile(*args, **kwargs)
    # endregion ALIASES
    # ----------------------------------------

    # ----------------------------------------
    # region SIZES
    def draw_size(self, size=None, color=1, name=''):
        if size is None:
            return warning('The {} size of "{}" was not specified in the dia info'.format(' {}'.format(name) if name else '', self.DUT.Name))
        c = get_last_canvas()
        cx, cy = self.find_center()
        x, y = size
        c.cd()
        Draw.box(cx - x / 2, cy - y / 2, cx + x / 2, cy + y / 2, line_color=color, width=3, name=name, fillstyle=4000)

    def draw_detector_size(self):
        self.draw_size(self.DUT.Size, color=432, name='detector')

    def draw_metal_size(self):
        size = [self.DUT.PadSize.n] * 2 if self.DUT.PadSize is not None else None
        self.draw_size(size, color=923, name='metal')

    def draw_guard_ring(self):
        size = [self.DUT.GuardRing] * 2 if self.DUT.GuardRing is not None else None
        self.draw_size(size, color=417, name='guard ring')

    def draw_all_sizes(self):
        self.draw_fid_cut()
        self.draw_detector_size()
        self.draw_metal_size()
        self.draw_guard_ring()

    @save_pickle('Center', sub_dir='Maps')
    def find_center(self, _redo=False, h=None, _no_save=False):
        h = choose(h, self.get_signal_map(cut='', _redo=_redo))
        px, py = h.ProjectionX(), h.ProjectionY()
        return array([mean([p.GetBinCenter(b) for b in [p.FindFirstBinAbove(p.GetMaximum() * .5), p.FindLastBinAbove(p.GetMaximum() * .5)]]) for p in [px, py]])
    # endregion SIZES
    # ----------------------------------------

    # ----------------------------------------
    # region EFFICIENCY
    def draw_efficiency(self, bin_size=None, show=True):
        self.Draw.efficiency(*self.get_tree_vec([self.get_t_var(), self.get_eff_var()], self.Cut()), self.Bins.get_time(bin_size), show=show, **self.get_t_args())

    def draw_efficiency_map(self, res=None, n_sigma=4, cut=None, show=True):
        cut_string = choose(self.Cut(cut) + self.Cut.get('tracks'), self.Cut.generate_custom(exclude=['fiducial']), decider=cut)
        e, x, y = self.get_tree_vec(var=[self.get_eff_var(n_sigma)] + self.get_track_vars(), cut=cut_string)
        p = self.Draw.prof2d(x, y, e * 100, Bins.get_global(res), 'Efficiency Map', x_tit='X [cm]', y_tit='Y [cm]', z_tit='Efficiency [%]', ncont=100, z_range=[0, 100], show=show)
        set_2d_ranges(p, dx=3, dy=3)
        self.draw_fid_cut(scale=10)
        self.draw_detector_size()
        self.Draw.save_plots('EffMap')
    # endregion SIZES
    # ----------------------------------------

    def draw_ph_pull(self, bin_size=None, fit=True, binning=None, **kwargs):
        p = self.draw_pulse_height(bin_size, show=False, save=False)[0]
        h = self.Draw.pull(p, binning, ret_h=True, title=f'Signal Bin{Bins.get_size(bin_size)} Distribution', **prep_kw(kwargs))
        format_statbox(h, all_stat=True, fit=fit)
        h.Fit('gaus', f'qs{"" if fit else 0}')
        self.Draw.save_plots(h.GetTitle().replace(" ", ""))
        return h

    def draw_track_length(self, show=True):
        h = TH1F('htd', 'Track Distance in Diamond', 200, self.DUT.Thickness, self.DUT.Thickness + 1)
        self.Tree.Draw('{}>>htd'.format(self.get_track_length_var()), 'n_tracks', 'goff')
        format_histo(h, x_tit='Distance [#mum]', y_tit='Entries', y_off=2, lw=2, stats=0, fill_color=Draw.FillColor, ndivx=405)
        self.Draw(h, show, lm=.16)
        return h

    def draw_signal_vs_angle(self, mode='x', bin_size=.1, show=True):
        p = TProfile('psa{}'.format(mode), 'Pulse Height vs. Angle in {}'.format(mode.title()), *self.Bins.get_angle(bin_size))
        self.Tree.Draw('{}:angle_{m}>>psa{m}'.format(self.get_signal_var(), m=mode), self.Cut(), 'goff')
        format_histo(p, x_tit='Track Angle [deg]', y_tit='Pulse Height [mV]', y_off=1.5)
        self.Draw(p, show=show, lm=.12, stats=set_entries())

    def draw_signal_vs_trigger_phase(self, dut=False, cut=None, show=True):
        x, y = self.get_tree_vec(['trigger_phase[{}]'.format([0, 1][dut]), self.get_signal_var()], self.Cut(cut))
        self.Draw.profile(x, y, make_bins(11), 'Signal vs. Trigger Phase', x_tit='Trigger Phase', y_tit='Pulse Height [mV]', show=show)

    # ----------------------------------------
    # region SIGNAL MAP
    @save_pickle('SMRange', sub_dir='Maps', suf_args='all')
    def find_sm_range(self, res=None, square=False, m=None, n=None, _redo=False):
        var, bins = self.get_track_vars() + [self.get_ph_var()],  Bins.get_global(res, square) if m is None else self.get_fid_bins(m, n)
        return ax_range(get_2d_hist_vec(self.Draw.prof2d(*self.get_tree_vec(var, self.Cut()), bins, show=False)), thresh=4)

    @save_pickle('SM', sub_dir='Maps', print_dur=True, suf_args='all')
    def get_signal_map(self, res=None, cut=None, fid=False, square=False, m=None, n=None, local=True, _redo=False):
        self.Tree.SetEstimate(self.Run.NEvents)
        var, bins = self.get_track_vars(local=local) + [self.get_ph_var()], Bins.get_global(res, square) if m is None else self.get_fid_bins(m, n)
        x, y, zz = self.get_tree_vec(var, self.Cut.generate_custom(exclude='fiducial', prnt=False) if not fid and cut is None else self.Cut(cut))
        tit, ztit = 'Pulse Height Map', 'Pulse Height [mV]'
        return self.Draw.prof2d(x, y, zz, bins, tit, **self.Tracks.ax_tits(), z_tit=ztit, z_range=self.find_sm_range(res, square, m, n, _redo=_redo), show=False, pal=53)

    def draw_signal_map(self, res=None, cut=None, fid=False, square=False, m=None, n=None, local=True, scale=False, redo=False, **kwargs):
        h = self.get_signal_map(res, cut, fid, square, m, n, local, _redo=redo)
        h.Scale(1 / self.get_pulse_height().n) if scale else do_nothing()
        rz = array([h.GetMinimum(), h.GetMaximum()]) * 1 / self.get_pulse_height().n if scale else None
        h = self.Draw.prof2d(h, **prep_kw(kwargs, centre=4, ncont=50, ndivy=510, ndivx=510, pal=53, z_tit='Relative Pulse Height' if scale else None, z_range=rz))
        if m is None:
            self.draw_fid_cut() if 'mirror' not in kwargs and 'rot' not in kwargs else do_nothing()
        self.Draw.save_plots('SignalMap2D', **prep_kw(kwargs, save=res is None and m is None and not scale))
        return h

    @save_pickle('HM', sub_dir='Maps', suf_args='all')
    def get_hitmap(self, res=None, cut='', _redo=False):
        x, y = self.get_tree_vec(self.get_track_vars(), self.Cut(cut))
        return self.Draw.histo_2d(x, y, Bins.get_global(res), 'Hit Map', x_tit='Track Position X [mm]', y_tit='Track Position Y [mm]', show=False)

    def draw_hitmap(self, res=None, cut='', redo=False, **kwargs):
        return self.Draw(self.get_hitmap(res, cut, _redo=redo), **prep_kw(kwargs, centre=4, title='DUT Hit Map'), leg=self.Cut.get_fid(), file_name='HitMap')

    def get_fid_bins(self, m, n):
        n = choose(n, m)
        x, y = self.Cut.load_fiducial() * 10
        x_bins = linspace(x[0], x[2], m + 1)
        y_bins = linspace(y[0], y[2], n + 1)
        return [m, x_bins, n, y_bins]

    def get_fid_bin_events(self, m, n, mi, ni):
        m, x, n, y = self.get_fid_bins(m, n)
        cut = self.Cut.generate_sub_fid('f{}{}{}{}'.format(m, n, mi, ni), *array([x[mi], x[mi + 1], y[ni], y[ni + 1]]) / 10)
        return self.get_sub_events(cut)

    def split_signal_map(self, m=2, n=None, redo=False, draw_n=False, grid=True, **kwargs):
        h = self.draw_signal_map(fid=True, m=m, n=n, redo=redo, **prep_kw(kwargs, stats=False))
        Draw.grid(*get_2d_bins(h, arr=True), width=2) if grid else do_nothing()
        self.Draw.save_plots('SplitSigMap')
        Draw.bin_numbers(h, draw_n)
        return h

    def draw_split_means(self, n=10):
        x = arange(1, n + 1).tolist()
        y = [mean_sigma(get_2d_hist_vec(h))[0] for h in [self.split_signal_map(i, i, show=False)[0] for i in x]]
        self.Draw.graph(x, y, title='Split Means', x_tit='Division', y_tit='Pulse Height [mV]', draw_opt='ap')

    def get_ph_bins(self, n=10, pmin=90, pmax=95, show=True):
        h = self.split_signal_map(n, n, show=show, grid=False)
        (x, y), v = get_2d_vecs(h)
        wx, wy = diff(x)[0] / 2, diff(y)[0] / 2
        points = array(meshgrid(x, y)).T[where((v >= pmin) & (v < pmax))]
        for ix, iy in points:
            Draw.box(ix - wx, iy - wy, ix + wx, iy + wy, line_color=840, width=4, show=show)
        return points, wx, wy

    def draw_ph_bin_disto(self, n=10, pmin=90, pmax=95, **dkw):
        ph, x, y = self.get_tree_vec(var=[self.get_ph_var()] + self.get_track_vars(), cut=self.Cut())
        points, wx, wy = self.get_ph_bins(n, pmin, pmax, show=False)
        cut = any([(x > ix - wx) & (x < ix + wx) & (y > iy - wy) & (x < iy + wy) for ix, iy in points], axis=0)
        return self.Draw.distribution(ph[cut], tit=f'Pulse Height of Areas in [{pmin}, {pmax}] mV', **prep_kw(dkw, x_tit='Pulse Height [mV]'))

    def draw_normal_distribution(self, m=20, n=30, show=True):
        ph, x, y = self.get_tree_vec(var=[self.get_ph_var()] + self.get_track_vars(), cut=self.Cut())
        ix, bx, iy, by = self.get_fid_bins(m, n)
        n = cumsum(histogram2d(x, y, [bx, by])[0].flatten().astype('i'))[:-1]  # get the number of events for each bin
        values = split(ph[lexsort((digitize(x, bx), digitize(y, by)))], n)  # bin x and y and sort then ph according to bins
        values = concatenate([lst / mean(lst) * mean(ph) for lst in values if lst.size > 2])  # normalise the values of each bin
        return self.Draw.distribution(values, self.Bins.get_pad_ph(Bins.find_width(values)), 'Signal Distribution Normalised by area mean', x_tit='Pulse Height [mV]', show=show)

    def draw_sig_map_disto(self, res=None, cut=None, fid=True, x_range=None, redo=False, normalise=False, ret_value=False, ph_bins=None, show=True, save=True):
        x = get_2d_hist_vec(self.draw_signal_map(res, cut, fid, redo=redo, show=False), err=False) / (self.get_pulse_height() if normalise else 1)
        x_range = choose(x_range, ax_range(x, fl=.1, fh=.1))
        h = self.Draw.distribution(x, Bins.make(*x_range, n=sqrt(x.size)), 'Signal Map Distribution', x_tit='Pulse Height [mV]', y_off=2, lm=.15, show=show, save=save)
        return mean_sigma(x) if ret_value else h

    def draw_split_map_disto(self, m, n=None, x_range=None, fit=True, normalise=False, show=True):
        x = get_2d_hist_vec(self.split_signal_map(m, n, show=False)[0]) / (self.get_pulse_height() if normalise else 1)
        h = self.Draw.distribution(x, Bins.make(*choose(x_range, ax_range(x, 0, .1, .4, thresh=3)), n=sqrt(x.size)), 'Signal Map Distribution',
                                   x_tit='{}Pulse Height{}'.format(*['Normalised ', ''] if normalise else ['', ' [mV]']), show=show)
        h.Fit('gaus', 'qs{}'.format('' if fit else '0'))
        format_statbox(h, all_stat=True, fit=fit)
        return mean_sigma(x)

    def get_sm_std(self, res=None, redo=False):
        return self.draw_sig_map_disto(show=False, res=res, redo=redo, normalise=True, ret_value=True, save=False)[1]

    def draw_sm_profile(self, mode='x', factor=1.5, cut=None, fid=False, hitmap=False, redo=False, show=True):
        s = self.draw_signal_map(factor, cut, fid, hitmap=hitmap, redo=redo, show=False)
        g = Draw.make_tgrapherrors('g_smp', 'Signal Map Profile')
        values = [[] for _ in range(s.GetNbinsX() if mode == 'x' else s.GetNbinsY())]
        for xbin in range(s.GetNbinsX()):
            for ybin in range(s.GetNbinsY()):
                value = s.GetBinContent(xbin, ybin)
                if value:
                    values[(xbin if mode == 'x' else ybin)].append(value)
        for i, lst in enumerate(values):
            m, sigma = mean_sigma(lst) if lst else (0., 0.)
            xval = s.GetXaxis().GetBinCenter(i) if mode == 'x' else s.GetYaxis().GetBinCenter(i)
            g.SetPoint(i, xval, m)
            g.SetPointError(i, 0, sigma)

        format_histo(g, x_tit='Track in {m} [cm]'.format(m=mode), y_tit='Pulse Height [au]', y_off=1.5, ndivx=515)
        self.Draw(g, 'SignalMapProfile', draw_opt='ap', lm=.14, show=show, gridx=True)

    def draw_error_signal_map(self, show=True):
        h = self.draw_signal_map(show=False, fid=True).ProjectionXY('', 'c=e')
        format_histo(h, name='hsme', title='Signal Map Errors', y_off=1.6)
        self.Draw(h, 'SignalMapErrors', lm=.12, rm=.11, show=show, draw_opt='colz')
        return h

    def get_signal_spread(self, min_percent=5, max_percent=99, prnt=True):
        """ Calculates the relative spread of mean signal response from the 2D signal response map. """
        pickle_path = self.make_pickle_path('SignalMaps', 'Spread', self.Run.Number, self.DUT.Number)

        def f():
            h = self.draw_sig_map_disto(show=False)
            q = array([min_percent, max_percent]) / 100.
            y = zeros(2, 'd')
            mean_error = mean([v.n for v in get_2d_hist_vec(self.draw_error_signal_map(show=False))])
            h.GetQuantiles(2, y, q)
            return ufloat(y[1], mean_error) / ufloat(y[0], mean_error) - 1
        ratio = do_pickle(pickle_path, f)
        self.info('Relative Signal Spread is: {:2.2f} %'.format(ratio * 100), prnt=prnt)
        return ratio

    def draw_low_sig_map(self, high=10, low=None, fid=False, cut=None, hit_map=True):
        low = '&&{}>{}'.format(self.get_signal_var(), low) if low is not None else ''
        fid_cut = self.Cut.generate_custom(exclude='fiducial') if cut is None else self.Cut(cut)
        kwargs = {'redo': True, 'cut': TCut('{}<{}{}'.format(self.get_signal_var(), high, low)) + (self.Cut(cut) if fid else fid_cut)}
        self.draw_hitmap(**kwargs) if hit_map else self.draw_signal_map(**kwargs)

    def correlate_sm(self, run, col=True, **kwargs):
        sm1, sm2 = [ana.get_signal_map(fid=True, **kwargs) for ana in [self, self.__class__(run, self.DUT.Number, self.TCString, verbose=False, prnt=False) if isint(run) else run]]
        c = correlate_all_maps(sm1, sm2)
        return c if col else c.T

    def find_best_sm_correlation(self, run):
        c = self.correlate_sm(run, res=1)
        self.info(f'best correlation: {max(c):.2f} at shift {unravel_index(argmax(c), c.shape)}')
        return c.max()

    def draw_sm_correlation_vs_shift(self, run, col=True, res=1, n=4, **kwargs):
        c = self.correlate_sm(run, col, res=res)
        iy, ix = unravel_index(argmax(c), c.shape)
        x = (arange(2 * n + 1) - n + ix)
        x = x % c.shape[1] if max(x) > c.shape[1] else x
        self.Draw.graph(x, c[iy][x], **prep_kw(kwargs, x_tit=f'Shift in {"X" if col else "Y"}', y_tit='Correlation Factor'))

    def draw_sm_correlation(self, run, m=10, show=True):
        x0, x1 = [get_2d_hist_vec(f.split_signal_map(m, show=False)[0], err=False) for f in [self, self.__class__(run, self.DUT.Number, self.TCString, prnt=False)]]
        g = self.Draw.histo_2d(x0, x1, self.Bins.get_pad_ph(2) * 2, 'Signal Map Correlation', show=show, x_tit='Pulse Height {} [mV]'.format(self.Run.Number), y_tit='Pulse Height {} [mV]'.format(run),
                               x_range=ax_range(x0, 0, .1, .1), y_range=ax_range(x1, 0, .1, .1))
        Draw.info('Correlation Factor: {:.2f}'.format(g.GetCorrelationFactor()))

    def draw_corr_coeff(self, run, show=True):
        x = [5, 10, 25, 50, 100, 200]
        ana = self.__class__(run, self.DUT.Number, self.TCString, prnt=False)
        y = [correlate(*[get_2d_hist_vec(f.split_signal_map(m, show=False)[0], err=False, zero_supp=0) for f in [self, ana]]) for m in x]
        self.Draw.graph(x, y, 'Signal Map Correlation', x_tit='Number of Bins', y_tit='Correlation Coefficient', show=show)
    # endregion SIGNAL MAP
    # ----------------------------------------

    def draw_uniformity(self, h=None, use_fwc=True, redo=False, **kwargs):
        noise = self.Pedestal.get_fwhm(raw=True, redo=redo) if h is None and hasattr(self, 'Pedestal') else 0
        h = choose(h, self.draw_raw_signal_distribution, redo=redo, **prep_kw(kwargs, normalise=True))
        format_histo(h, **prep_kw(kwargs, x_range=ax_range(h=h, thresh=h.GetMaximum() * .02, fl=.2, fh=.6)))
        format_statbox(h, w=.35, all_stat=True)
        (low, high), m = get_fwhm(h, ret_edges=True), get_fw_center(h) if use_fwc else ufloat(h.GetMean(), h.GetMeanError())
        fwhm, half_max = high - low, h.GetMaximum() / 2
        if noise > fwhm:
            return fwhm, noise, 0
        li = Draw.vertical_line(m.n, style=7, w=2)
        a = Draw.arrow(low.n, high.n, half_max, half_max, col=2, width=3, opt='<>', size=.02)
        Draw.legend([a, li], ['FWHM', 'FWC' if use_fwc else 'Mean'], 'l', y2=.72, w=.2)
        fwhm = usqrt(fwhm ** 2 - noise ** 2)  # correct fwhm for noise
        fwhm += ufloat(0, 1 if fwhm > 15 else -.3 * fwhm + 5.3)  # add error for false estimate
        value = fwhm / m
        Draw.add_stats_entry(h, f'FWHM/{"FWC" if use_fwc else "Mean"}', value, line=3)
        self.info(f'Uniformity: {value:.2f}')
        self.Draw.save_plots('Uniformity', **kwargs)
        return m, fwhm, value

    def model_trap_number(self, f=1000, t=1, max_traps=10000, steps=20, show=True):
        filled_traps = zeros(steps, dtype=int)
        decay = vectorize(self.decay)
        n_traps = []
        for i in range(steps):
            filled_traps[i] = f
            filled_traps = decay(filled_traps, t)
            n_traps.append(min(sum(filled_traps), max_traps))
        g = Draw.make_tgrapherrors(x=arange(steps), y=n_traps)
        format_histo(g, title='Number of Filled Traps', x_tit='Time [s]', y_tit='Number of Traps', y_off=1.7)
        self.Draw(g, draw_opt='ap', lm=.13, show=show)
        return n_traps[-1]

    def draw_n_traps(self, t=1, max_traps=1e5, nbins=20):
        x, y = log_bins(nbins, 100, 1e6)[1], []
        self.PBar.start(x.size)
        for f in x:
            y.append(self.model_trap_number(f, t, int(max_traps), show=False))
            self.PBar.update()
        g = Draw.make_tgrapherrors(x=x / 1000, y=y)
        format_histo(g, title='Number of Filled Traps vs Flux', x_tit='Flux [kHz/cm^{2}]', y_tit='Number of Filled Traps', y_off=1.7)
        self.Draw(g, draw_opt='ap', lm=.13, logx=True)

    # ----------------------------------------
    # region HELPERS
    def reload_tree_(self):
        self.Tree = self.Run.load_rootfile(prnt=False)
        for field in self.__dict__.values():
            if hasattr(field, 'Tree'):
                field.Tree = self.Tree

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8, fit_range=None):
        h = self.draw_signal_distribution(show=False) if h is None and hasattr(self, 'draw_signal_distribution') else h
        self.Draw.histo(h, show=show)
        fit = Langau(h, nconv, fit_range)
        fit.get_parameters()
        fit(show=show)
        update_canvas()
        if fit.get_chi2() > chi_thresh and nconv < 80:
            self.info('Chi2 too large ({c:2.2f}) -> increasing number of convolutions by 5'.format(c=fit.get_chi2()))
            fit = self.fit_langau(h, nconv + Draw.get_count('langau') * 5, chi_thresh=chi_thresh, show=show)
        print('MPV: {:1.1f}'.format(fit.get_mpv()))
        format_statbox(h, fit=True, all_stat=True)
        Draw.reset_count('langau')
        self.Draw.add(fit)
        return fit

    @staticmethod
    def decay(n, t):
        return count_nonzero(rand(n) <= exp(-1. / t))
    # endregion HELPERS
    # ----------------------------------------


if __name__ == '__main__':
    pargs = init_argparser(run=88, has_verbose=True, tree=True, dut=True)
    z = DUTAnalysis(pargs.run, pargs.dut, pargs.testcampaign, pargs.tree, verbose=True)
