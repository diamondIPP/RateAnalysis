#! /usr/bin/env python
# --------------------------------------------------------
#       DIAMOND RATE SCANS
# created on June 24th 2016 by M. Reichmann
# --------------------------------------------------------

from json import dump
from numpy import insert

import plotting.latex as latex
from analyse import collection_selector
from pad.collection import AnalysisCollection, PadCollection, fname
from pixel.collection import PixCollection
from plotting.draw import get_graph_y, ax_range, markers, TMultiGraph, mean_sigma, FitRes, set_statbox, make_ufloat
from src.analysis import *
from src.binning import Bins
from src.dut import PixelDUT
from src.run_selection import RunSelector
from src.voltage_scan import VoltageScan


class DiaScans(Analysis):

    Dir = Dir.joinpath(Analysis.MainConfig.get('SELECTION', 'dir'))
    Selections = load_json(Dir.joinpath(Analysis.MainConfig.get('SELECTION', 'runplan selection file')))

    def __init__(self, selection_name=None, verbose=False):

        self.Name = selection_name
        Analysis.__init__(self, verbose=verbose, results_dir='selections', sub_dir=selection_name if selection_name is not None else '')

        self.print_start(run=selection_name, tc=False, prnt=verbose)

        # Main
        self.Selection = self.load_selection(selection_name)

        # Config
        self.DUTParser = load_parser(join(self.Dir, 'config', 'DiamondAliases.ini'))

        # Info
        self.RS = RunSelector()  # dummy for information
        self.RunPlans = self.load_runplans()
        self.TestCampaigns = self.load_test_campaigns()
        self.Info = self.load_selection_info()
        self.DUTName = self.load_dut_name()
        self.NPlans = len(self.Info) if self.Info else None
        self.Ana = self.load_dummy()

        if self.Info:
            self.show_selection()
        self.print_finished(prnt=verbose)

    def __str__(self):
        return self.Name

    def __repr__(self):
        return f'{self.Name} ({self.DUTName})\n{self.show_selection(ret=True)}'

    @property
    def is_volt_scan(self):
        return all([i.is_volt_scan for i in self.Info])

    def insert_dummy(self, i):
        self.Info = insert(self.Info, i, array(self.Info)[i]).tolist()
        self.NPlans = len(self.Info)

    # ----------------------------------------
    # region RATE DEPENDENCE PARS
    def pol0_chi2s(self, redo=False):
        return array(self.get_values(self.Ana.pol0_chi2, PickleInfo('Pol0Chi2'), _redo=redo))

    def rel_spreads(self, redo=False):
        return array(self.get_values(self.Ana.rel_spread, PickleInfo('RelSpread'), _redo=redo), dtype=object) * 100

    def mean_rel_spread(self, **dkw):
        y = [v.flip_errors if hasattr(v, 'flip_errors') else v for v in self.rel_spreads()]  # asym errors are flipped for latex
        return FitRes(self.Draw.graph(arange(self.NPlans), y, **prep_kw(dkw, show=False, y_tit='Rel Spread [%]')).Fit('pol0', 'qs'))[0]

    def print_mean_rd(self, tex=False, prec=1):
        chi, s = [mean_sigma(lst)[0] for lst in [self.pol0_chi2s(), self.rel_spreads()]]
        if tex:
            print(latex.table(False, [[self.DUTName] + latex.si(self.get_bias_voltages()[0], fmt='.0f', unt='V') + latex.si(chi) + latex.si(s, unt='percent')]))
        else:
            info(f'Mean chi²:       {chi:.{prec}f}')
            info(f'Mean rel spread: {s:.{prec}f} %')

    def print_rate_dependencies(self):
        for i, s1, s2 in zip(self.Info, self.pol0_chi2s(), self.rel_spreads()):
            print(i)
            print(f'  Chi²:       {s1:.1f} %')
            print(f'  Rel spread: {s2:.1f} %')

    def print_rd_table(self):
        data = array([self.pol0_chi2s(), self.rel_spreads(), self.max_fluxes(avrg=True) / 1e3]).T.reshape((self.NPlans // 2, -1))
        rows = [latex.num(self.Info[2 * i].Irradiation.n, fmt='.1e', rm='+') + latex.num(s1, s2) + latex.si(r1, r2) + latex.num(f1.n, f2.n) for i, (s1, r1, f1, s2, r2, f2) in enumerate(data)]
        rows = [array(row)[[0, 1, 3, 5, 2, 4, 6]] for row in rows]
        print(latex.table(None, rows))
    # endregion RATE DEPENDENCE PARS
    # ----------------------------------------

    # ----------------------------------------
    # region INIT
    def load_selection(self, name):
        if name is None:
            return
        if name not in self.Selections.keys():
            warning('"{} does not exist in:\n{} '.format(name, sorted(self.Selections.keys())))
            self.Name = None
            return
        return self.Selections[name]

    def load_dut_name(self, name=None):
        name = choose(name, self.Name)
        if name is not None:
            if any([word in name.lower() for word in ['all', 'test']]):
                return name.title()
            return self.Info[0].DUTName

    def load_runplans(self):
        with open(self.RS.RunPlanPath) as f:
            return load(f, object_pairs_hook=OrderedDict)

    def load_test_campaigns(self):
        return sorted(self.RunPlans.keys() if self.Name is None else list(set(self.Selection.keys())))

    def load_selection_info(self):
        selections = []
        if self.Selection is None:
            return
        for tc, rps in self.Selection.items():
            for rp, dut_nrs in rps.items():
                for dut_nr in array([dut_nrs]).flatten():
                    selections.append(SelectionInfo(RunSelector(tc).select_runs_from_runplan(rp, dut_nr, unselect=True)))
        return selections

    def load_dummy(self):
        return AnalysisCollection if not self.Info else PixCollection if self.Info[0].IsPixel else PadCollection
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    @staticmethod
    def get_rp_diamonds(tc, rp):
        sel = RunSelector(tc).select_runs_from_runplan(rp, unselect=True)
        return sel.get_dut_names()

    def get_first_run(self, tc, rp):
        return self.get_runs(rp, tc)[0]

    def get_runs(self, rp, tc):
        return self.RunPlans[tc][rp]['runs']

    def get_dut_names(self):
        return [sel.DUTName for sel in self.Info]

    @property
    def run_types(self):
        return [sel.Type.lower().replace('rate', 'normal') for sel in self.Info]

    def get_irradiations(self, string=True):
        return array([irr2str(sel.Irradiation) if string else sel.Irradiation / 1e15 for sel in self.Info])

    def get_bias_voltages(self):
        return [sel.Bias for sel in self.Info]

    def get_bias_strs(self, root=True):
        return bias2str(self.get_bias_voltages(), root=root)

    def get_bias_str(self):
        return ' at {bias} V'.format(bias=self.Info[0].Bias) if len(set(self.get_bias_voltages())) == 1 else ''

    def get_biases(self):
        return [array([float(i.RunInfo[run][f'dia{i.DUTNr}hv']) for run in i.Runs]) for i in self.Info]

    def get_e_fields(self):
        return [array([i.DUT.get_e_field(float(i.RunInfo[run][f'dia{i.DUTNr}hv'])) for run in i.Runs]) for i in self.Info]

    def get_rp_values(self, sel, f, pickle_info=None, redo=False, load_tree=True, *args, **kwargs):
        pickle_path = pickle_info.path(sel) if pickle_info else ''
        if file_exists(pickle_path) and not redo:
            return do_pickle(pickle_path, do_nothing())
        self.info(f'Did not find {pickle_path}', prnt=pickle_path)
        with collection_selector(sel.RunPlan, sel.DUTNr, sel.TCString, load_tree) as ana:
            if 'redo' in signature(f).parameters:
                pf = partial(f, ana, redo=redo, *args, **kwargs)
                return do_pickle(pickle_path, pf, redo=redo) if pickle_info else pf()
            else:
                pf = partial(f, ana, *args, **kwargs)
                return do_pickle(pickle_path, pf, redo=redo) if pickle_info else pf()

    def get_values(self, f, pickle_info=None, redo=False, load_tree=True, *args, **kwargs):
        return [self.get_rp_values(sel, f, pickle_info, redo, load_tree, *args, **kwargs) for sel in self.Info]

    def get_pulse_heights(self, avrg=False, err=True, redo=False):
        return self.get_values(self.Ana.get_pulse_heights, PickleInfo('PHVals', avrg, err), redo=redo, avrg=avrg, err=err)

    def get_mean_ph(self, redo=False):
        return array([mean_sigma(i)[0] for i in self.get_pulse_heights(redo=redo)])

    def get_scaled_pulse_heights(self, avrg=False, redo=False):
        """ scale the ph to the mean of the pulse heights in the 'min flux range' """
        def scale_ph(x, y):
            xmin, xmax = self.MainConfig.get_list('MAIN', 'min flux range')
            y0 = y[where((x >= xmin) & (x <= xmax))]
            return y / mean(y0).n
        return [scale_ph(x, y) for x, y in zip(self.get_fluxes(avrg), self.get_pulse_heights(avrg, redo))]

    def get_rp_pulse_heights(self, sel, corr=True, redo=False):
        return self.get_rp_values(sel, self.Ana.get_pulse_heights, PickleInfo('PHVals', corr), redo=redo, corr=corr)

    def get_noise_spread(self, redo=False):
        noise = array(self.get_pedestals(redo), object)[:, 1]
        return mean_sigma([v - mean_sigma(lst)[0] for lst in noise for v in lst])[1]

    def get_mean_uniformities(self, use_fcw=True, redo=False, low=False, high=False):
        pickle_info = PickleInfo('Uni', low, high, use_fcw)
        return self.get_values(self.Ana.get_mean_uniformity, pickle_info, redo=redo, high_flux=high, low_flux=low)

    def get_uniformities(self, use_fcw=True, redo=False, low=False, high=False):
        pickle_info = PickleInfo('Uniformity', 'UniSMSTD', '{}{}{}'.format(int(low), int(high), int(use_fcw)))
        return self.get_values(self.Ana.get_uniformities, pickle_info, redo=redo, high_flux=high, low_flux=low, use_fcw=use_fcw)

    def get_currents(self):
        return self.get_values(self.Ana.get_currents, PickleInfo('CurrentsVals'))

    def get_fluxes(self, avrg=False, redo=False):
        return self.get_values(self.Ana.get_fluxes, PickleInfo('FluxVals', avrg), avrg=avrg, redo=redo)

    def max_fluxes(self, avrg=False, redo=False):
        return array([arr.max() for arr in self.get_fluxes(avrg, redo)])

    def times(self):
        return self.get_values(self.Ana.get_times, PickleInfo('Times'))

    def mean_times(self):
        t = array([[mean(t).n, t[0].n + t[0].s, t[-1].n + t[-1].s] for t in self.times()])
        return array([[it, it - s, e - it] for it, s, e in t])

    def get_x(self, avrg=False, e_field=False, irr=False, t=False, redo=False):
        if self.is_volt_scan:
            return self.get_e_fields() if e_field else self.get_biases()
        return self.times() if t else self.get_irradiations(string=False) if irr else self.get_fluxes(avrg, redo)

    def get_all_infos(self):
        return [sel for tc in self.RunPlans.keys() for sel in self.get_tc_infos(tc)]

    def get_dia_infos(self, dut_name):
        dut_name = self.RS.Run.translate_dia(dut_name)
        return [sel for tc in self.RunPlans.keys() for sel in self.get_tc_infos(tc) if dut_name == sel.DUTName]

    def get_tc_infos(self, tc):
        rs = RunSelector(tc)
        return [SelectionInfo(rs.select_runs_from_runplan(rp, dut + 1, unselect=True)) for rp in sorted(self.RunPlans[tc]) for dut in range(rs.get_n_duts(run_plan=rp))]

    def get_all_ana_strings(self, dut=None, tc=None, redo=False):
        selections = self.get_all_infos() if dut is None and tc is None else self.get_dia_infos(dut) if tc is None else self.get_tc_infos(tc)
        selections = [sel for sel in selections if self.RS.Run.translate_dia(dut) == sel.DUTName] if dut is not None else selections
        redo = ' -rd' if redo else ''
        return '_'.join(['analyse {} {} -c -tc {} -d{}'.format(sel.RunPlan, sel.DUTNr, sel.TCString, redo) for sel in selections])

    def get_savename(self, name):
        return '{}{}'.format(name, self.Name.title().replace('-', '').replace('_', ''))

    def get_signal_maps(self, fid=False, res=.2, square=True, scale=True, redo=False):
        pickle_info = PickleInfo('SM', make_suffix(self, fid, res, square, scale))
        return self.get_values(self.Ana.draw_signal_map, pickle_info, fid=fid, res=res, square=square, scale=scale, show=False, redo=redo)

    def get_cluster_size(self, avrg=False, redo=False):
        return self.get_values(PixCollection.get_cluster_sizes, PickleInfo('CS', avrg), redo=redo, avrg=avrg)

    def get_efficiency(self, avrg=False, redo=False):
        return self.get_values(PixCollection.get_efficiencies, PickleInfo('Eff', avrg), redo=redo, avrg=avrg)

    def get_x_args(self, vs_time=False, rel_time=False, vs_irrad=False, draw=True, e_field=False, **kwargs):
        return (VoltageScan if self.is_volt_scan else self.Ana).get_x_args(vs_time, rel_time, vs_irrad, draw, e_field=e_field, **kwargs)

    def flux_splits(self, redo=False):
        return self.get_values(self.Ana.get_flux_splits, PickleInfo('FluxSplit'), redo=redo)

    def flux_avrg(self, x):
        x = [ix[f.argsort()] for ix, f in zip(x, self.get_fluxes())]  # sort by ascending fluxes
        return [array([mean_sigma(lst)[0] for lst in split(ix, s)]) for ix, s in zip(x, self.flux_splits())]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region LEGEND
    def leg_titles(self, irr=None, dut=None, tc=None, bias=True):
        biases = lambda x: '' if self.is_volt_scan or not bias or len(set(self.get_bias_voltages())) == 1 else bias2rootstr(x.Bias)
        irrads = lambda x: irr2str(x.Irradiation) if irr else ''
        dut = choose(dut, len(set(self.get_dut_names())) > 1)
        duts = lambda x: x.DUT.full_name(x.TCString) if dut else ''
        tcs = lambda x: tc2str(x.TCString, short=False) if tc else ''
        return array([w for i in self.Info for w in [duts(i), biases(i), tcs(i), irrads(i)] if w])

    def make_legend(self, g, dut=None, tc=None, irr=False, bias=True, custom=False, **kw):
        tits = self.leg_titles(irr, dut, tc, bias)
        cols = len(tits) // len(g)
        styles = alternate(['p'] * len(g), zeros((cols - 1, len(g)), 'S'))
        if custom:
            return self.Draw.legend(g, custom, show=False, **prep_kw(kw, styles='p'))
        return self.Draw.legend(alternate(g, zeros((cols - 1, len(g)), 'i')), tits, show=False, **prep_kw(kw, scale=1, cols=cols, w=(.2 if dut else .25) + .15 * (cols - 1), styles=styles))

    def make_info_legend(self, i, g, irr, c=None, bias=True):
        tits = self.leg_titles(irr, bias=bias)
        n = tits.size // self.NPlans - 1
        w = len(max(tits, key=len)) * .011 * (n + 1)
        Draw.legend([g] + n * [''], tits.reshape((-1, n + 1))[i], ['pe'] + n * [''], w=w, ts=.2, scale=5, c=c, cols=n + 1)

    def make_full_legend(self, graphs, irr=True):
        same_bias = len(set(self.get_bias_voltages())) == 1
        cols = 1 + (not same_bias) + irr
        legend = Draw.make_legend(y2=.4, w=.12 * cols, nentries=4, cols=cols)
        for i, (g, sel) in enumerate(zip(graphs, self.Info)):
            legend.AddEntry(g, '{} - {}'.format(sel.RunPlan, make_tc_str(sel.TCString, long_=False)), 'lp')
            legend.AddEntry(0, irr2str(sel.Irradiation), '') if irr else do_nothing()
            legend.AddEntry(0, bias2rootstr(sel.Bias), '') if not same_bias else do_nothing()
        return legend
    # endregion LEGEND
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW
    def show_selections(self, dut=None):
        header = ['Name', 'Diamond', 'Campaigns']
        rows = []
        old_sel = deepcopy(self.Name)
        for name in [key for key in self.Selections if dut is None or dut.lower() in key.lower()]:
            self.set_selection_name(name)
            row = [name, self.load_dut_name(), ', '.join(str(tc) for tc in self.load_test_campaigns())]
            rows.append(row)
        self.set_selection_name(old_sel)
        print_table(rows, header)

    def show_selection(self, ret=False):
        """ Gives detailed information about the chosen selection """
        t = print_table(rows=[sel() for sel in self.Info], header=['TC', 'RunPlan', 'DUT', 'Nr', 'Runs', 'Bias', 'Type', 'Irrad'], prnt=not ret) if self.Info else warning('Selection is empty!')
        return t if ret else None

    def show_all_runplans(self):
        old_sel = deepcopy(self.Name)
        self.set_selection(None)
        for tc, runplans in sorted(self.RunPlans.items()):
            if runplans:
                print_small_banner(tc, color='yellow')
                header = ['Runplan', 'Runs', 'Diamonds']
                rows = []
                for rp, dic in sorted(runplans.items()):
                    runs = dic['runs']
                    rows.append([rp, '{:03d}-{:03d}'.format(runs[0], runs[-1]), ', '.join(self.get_rp_diamonds(tc, rp))])
                print_table(rows, header)
        self.set_selection_name(old_sel)

    def show_pulse_heights(self):
        rows = [[sel.TCString, sel.RunPlan] + ['{:2.2f}'.format(i) for i in mean_sigma([dic['ph'] for dic in phs.values()])] for phs, sel in zip(self.get_pulse_heights(), self.Info)]
        print_table(rows, header=['Test Campaign', 'Run Plan', 'PH [mv]', 'STD [mV]'])
    # endregion SHOW
    # ----------------------------------------

    # ----------------------------------------
    # region SELECTION
    def set_selection(self, name):
        self.info('Set Selection {0}'.format(name))
        self.Name = name
        self.Selection = self.load_selection(name)
        self.TestCampaigns = self.load_test_campaigns()
        self.Info = self.load_selection_info()
        self.DUTName = self.load_dut_name()
        self.Draw.set_results_dir(join('Results', 'selections', self.Name)) if name else do_nothing()

    def set_selection_name(self, name):
        self.Name = name
        self.Selection = self.load_selection(name)

    def clear_selection(self):
        self.Selection = {}

    def select_runplan(self, runplan, dut=1, testcampaign=None):
        rp = rp2str(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        if rp in self.RunPlans[tc]:
            if tc not in self.Selection:
                self.Selection[tc] = {}
            self.Selection[tc][rp] = dut
        else:
            warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def unselect_runplan(self, runplan, testcampaign=None):
        rp = rp2str(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        try:
            self.Selection[tc].pop(rp)
        except KeyError:
            warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def add_selection(self):
        name = input('Enter the name of the selection: ')
        self.clear_selection()
        self.Selections[name] = {} if name not in self.Selections else self.Selections[name]
        run_plan = input('Enter test campaign, run plan number, channel: ')
        while run_plan:
            tc, rp, ch = [string.strip(' ') for string in run_plan.split(',')]
            self.select_runplan(rp, int(ch), tc)
            run_plan = input('Enter test campaign, run plan number, channel (leave blank to finish): ')
        self.save_selection(name)

    def save_selection(self, name=None):
        name = input('Enter the name of the selection: ') if name is None else name
        if not self.Selection:
            warning('Selection is empty!')
            return
        with open(DiaScans.Dir.joinpath(self.MainConfig.get('SELECTION', 'runplan selection file')), 'w') as f:
            if name in self.Selections:
                query = input('{} does already exist. Do you want to overwrite it? (y/n) '.format(name))
                if query.lower().strip() in ['no', 'n']:
                    return
            self.Selections[name] = self.Selection
            dump(self.Selection, f, indent=2, sort_keys=True)
            self.info('Saved {} to selections'.format(name))
    # endregion SELECTION
    # ----------------------------------------

    # ----------------------------------------
    # region DRAWING
    def ph_spreads(self, avrg=True, rel=True):
        x = self.get_scaled_pulse_heights(avrg) if rel else self.get_pulse_heights(avrg)
        return array([max(i) - min(i) for i in x])

    def draw_pulse_heights(self, avrg=False, ef_ax=False, redo=False, **dkw):
        g = [self.Draw.graph(x, y, title='PH', y_tit=self.Ana.PhTit) for x, y in zip(self.get_x(avrg, redo=redo), self.get_pulse_heights(avrg, redo))]
        mg = self.Draw.multigraph(g, 'Pulse Heights', leg=self.make_legend(g, **dkw), **prep_kw(dkw, **self.get_x_args(draw=True), draw_opt='pl', tm=.116 if ef_ax else None))
        self.draw_ef_axis(mg, ef_ax)
        self.Draw.save_plots(fname('PH', avrg))

    def draw_scaled_pulse_heights(self, avrg=False, redo=False, yr=None, **dkw):
        x, y = self.get_x(avrg, redo=redo), self.get_scaled_pulse_heights(avrg, redo)
        g = [self.Draw.graph(ix, iy, title='PH', y_tit='Relative Pulse Height', marker=markers(i), show=False) for i, (ix, iy) in enumerate(zip(x, y))]
        f, yr = fname('NormalPH', avrg), None if yr is None else [1 - yr, 1 + yr]
        return self.Draw.multigraph(g, 'Scaled Pulse Heights', leg=self.make_legend(g, **dkw), **prep_kw(dkw, **self.get_x_args(draw=True), gridy=True, draw_opt='pl', file_name=f, y_range=yr))

    def draw_dia_rate_scans(self, redo=False, irr=True, corr=True):
        mg = TMultiGraph('mg_ph', '{dia} Rate Scans{b};Flux [kHz/cm^{{2}}]; Pulse Height [mV]'.format(dia=self.DUTName, b=self.get_bias_str()))
        mgs = self.get_values(self.Ana.draw_pulse_heights, PickleInfo('PHMG', '10000_{}'.format(corr)), redo=redo, show=False)
        for i, (mgi, sel) in enumerate(zip(mgs, self.Info)):
            for g in mgi.GetListOfGraphs():
                format_histo(g, color=self.Draw.get_color(self.NPlans, i), markersize=1.5, lw=2)
                if g.GetName() == 'gFirst':
                    format_histo(g, color=1, marker=26, markersize=2)
                elif g.GetName() == 'gLast':
                    format_histo(g, color=1, marker=23, markersize=2)
            mg.Add(mgi)
        legend = self.make_full_legend([mgi.GetListOfGraphs()[0] for mgi in mgs], irr)
        y = concatenate([get_graph_y(g) for g in mg.GetListOfGraphs()])
        format_histo(mg, draw_first=True, y_tit='Pulse Height [au]', y_range=[0, max(y).n * 1.1], tit_size=.05, lab_size=.05, y_off=.91, x_off=1.2, x_range=Bins.FluxRange)
        self.Draw(mg, 'DiaScans{dia}'.format(dia=make_dia_str(self.DUTName)), draw_opt='ap', logx=True, leg=legend, w=1.6, lm=.092, bm=.12, gridy=True)

    def draw_title_pad(self, h, x0, lm, c_height):
        if Draw.Title:
            biases = list(set(self.get_bias_voltages()))
            bias_str = ' at {b}'.format(b=bias2str(biases[0])) if len(biases) == 1 else ''
            Draw.tpad('p0', pos=[x0, 1 - h / c_height, 1, 1], margins=[0, 0, 0, 0], transparent=True)
            Draw.tpavetext('{dia} Rate Scans{b}'.format(dia=self.DUTName, b=bias_str), lm, 1, 0, 1, font=62, align=13, size=.5, margin=0)
            get_last_canvas().cd()

    def draw_scaled_rate_scans(self, irr=False, yr=.07, pad_height=.18, avrg=False, **dkw):
        both_pol = len(set(self.get_bias_voltages())) == 2
        title_height = pad_height / 2 if Draw.Title else .03                                # half of a pad for title
        c_height = (self.NPlans / (2 if both_pol else 1) + .5) * pad_height + title_height  # half of a pad for the x-axis
        c_width = 1.3 * pad_height / .2                                                     # keep aspect ratio for standard pad_height
        c = Draw.canvas(w=c_width, h=c_height, transp=True, logx=True, gridy=True)
        lm, rm, x0, size = .07, .02, .08, .22
        self.draw_title_pad(title_height, x0, lm, c_height)
        Draw.tpad('p1', pos=[0, 0, x0, 1], margins=[0, 0, 0, 0], transparent=True)           # info pad
        Draw.tpavetext('Scaled Pulse Height', 0, 1, 0, 1, align=22, size=.5, angle=90, margin=0)   # y-axis title

        g = self.get_values(self.Ana.draw_scaled_pulse_heights, PickleInfo('PHScaledGraph', avrg), avrg=avrg, irr=False, show=False)
        g = array(g).reshape((-1, 2)) if both_pol else g
        for i, ig in enumerate(g):
            c.cd()
            y0, y1 = [(c_height - title_height - pad_height * (i + j)) / c_height for j in [1, 0]]
            p = Draw.tpad(pos=[x0, y0, 1, y1], margins=[lm, rm, 0, 0], logx=True, gridy=True, gridx=True, fix=True)
            i, ig = i * (2 if both_pol else 1), make_list(ig)
            [format_histo(jg, color=self.Draw.get_color(self.NPlans), marker=markers(i + j, both_pol), markersize=1.5) for j, jg in enumerate(ig)]
            self.Draw.multigraph(ig, canvas=p, **prep_kw(dkw, **self.get_x_args(), show=False, color=None, gridy=True, lab_size=size, ndivy=504, x_ticks=.15, y_range=[1 - yr, 1 + yr]))
            self.make_info_legend(i, ig[0], irr, p.cd(), bias=False)
            c.cd()

        Draw.tpad('p2', pos=[x0, 0, 1, pad_height / 2 / c_height], margins=[lm, rm, 0, 0], transparent=True)  # x-axis pad
        Draw.x_axis(1, lm, 1 - rm, 'Flux [kHz/cm^{2}]', Bins.FluxRange, opt='', log=True, tick_size=0, lab_size=size * 2, tit_size=size * 2, off=1.1)
        self.Draw.save_plots(fname('ScaledPH', avrg))

    def draw_scaled_distribution(self, excluded=None):
        values = concatenate(([vals / mean_sigma(vals)[0] for i, vals in enumerate(self.get_pulse_heights()) if i not in make_list(excluded)]))
        self.Draw.distribution(values, [40, .9, 1.1], 'Scaled Pulse Height Distribution', x_tit='Scaled Pulse Height')
        return values

    def draw_currents(self, align=False, show=True):
        mg = TMultiGraph('mgc', 'Leakage Current vs. Flux')
        legend = Draw.make_legend(nentries=self.NPlans, w=.4, x2=.52)
        currents = self.get_currents()
        fluxes = self.get_fluxes()
        for i, (x, y) in enumerate(zip(fluxes, currents)):
            g = Draw.make_tgraph(x, y)
            if align:
                fit = g.Fit('pol1', 'qs0')
                g = Draw.make_tgraph(fluxes[i], array(currents[i]) - fit.Parameter(0) + .1)
            format_histo(g, color=self.Draw.get_color(self.NPlans))
            legend.AddEntry(g, '{tc} - {hv}'.format(tc=self.Info[i].TCString, hv=self.get_rp_values(self.Info[i], self.Ana.get_hv_name, load_tree=False)), 'pl')
            mg.Add(g)
        format_histo(mg, draw_first=True, y_tit='Current [nA]', x_tit='Flux [kHz/cm^{2}]', y_range=[.1, max(concatenate(currents)).n * 2], x_range=Bins.FluxRange)
        self.Draw(mg, 'CurrentFlux{}'.format(self.Name), draw_opt='ap', logx=True, logy=True, leg=legend, bm=.17, show=show)

    def set_bin_labels(self, h, rp=True):
        for i, sel in enumerate(self.Info):
            h.GetXaxis().SetBinLabel(h.GetXaxis().FindBin(i), f'{make_tc_str(sel.TCString, 0)}{f" - {sel.RunPlan}" if rp else ""}')

    def draw_means(self, **dkw):
        x, y = self.mean_times(), array([ufloat(*mean_sigma(ph_list, err=False)) for ph_list in self.get_pulse_heights()])
        x_args = self.get_x_args(vs_time=True, x_tit='Time [YY/MM]', tform='%y/%m')
        return self.Draw.graph(x, y, **prep_kw(dkw, title='MeanPH', **x_args, **Draw.mode(2), gridy=True, y_tit='Mean Pulse Height [mV]', file_name='MeanPH'))

    def draw_sigmas(self, y_range=None, show=True):
        y = array([mean_sigma(ph_list)[1] for ph_list in self.get_pulse_heights()])
        g = Draw.make_tgraph(arange(y.size), y=y, title='Pulse Height STD Evolution', x_tit='Run Plan', y_tit='Pulse Height Standard Deviation [%]')
        format_histo(g, y_off=1.2, x_range=ax_range(0, y.size - 1, .3, .3), x_off=2.5, y_range=y_range)
        self.set_bin_labels(g)
        self.Draw(g, self.get_savename('Sigmas'), show, draw_opt='ap', bm=.2, w=1.5, h=.75, gridy=True)

    def draw_uniformity(self, arg=2, use_fcw=True, redo=False, low=False, high=False):
        ph_var = 'FWC' if use_fcw else 'Mean'
        var, unit, tit = [ph_var, 'FWHM', 'FWHM/{}'.format(ph_var)][arg], ' [mV]' if arg < 2 else '', [ph_var, 'Full Width Half Maximum', 'Uniformity'][arg]
        y_values = [v[arg][0] for v in self.get_mean_uniformities(redo=redo, low=low, high=high)]
        x_values = array([ufloat(v, v * .2) for v in self.get_irradiations(string=False)])
        return self.Draw.graph(x_values, y_values, title=tit, x_tit='Irradiation [10^{15}n/cm^{2}]', y_tit='{}{}'.format(var, unit), draw_opt='ap', x_off=1.2)

    def draw_peak_flux(self, show=True):
        leg = Draw.make_legend(x2=.45, w=.3)
        mg = TMultiGraph('mgph', 'Peak Flux vs. FAST-OR Flux')
        values_list = self.get_values(PadCollection.get_peak_fluxes, PickleInfo('PeakFlux'))
        flux_list = self.get_fluxes()
        for sel, values, fluxes in zip(self.Info, values_list, flux_list):
            g = Draw.make_tgraph('g{}'.format(sel.RunPlan), '', x=fluxes, y=values)
            format_histo(g, color=self.Draw.get_color())
            leg.AddEntry(g, '{} @ {:+1.0f}V'.format(sel.DUTName, sel.Bias), 'pl')
            mg.Add(g, 'p')
        x, y = concatenate(flux_list), concatenate(values_list)
        x_range, y_range = [.5 * min(x).n, 1.2 * max(x).n], [.5 * min(y).n, 1.2 * max(y).n]
        format_histo(mg, draw_first=True, y_tit='Peak Flux [kHz/cm^{2}] ', x_tit='FAST-OR Flux [kHz/cm^{2}]', x_range=x_range, y_off=1.8, y_range=y_range)
        self.Draw(mg, 'PeakFluxes{}'.format(self.Name), draw_opt='a', leg=leg, show=show, lm=.13)

    def draw_purity(self, avrg=False, redo=False, **dkw):
        g = self.get_values(PadCollection.draw_purity, PickleInfo('Purity'), redo=redo, avrg=avrg, show=False)
        mg = self.Draw.multigraph(g, 'Purity', show=False)
        self.Draw(mg, **prep_kw(dkw, leg=self.make_legend(g, **dkw), **self.get_x_args(), y_range=[0, 105], file_name='Purity'))
    # endregion DRAWING
    # ----------------------------------------

    # ----------------------------------------
    # region PAD
    def draw_bucket_ratios(self, fit=False, avrg=False, redo=False, **dkw):
        g = self.get_values(PadCollection.draw_bucket_ratio, PickleInfo('BucketRatio', avrg, fit), redo=redo, avrg=avrg, show=False, stats=False, fit=fit, all_cuts=False)
        self.Draw.multigraph(g, 'Bucket Ratios', **prep_kw(dkw, **self.Ana.get_x_args(draw=True), wleg=.3))
        self.make_legend(g, dut=True, tc=True, **dkw).Draw('same')
        self.Draw.save_plots('BucketRatios')

    def draw_tp_ratios(self, redo=False, **dkw):
        g = self.get_values(PadCollection.draw_avrg_tp_ratio, PickleInfo('TPRatio'), redo=redo, show=False, fit=False)
        self.Draw.multigraph(g, 'TPRatios', **prep_kw(dkw, **self.Ana.get_x_args(draw=True)))
        self.make_legend(g, dut=True, tc=True, **dkw).Draw('same')
        self.Draw.save_plots('TPRatios')

    def snr(self, redo=False):
        return array(self.get_values(PadCollection.get_snr, PickleInfo('SNR'), redo=redo))[:, 0]

    def draw_snr(self, redo=False, **dkw):
        x, y = self.get_x(irr=True), self.snr(redo)
        return self.Draw.graph(x, y[:, 0], **prep_kw(dkw, **self.get_x_args(vs_irrad=True), y_tit='SNR', file_name='SNR'))

    def draw_sat_diff(self, redo=False, **dkw):
        g = self.get_values(PadCollection.draw_sat_ph_diff, PickleInfo('SatDiff'), redo=redo, show=False)
        mg = self.Draw.multigraph(g, 'Sat PH Diff')
        self.Draw(mg, **prep_kw(dkw, leg=self.make_legend(g, **dkw), **self.get_x_args(), file_name='SatDiff'))
    # endregion PAD
    # ----------------------------------------

    # ----------------------------------------
    # region PIXEL
    def draw_ef_axis(self, mg, draw=True, **dkw):
        if draw:
            x, y = array([getattr(mg.GetXaxis(), f'GetX{i}')() for i in ['min', 'max']]), get_last_canvas().GetUymax()
            Draw.x_axis(y, *x, 'Electric Field [V/#mum]', self.Info[0].DUT.get_e_field(x), **prep_kw(dkw, off=1.1, tit_size=.05, lab_size=.045, opt='-'))

    def draw_efficiency(self, avrg=False, e_field=False, ef_ax=False, redo=False, **dkw):
        g = [self.Draw.graph(x, y, 'Efficiency', y_tit='Hit Efficiency [%]', marker=markers(i), show=False) for i, (x, y) in enumerate(zip(self.get_x(avrg, e_field), self.get_efficiency(avrg, redo)))]
        mg = self.Draw.multigraph(g, 'Eff', leg=self.make_legend(g, **dkw), **prep_kw(dkw, **self.get_x_args(draw=True, e_field=e_field), draw_opt='pl', y_range=[0, 105], tm=.116 if ef_ax else None))
        self.draw_ef_axis(mg, ef_ax)
        self.Draw.save_plots(fname('Efficiency', avrg))

    def draw_cluster_size(self, avrg=False, e_field=False, dut=False, redo=False, **dkw):
        g = [self.Draw.graph(x, y[:, 0], title='Cluster Sizes', y_tit='Cluster Size') for x, y in zip(self.get_x(avrg, e_field), self.get_cluster_size(avrg, redo))]
        leg = self.make_legend(g, dut=dut, scale=1, **dkw)
        return self.Draw.multigraph(g, 'Cluster Sizes', leg=leg, **prep_kw(dkw, **self.get_x_args(draw=True, e_field=e_field), file_name='ClusterSize', draw_opt='pl'))
    # endregion PIXEL
    # ----------------------------------------

    # ----------------------------------------
    # region PULSER
    def get_pulser_pulse_heights(self, avrg=False, redo=False):
        return self.get_values(PadCollection.get_pulser_pulse_heights, PickleInfo('PulserPH', f'{avrg:d}'), redo=redo, avrg=avrg)

    def pulser_rates(self, avrg=False, redo=False):
        return self.get_values(PadCollection.get_pulser_rates, PickleInfo('PulserRate', f'{avrg:d}'), redo=redo, avrg=avrg)

    def get_pulser_stability(self):
        """ returns: relative standard deviation of the pulser """
        d = [s / m * 100 for m, s in [mean_sigma(v) for v in self.get_pulser_pulse_heights()]]
        for i, v in enumerate(d):
            print(f'Pulser stability for {self.Info[i]}: {v:.3f}')
        return d

    def draw_pulser_pulse_heights(self, avrg=False, scaled=True, ym=.01, redo=False, **dkw):
        data = zip(self.get_fluxes(avrg), self.get_pulser_pulse_heights(avrg, redo))
        g = [self.Draw.graph(x, y / (mean_sigma(y)[0].n if scaled else 1), y_tit='Pulser Pulse Height [mV]', markersize=1.2, show=False) for x, y in data]
        leg_tit = [f'{i.PulserType}al @ {bias2str(i.Bias)}' for i in self.Info]
        return self.Draw.multigraph(g, 'Pulser PH', leg_tit, **prep_kw(dkw, **self.Ana.get_x_args(draw=True), wleg=.3, y_range=[1 - ym, 1 + ym]), file_name='PulserPH')

    def draw_pulser_rates(self, **dkw):
        x, y = concatenate(self.get_fluxes()), concatenate(self.pulser_rates())
        self.Draw.graph(x, y, 'Pulser Rates', **prep_kw(dkw, y_tit='Pulser Fraction [%]', **self.get_x_args(), **Draw.mode(2), file_name='PulserRates', y_range=[0, max(y[:, 0]) * 1.2]))
    # endregion PULSER
    # ----------------------------------------

    # ----------------------------------------
    # region PEDESTAL
    def get_pedestals(self, avrg=True, redo=False):
        return self.get_values(PadCollection.get_pedestals, PickleInfo('PedVals'), redo=redo, avrg=avrg)

    def get_noise(self, avrg=True, redo=False):
        return self.get_values(PadCollection.get_noise, PickleInfo('NoiseVals'), redo=redo, avrg=avrg)

    def get_ped_spreads(self, avrg=True, redo=False):
        x = self.get_pedestals(avrg, redo)
        return array([max(i) - min(i) for i in x])

    def get_noise_spreads(self, avrg=True, redo=False):
        x = self.get_noise(avrg, redo)
        return array([max(i) - min(i) for i in x])

    def draw_ped_spreads(self, avrg=True, redo=False, **dkw):
        x, y = self.get_x(avrg, irr=True), self.get_ped_spreads(avrg, redo)
        return self.Draw.graph(x, y, 'Ped Spreads', **prep_kw(dkw, **self.get_x_args(vs_irrad=True, draw=True), y_tit='Pedestal Spread [mV]', file_name='PedSpread'))

    def draw_rel_ped_change(self, avrg=True, redo=False, **dkw):
        x, y, s = self.get_x(avrg, irr=True), self.get_ped_spreads(avrg, redo), self.get_mean_ph(redo)
        return self.Draw.graph(x, y / s * 100, 'Rel Ped Spread', **prep_kw(dkw, **self.get_x_args(vs_irrad=True, draw=True), y_tit='Pedestal Spread / Mean Pulse Height [%]', file_name='RelPedSpread'))

    def draw_noise_spreads(self, avrg=True, redo=False, **dkw):
        x, y = self.get_x(avrg, irr=True), self.get_noise_spreads(avrg, redo)
        return self.Draw.graph(x, y, 'Noise Spreads', **prep_kw(dkw, **self.get_x_args(vs_irrad=True, draw=True), y_tit='Noise Spread [mV]', file_name='NoiseSpread'))

    def draw_pedestals(self, rel=False, redo=False, show=True, irr=True):
        mg = TMultiGraph('mg_ph', '{dia} Pedestals{b};Flux [kHz/cm^{{2}}]; Pulse Height [mV]'.format(dia=self.DUTName, b=self.get_bias_str()))
        for i, (values, sel, fluxes) in enumerate(zip(self.get_pedestals(redo), self.Info, self.get_fluxes())):
            pedestals = array([make_ufloat(*tup) for tup in array(values).T])
            if rel:
                pedestals /= array([dic['ph'] for dic in self.get_rp_pulse_heights(sel, redo).values()]) * .01
            g = Draw.make_tgraph(fluxes, pedestals, color=self.Draw.get_color(self.NPlans))
            mg.Add(g)
        legend = self.make_full_legend(mg.GetListOfGraphs(), irr)
        format_histo(mg, draw_first=True, y_tit='Pulse Height [au]', tit_size=.05, lab_size=.05, y_off=.91, x_off=1.2, x_range=Bins.FluxRange)
        self.Draw(mg, '{}Pedestals'.format(self.Name), draw_opt='ap', logx=True, leg=legend, w=1.6, lm=.092, bm=.12, gridy=True, show=show)

    def draw_mean_pedestals(self, sigma=False, irr=False, redo=False, show=True):
        y = array([mean_sigma(tc_values[1 if sigma else 0])[0] for tc_values in self.get_pedestals(redo)])
        x = self.get_irradiations(string=False) if irr else arange(y.size)
        g = Draw.make_tgraph(x, y, title='Mean {}'.format('Noise' if sigma else 'Pedestals'), x_tit='Irradation [10^{14} n/cm^{2}]' if irr else 'Run Plan', y_tit='Pulse Height [mV]')
        format_histo(g, y_off=1.2, x_range=ax_range(x, fl=.1, fh=.1) if irr else ax_range(0, y.size - 1, .3, .3), x_off=2.5)
        self.set_bin_labels(g) if not irr else do_nothing()
        self.Draw(g, self.get_savename('PedestalMeans'), show, draw_opt='ap', bm=.2, w=1.5, h=.75, gridy=True)

    def draw_mean_noise(self, irr=False, redo=False, show=True):
        self.draw_mean_pedestals(True, irr, redo, show)
    # endregion PEDESTAL
    # ----------------------------------------

    # ----------------------------------------
    # region SYS ERROR
    def subtract_means(self, err=True, ey=0):
        x, y = self.get_fluxes(avrg=True), self.flux_avrg([add_err(y, ey) for y in self.get_pulse_heights(err=err)])  # rate scans must have the same flux points
        return concatenate(array(x).T), array([iy - mean_sigma(g)[0].n for g in array(y).T for iy in g])

    def mean_max_spread(self):
        x = self.get_mean_ph()
        return (max(x) - min(x)) / mean(x)

    def max_spread(self):
        """ returns: maximum spread of averaged points of similar flux for two rate scans. Both rate scans must have the same flux points."""
        x = self.get_pulse_heights(err=False, avrg=True)
        return max(abs(diff(x, axis=0)[0] / mean(x, axis=0)))

    @save_pickle('RelSErr', field='Name')
    def calc_rel_sys_error(self, _redo=False):
        """ vary the error of the indiviual ph points such that the avrg measurements of two rate scans agree. """
        f = Draw.make_tf1('chi', lambda x: FitRes(Draw.make_tgraph(*self.subtract_means(err=False, ey=x)).Fit('pol0', 'qs')).get_chi2(), 0, 1)
        v = f.GetX(1)  # get the error size such that the chi2 is 1
        return ufloat(v, abs(v - f.GetX(1.1))) / mean_sigma(self.get_mean_ph())[0].n

    def print_err(self, redo=False):
        h = latex.bold('Detector', f'Rel. Sys. Error {latex.unit("percent")}', f'Bias {latex.unit("V")}')
        print(latex.table(h, [[self.DUTName] + latex.si(self.calc_rel_sys_error(_redo=redo) * 100, fmt='.2f') + latex.si(self.Info[0].Bias, fmt='+.0f')], align_header=True))

    def draw_0fit(self, err=True, ey=0, **dkw):
        x, y = self.subtract_means(err, ey)
        g = self.Draw.graph(x, y, y_tit='Shifted Pulse Heigt [mV]')
        g.Fit('pol0', 'qs')
        return self.Draw(g, **prep_kw(dkw, **self.get_x_args(), stats=set_statbox(fit=True), file_name='ShiftedPHFit'))
    # endregion ERROR
    # ----------------------------------------


class SelectionInfo:
    def __init__(self, sel: RunSelector):
        self.RunInfo = sel.RunInfos
        self.Run = sel.Run
        self.TCString = sel.TCString
        self.RunPlan = sel.SelectedRunplan
        self.DUT = sel.SelectedDUT
        self.DUTName = self.DUT.Name
        self.DUTNr = self.DUT.Number
        self.IsPixel = 'pix' in (self.DUT.Type if type(self.DUT.Type) is str else self.DUT.Type[self.TCString])
        self.Verbose = sel.Run.Verbose
        self.Bias = self.DUT.Bias
        self.Irradiation = self.DUT.get_irradiation(self.TCString)
        self.Type = sel.SelectedType.lower()
        self.Runs = sel.get_selected_runs()
        self.PulserType = sel.PulserType
        self.DUT = PixelDUT(self.DUT.Number, self.RunInfo[self.Runs[0]])

    def __str__(self):
        return f'Selection instance: {self.TCString:<8} {self.RunPlan:<4} {self.DUTName}, {bias2str(self.Bias)}'

    def __repr__(self):
        return self.__str__()

    def __call__(self):
        return [self.TCString, self.RunPlan, self.DUTName, self.DUTNr, f'{self.Runs[0]:03d}-{self.Runs[-1]:03d}', f'{self.Bias:+4.0f}V', self.Type, f'{self.Irradiation.n:.1e}']

    @property
    def is_volt_scan(self):
        return 'voltage' in self.Type


class PickleInfo:
    def __init__(self, name=None, *suf):
        self.SubDir = 'selections'
        self.Name = name if name else None
        self.Suffix = make_suffix(None, *suf)

    def __call__(self, *args, **kwargs):
        return self.SubDir is not None

    def path(self, sel: SelectionInfo):
        return sel.Run.make_pickle_path(self.SubDir, self.Name, sel.RunPlan, sel.DUT.Number, self.Suffix, sel.TCString) if self() else ''


if __name__ == '__main__':

    aparser = ArgumentParser()
    aparser.add_argument('sel', nargs='?', default='test-a', help='name of the selection')
    aparser.add_argument('-d', nargs='?', default=None, help='dut')
    aparser.add_argument('-tc', nargs='?', default=None, help='test campaign')
    aparser.add_argument('-v', '--verbose', action='store_false')
    aparser.add_argument('-p', action='store_true', help='print analysis strings')
    aparser.add_argument('-r', action='store_true', help='redo')
    aparser.add_argument('-la', action='store_true', help='active show all selections')
    aparser.add_argument('-ls', nargs='?', default=None, help='active show all selections for given DUT')
    pargs = aparser.parse_args()

    z = DiaScans(pargs.sel, pargs.verbose and not (pargs.ls or pargs.la))
    if pargs.p:
        print(z.get_all_ana_strings(pargs.d, pargs.tc, pargs.r))
    if pargs.ls or pargs.la:
        z.show_selections(pargs.ls)
        critical('finish')
