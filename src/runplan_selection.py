#! /usr/bin/env python
# --------------------------------------------------------
#       DIAMOND RATE SCANS
# created on June 24th 2016 by M. Reichmann
# --------------------------------------------------------

from json import dump
from re import split as splitname
from ROOT import TMultiGraph
from src.binning import Bins
from pad.collection import AnalysisCollection, PadCollection
from src.run_selection import RunSelector
from src.analysis import *
from analyse import collection_selector
from plotting.draw import get_graph_y, ax_range, markers
from inspect import signature


class DiaScans(Analysis):
    def __init__(self, selection_name=None, verbose=False):

        self.Name = selection_name
        Analysis.__init__(self, verbose=verbose, results_dir='selections', sub_dir=selection_name if selection_name is not None else '')

        self.print_start(run=selection_name, tc=False, prnt=verbose)

        # Main
        self.Selections = self.load_selections()
        self.Selection = self.load_selection(selection_name)

        # Config
        self.DUTParser = load_parser(join(self.Dir, 'config', 'DiamondAliases.ini'))

        # Info
        self.RS = RunSelector()  # dummy for information
        self.DUTName = self.load_dut_name()
        self.RunPlans = self.load_runplans()
        self.TestCampaigns = self.load_test_campaigns()
        self.Info = self.load_selection_info()
        self.NPlans = len(self.Info) if self.Info else None

        self.print_finished(prnt=verbose)

    def __str__(self):
        return self.Name

    # ----------------------------------------
    # region INIT
    def load_selections(self):
        with open(join(self.Dir, self.MainConfig.get('MISC', 'runplan selection file'))) as f:
            return load(f, object_pairs_hook=OrderedDict)

    def load_selection(self, name):
        if name is None:
            return
        if name not in self.Selections.keys():
            warning('"{} does not exist in:\n{} '.format(name, sorted(self.Selections.keys())))
            self.Name = None
            return
        return self.Selections[name]

    def load_dut_name(self, name=None):
        name = self.Name if name is None else name
        if name is None:
            return
        if any([word in name.lower() for word in ['all', 'test']]):
            return name.title()
        dut = name.lower()
        if dut not in self.DUTParser.options('ALIASES'):
            dut = splitname('[-_]', name)[0]
        if dut.lower() not in self.DUTParser.options('ALIASES'):
            dut = '-'.join(splitname('[-_]', name)[:2])
        if dut.lower() not in self.DUTParser.options('ALIASES'):
            warning('No diamond name found in "{}". Please choose one from \n{}'.format(name, ', '.join(self.DUTParser.options('ALIASES'))))
        return self.DUTParser.get('ALIASES', dut)

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

    def get_run_types(self):
        return [sel.Type.lower() for sel in self.Info]

    def get_irradiations(self, string=True):
        return array([make_irr_string(sel.Irradiation) if string else float(sel.Irradiation) for sel in self.Info])

    def get_bias_voltages(self):
        return [sel.Bias for sel in self.Info]

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

    def get_pulse_heights(self, avrg=False, redo=False):
        return self.get_values(AnalysisCollection.get_pulse_heights, PickleInfo('PHVals', f'{avrg:d}'), redo=redo, avrg=avrg)

    def get_rate_dependencies(self, redo=False):
        return array(self.get_values(AnalysisCollection.get_rate_dependence, PickleInfo('PHRD'), redo=redo))

    def print_mean_rd(self, latex=False):
        s, r = [mean_sigma(i)[0] * 100 for i in self.get_rate_dependencies().T]
        info(f'Mean rel STD:    {si(s, unit="percent") if latex else f"{s:.1f} %"}')
        info(f'Mean rel Spread: {si(r, unit="percent") if latex else f"{r:.1f} %"}')

    def print_rate_dependencies(self):
        for i, (s1, s2) in zip(self.Info, self.get_rate_dependencies() * 100):
            print(i)
            print(f'  Rel STD:    {s1.n:.1f} %')
            print(f'  Rel Spread: {s2:.1f} %')

    def print_rd_table(self):
        data = 100 * self.get_rate_dependencies().reshape((self.NPlans // 2, 4))
        rows = [[f'\\SI{{{self.Info[2 * i].Irradiation}}}{{}}'] + [f'{s.n:.1f}' for s in [s1, s2]] + [si(r, '.1f') for r in [r1, r2]] for i, (s1, r1, s2, r2) in enumerate(data)]
        print(make_latex_table([], rows))

    def get_rp_pulse_heights(self, sel, corr=True, redo=False):
        return self.get_rp_values(sel, AnalysisCollection.get_pulse_heights, PickleInfo('PHVals', corr), redo=redo, corr=corr)

    def get_pedestals(self, redo=False):
        return self.get_values(PadCollection.get_pedestals, PickleInfo('PedVals'), redo=redo)

    def get_noise_spread(self, redo=False):
        noise = array(self.get_pedestals(redo), object)[:, 1]
        return mean_sigma([v - mean_sigma(lst)[0] for lst in noise for v in lst])[1]

    def get_mean_uniformities(self, use_fcw=True, redo=False, low=False, high=False):
        pickle_info = PickleInfo('Uni', low, high, use_fcw)
        return self.get_values(AnalysisCollection.get_mean_uniformity, pickle_info, redo=redo, high_flux=high, low_flux=low)

    def get_uniformities(self, use_fcw=True, redo=False, low=False, high=False):
        pickle_info = PickleInfo('Uniformity', 'UniSMSTD', '{}{}{}'.format(int(low), int(high), int(use_fcw)))
        return self.get_values(AnalysisCollection.get_uniformities, pickle_info, redo=redo, high_flux=high, low_flux=low, use_fcw=use_fcw)

    def get_currents(self):
        return self.get_values(AnalysisCollection.get_currents, PickleInfo('CurrentsVals'))

    def get_fluxes(self, avrg=False, redo=False):
        return self.get_values(AnalysisCollection.get_fluxes, PickleInfo('FluxVals', avrg), avrg=avrg, redo=redo)

    def get_all_infos(self):
        return [sel for tc in self.RunPlans.keys() for sel in self.get_tc_infos(tc)]

    def get_dia_infos(self, dut_name):
        dut_name = self.RS.Run.translate_dia(dut_name)
        return [sel for tc in self.RunPlans.keys() for sel in self.get_tc_infos(tc) if dut_name == sel.DUTName]

    def get_tc_infos(self, tc):
        rs = RunSelector(tc)
        return [SelectionInfo(rs.select_runs_from_runplan(rp, dut + 1, unselect=True)) for rp in sorted(self.RunPlans[tc]) for dut in range(rs.get_n_duts(run_plan=rp))]

    def get_bias_str(self):
        return ' at {bias} V'.format(bias=self.Info[0].Bias) if len(set(self.get_bias_voltages())) == 1 else ''

    def get_all_ana_strings(self, dut=None, tc=None, redo=False):
        selections = self.get_all_infos() if dut is None and tc is None else self.get_dia_infos(dut) if tc is None else self.get_tc_infos(tc)
        selections = [sel for sel in selections if self.RS.Run.translate_dia(dut) == sel.DUTName] if dut is not None else selections
        redo = ' -rd' if redo else ''
        return '_'.join(['analyse {} {} -c -tc {} -d{}'.format(sel.RunPlan, sel.DUTNr, sel.TCString, redo) for sel in selections])

    def get_savename(self, name):
        return '{}{}'.format(name, self.Name.title().replace('-', '').replace('_', ''))

    def get_signal_maps(self, fid=False, res=.2, square=True, scale=True, redo=False):
        pickle_info = PickleInfo('SM', make_suffix(self, [fid, res, square, scale]))
        return self.get_values(AnalysisCollection.draw_signal_map, pickle_info, fid=fid, res=res, square=square, scale=scale, show=False, redo=redo)
    # endregion GET
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

    def show_selection(self):
        """ Gives detailed information about the chosen selection """
        print_table(rows=[sel() for sel in self.Info], header=['TC', 'RunPlan', 'DUT', 'Nr', 'Runs', 'Bias', 'Type', 'Irrad']) if self.Info else warning('Selection is empty!')

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
        self.DUTName = self.load_dut_name()
        self.Selection = self.load_selection(name)
        self.TestCampaigns = self.load_test_campaigns()
        self.Info = self.load_selection_info()
        self.Draw.set_results_dir(join('Results', 'selections', self.Name)) if name else do_nothing()

    def set_selection_name(self, name):
        self.Name = name
        self.Selection = self.load_selection(name)

    def clear_selection(self):
        self.Selection = {}

    def select_runplan(self, runplan, dut=1, testcampaign=None):
        rp = make_runplan_string(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        if rp in self.RunPlans[tc]:
            if tc not in self.Selection:
                self.Selection[tc] = {}
            self.Selection[tc][rp] = dut
        else:
            warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def unselect_runplan(self, runplan, testcampaign=None):
        rp = make_runplan_string(runplan)
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
        with open(join(self.Dir, self.MainConfig.get('MISC', 'runplan selection file')), 'w') as f:
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
    def draw_dia_rate_scans(self, redo=False, irr=True, corr=True):
        mg = TMultiGraph('mg_ph', '{dia} Rate Scans{b};Flux [kHz/cm^{{2}}]; Pulse Height [mV]'.format(dia=self.DUTName, b=self.get_bias_str()))
        mgs = self.get_values(AnalysisCollection.draw_pulse_heights, PickleInfo('PHMG', '10000_{}'.format(corr)), redo=redo, show=False)
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

    def draw_pedestals(self, rel=False, redo=False, show=True, irr=True):
        mg = TMultiGraph('mg_ph', '{dia} Pedestals{b};Flux [kHz/cm^{{2}}]; Pulse Height [mV]'.format(dia=self.DUTName, b=self.get_bias_str()))
        for i, (values, sel, fluxes) in enumerate(zip(self.get_pedestals(redo), self.Info, self.get_fluxes())):
            pedestals = array([make_ufloat(*tup) for tup in array(values).T])
            if rel:
                pedestals /= array([dic['ph'] for dic in self.get_rp_pulse_heights(sel, redo).values()]) * .01
            g = Draw.make_tgrapherrors(fluxes, pedestals, color=self.Draw.get_color(self.NPlans))
            mg.Add(g)
        legend = self.make_full_legend(mg.GetListOfGraphs(), irr)
        format_histo(mg, draw_first=True, y_tit='Pulse Height [au]', tit_size=.05, lab_size=.05, y_off=.91, x_off=1.2, x_range=Bins.FluxRange)
        self.Draw(mg, '{}Pedestals'.format(self.Name), draw_opt='ap', logx=True, leg=legend, w=1.6, lm=.092, bm=.12, gridy=True, show=show)

    def make_full_legend(self, graphs, irr=True):
        same_bias = len(set(self.get_bias_voltages())) == 1
        cols = 1 + (not same_bias) + irr
        legend = Draw.make_legend(y2=.4, w=.12 * cols, nentries=4, cols=cols)
        for i, (g, sel) in enumerate(zip(graphs, self.Info)):
            legend.AddEntry(g, '{} - {}'.format(sel.RunPlan, make_tc_str(sel.TCString, long_=False)), 'lp')
            legend.AddEntry(0, make_irr_string(sel.Irradiation), '') if irr else do_nothing()
            legend.AddEntry(0, get_bias_root_string(sel.Bias), '') if not same_bias else do_nothing()
        return legend

    def draw_title_pad(self, h, x0, lm, c_height):
        if Draw.Title:
            biases = list(set(self.get_bias_voltages()))
            bias_str = ' at {b}'.format(b=make_bias_str(biases[0])) if len(biases) == 1 else ''
            Draw.tpad('p0', pos=[x0, 1 - h / c_height, 1, 1], margins=[0, 0, 0, 0], transparent=True)
            Draw.tpavetext('{dia} Rate Scans{b}'.format(dia=self.DUTName, b=bias_str), lm, 1, 0, 1, font=62, align=13, size=.5, margin=0)
            get_last_canvas().cd()

    def draw_scaled_rate_scans(self, irr=False, y_range=.07, pad_height=.18, avrg=False, **dkw):
        title_height = pad_height / 2 if Draw.Title else .03        # half of a pad for title
        c_height = (self.NPlans + .5) * pad_height + title_height   # half of a pad for the x-axis
        c_width = 1.3 * pad_height / .2                             # keep aspect ratio for standard pad_height
        c = Draw.canvas(w=c_width, h=c_height, transp=True, logx=True, gridy=True)
        lm, rm, x0, size = .07, .02, .08, .22
        self.draw_title_pad(title_height, x0, lm, c_height)
        Draw.tpad('p1', pos=[0, 0, x0, 1], margins=[0, 0, 0, 0], transparent=True)           # info pad
        Draw.tpavetext('Scaled Pulse Height', 0, 1, 0, 1, align=22, size=.5, angle=90, margin=0)   # y-axis title

        for i, g in enumerate(self.get_values(AnalysisCollection.draw_scaled_pulse_heights, PickleInfo('PHScaledGraph', avrg), avrg=avrg, show=False)):
            c.cd()
            y0, y1 = [(c_height - title_height - pad_height * (i + j)) / c_height for j in [1, 0]]
            p = Draw.tpad(pos=[x0, y0, 1, y1], margins=[lm, rm, 0, 0], logx=True, gridy=True, gridx=True, fix=True)
            self.Draw(g, title=' ', **prep_kw(dkw, draw_opt='ap', canvas=p, logx=True, gridy=True, color=self.Draw.get_color(self.NPlans), marker=markers(i), markersize=1.5,
                      y_range=[1 - y_range, 1 + y_range], lab_size=size, ndivy=504, x_ticks=.15, x_range=Bins.FluxRange))
            self.draw_legend(i, g, irr, p)
            c.cd()

        Draw.tpad('p2', pos=[x0, 0, 1, pad_height / 2 / c_height], margins=[lm, rm, 0, 0], transparent=True)  # x-axis pad
        Draw.x_axis(1, lm, 1 - rm, 'Flux [kHz/cm^{2}]', Bins.FluxRange, opt='', log=True, tick_size=0, lab_size=size * 2, tit_size=size * 2, off=1.1)
        self.Draw.save_plots('ScaledDiaScans{dia}'.format(dia=make_dia_str(self.DUTName)))

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
            g = Draw.make_tgrapherrors(x, y)
            if align:
                fit = g.Fit('pol1', 'qs0')
                g = Draw.make_tgrapherrors(fluxes[i], array(currents[i]) - fit.Parameter(0) + .1)
            format_histo(g, color=self.Draw.get_color(self.NPlans))
            legend.AddEntry(g, '{tc} - {hv}'.format(tc=self.Info[i].TCString, hv=self.get_rp_values(self.Info[i], AnalysisCollection.get_hv_name, load_tree=False)), 'pl')
            mg.Add(g)
        format_histo(mg, draw_first=True, y_tit='Current [nA]', x_tit='Flux [kHz/cm^{2}]', y_range=[.1, max(concatenate(currents)).n * 2], x_range=Bins.FluxRange)
        self.Draw(mg, 'CurrentFlux{}'.format(self.Name), draw_opt='ap', logx=True, logy=True, leg=legend, bm=.17, show=show)

    def get_titles(self, irr=False):
        if len(set(self.get_dut_names())) > 1:
            return [f'{dut}{" - " if irr else ""}{irr}' for dut, irr in zip(self.get_dut_names(), self.get_irradiations() if irr else [''] * self.NPlans)]
        tits = self.get_irradiations() if irr else [make_tc_str(tc) for tc in self.TestCampaigns]
        if any(['rand' in word for word in self.get_run_types()]):
            for i, sel in enumerate(self.Info):
                tits[i] += ' (random)' if 'rand' in sel.Type.lower() else '         '
        return [[t] + ([] if len(set(self.get_bias_voltages())) == 1 else make_bias_str(bias)) for t, bias in zip(tits, self.get_bias_voltages())]

    def draw_legend(self, i, gr, irr, c):
        tits = self.get_titles(irr)
        w = max([sum(len(t) for t in tit) for tit in tits]) * (.011 if irr else .022)
        Draw.legend([gr] + ([] if len(tits[i]) == 1 else ['']), tits[i], ['pe', ''], w=w, ts=.22, d=0, nentries=1.2, scale=5, c=c, cols=len(tits[i]))

    def set_bin_labels(self, h):
        for i, sel in enumerate(self.Info):
            h.GetXaxis().SetBinLabel(h.GetXaxis().FindBin(i), '{} - {}'.format(make_tc_str(sel.TCString, 0), sel.RunPlan))

    def draw_mean_pedestals(self, sigma=False, irr=False, redo=False, show=True):
        y = array([mean_sigma(tc_values[1 if sigma else 0])[0] for tc_values in self.get_pedestals(redo)])
        x = self.get_irradiations(string=False) / 1e14 if irr else arange(y.size)
        g = Draw.make_tgrapherrors(x, y, title='Mean {}'.format('Noise' if sigma else 'Pedestals'), x_tit='Irradation [10^{14} n/cm^{2}]' if irr else 'Run Plan', y_tit='Pulse Height [mV]')
        format_histo(g, y_off=1.2, x_range=ax_range(x, fl=.1, fh=.1) if irr else ax_range(0, y.size - 1, .3, .3), x_off=2.5)
        self.set_bin_labels(g) if not irr else do_nothing()
        self.Draw(g, self.get_savename('PedestalMeans'), show, draw_opt='ap', bm=.2, w=1.5, h=.75, gridy=True)

    def draw_mean_noise(self, irr=False, redo=False, show=True):
        self.draw_mean_pedestals(True, irr, redo, show)

    def draw_means(self, y_range=None, show=True):
        y = array([make_ufloat(mean_sigma(ph_list, err=False)) for ph_list in self.get_pulse_heights()])
        g = Draw.make_tgrapherrors(arange(y.size), y, title='Pulse Height Evolution', x_tit='Run Plan', y_tit='Mean Pulse Height [mV]')
        format_histo(g, y_off=1.2, x_range=ax_range(0, y.size - 1, .3, .3), x_off=2.5, y_range=y_range)
        self.set_bin_labels(g)
        self.Draw(g, self.get_savename('Means'), show, draw_opt='ap', bm=.2, w=1.5, h=.75, gridy=True)

    def draw_sigmas(self, y_range=None, show=True):
        y = array([mean_sigma(ph_list)[1] for ph_list in self.get_pulse_heights()])
        g = Draw.make_tgrapherrors(arange(y.size), y=y, title='Pulse Height STD Evolution', x_tit='Run Plan', y_tit='Pulse Height Standard Deviation [%]')
        format_histo(g, y_off=1.2, x_range=ax_range(0, y.size - 1, .3, .3), x_off=2.5, y_range=y_range)
        self.set_bin_labels(g)
        self.Draw(g, self.get_savename('Sigmas'), show, draw_opt='ap', bm=.2, w=1.5, h=.75, gridy=True)

    def draw_uniformity(self, arg=2, use_fcw=True, redo=False, low=False, high=False):
        ph_var = 'FWC' if use_fcw else 'Mean'
        var, unit, tit = [ph_var, 'FWHM', 'FWHM/{}'.format(ph_var)][arg], ' [mV]' if arg < 2 else '', [ph_var, 'Full Width Half Maximum', 'Uniformity'][arg]
        y_values = [v[arg][0] for v in self.get_mean_uniformities(redo=redo, low=low, high=high)]
        x_values = array([ufloat(v, v * .2) for v in self.get_irradiations(string=False)]) / 1e15
        return self.Draw.graph(x_values, y_values, title=tit, x_tit='Irradiation [10^{15}n/cm^{2}]', y_tit='{}{}'.format(var, unit), draw_opt='ap', x_off=1.2)

    def draw_peak_flux(self, show=True):
        leg = Draw.make_legend(x2=.45, w=.3)
        mg = TMultiGraph('mgph', 'Peak Flux vs. FAST-OR Flux')
        values_list = self.get_values(PadCollection.get_peak_fluxes, PickleInfo('PeakFlux'))
        flux_list = self.get_fluxes()
        for sel, values, fluxes in zip(self.Info, values_list, flux_list):
            g = Draw.make_tgrapherrors('g{}'.format(sel.RunPlan), '', x=fluxes, y=values)
            format_histo(g, color=self.Draw.get_color())
            leg.AddEntry(g, '{} @ {:+1.0f}V'.format(sel.DUTName, sel.Bias), 'pl')
            mg.Add(g, 'p')
        x, y = concatenate(flux_list), concatenate(values_list)
        x_range, y_range = [.5 * min(x).n, 1.2 * max(x).n], [.5 * min(y).n, 1.2 * max(y).n]
        format_histo(mg, draw_first=True, y_tit='Peak Flux [kHz/cm^{2}] ', x_tit='FAST-OR Flux [kHz/cm^{2}]', x_range=x_range, y_off=1.8, y_range=y_range)
        self.Draw(mg, 'PeakFluxes{}'.format(self.Name), draw_opt='a', leg=leg, show=show, lm=.13)
    # endregion DRAWING
    # ----------------------------------------

    # ----------------------------------------
    # region PULSER
    def get_pulser_pulse_heights(self, avrg=False, redo=False):
        return self.get_values(PadCollection.get_pulser_pulse_heights, PickleInfo('PulserPH', f'{avrg:d}'), redo=redo, avrg=avrg)

    def get_pulser_stability(self):
        """ returns: relative standard deviation of the pulser """
        d = [s / m * 100 for m, s in [mean_sigma(v) for v in self.get_pulser_pulse_heights()]]
        for i, v in enumerate(d):
            print(f'Pulser stability for {self.Info[i]}: {v:.3f}')
        return d

    def draw_pulser_pulse_heights(self, avrg=False, scaled=True, ym=.01):
        data = zip(self.get_fluxes(avrg), self.get_pulser_pulse_heights(avrg))
        g = [self.Draw.graph(x, y / (mean_sigma(y)[0].n if scaled else 1), y_tit='Pulser Pulse Height [mV]', show=False) for x, y in data]
        mg = self.Draw.multigraph(g, 'Pulser Pulse Heights', [f'{i.PulserType}al @ {make_bias_str(i.Bias)}' for i in self.Info], **AnalysisCollection.get_x_args(draw=True), wleg=.3)
        format_histo(mg, **AnalysisCollection.get_x_args(), y_range=1 + array([-ym, ym]))
    # endregion PULSER
    # ----------------------------------------

    def draw_bucket_ratios(self, avrg=False, redo=False, **kwargs):
        g = self.get_values(PadCollection.draw_bucket_ratio, PickleInfo('BucketRatio', f'{avrg:d}'), redo=redo, avrg=avrg, show=False, stats=False)
        mg = self.Draw.multigraph(g, 'Bucket Ratios', self.get_titles(irr=True), **AnalysisCollection.get_x_args(draw=True), wleg=.3)
        format_histo(mg, **prep_kw(kwargs, **AnalysisCollection.get_x_args()))


class SelectionInfo:
    def __init__(self, sel: RunSelector):
        self.Run = sel.Run
        self.TCString = sel.TCString
        self.RunPlan = sel.SelectedRunplan
        self.DUT = sel.SelectedDUT
        self.DUTName = self.DUT.Name
        self.DUTNr = self.DUT.Number
        self.Verbose = sel.Run.Verbose
        self.Bias = self.DUT.Bias
        self.Irradiation = self.DUT.get_irradiation(self.TCString)
        self.Type = sel.SelectedType.lower()
        self.Runs = sel.get_selected_runs()
        self.PulserType = sel.PulserType

    def __str__(self):
        return f'Selection instance: {self.TCString:<8} {self.RunPlan:<4} {self.DUTName}, {make_bias_str(self.Bias)}'

    def __repr__(self):
        return self.__str__()

    def __call__(self):
        return [self.TCString, self.RunPlan, self.DUTName, self.DUTNr, '{:03d}-{:03d}'.format(self.Runs[0], self.Runs[-1]), '{:+4.0f}V'.format(self.Bias), self.Type, self.Irradiation]


class PickleInfo:
    def __init__(self, name=None, *suf):
        self.SubDir = 'selections'
        self.Name = name if name else None
        self.Suffix = '_'.join(str(int(i) if type(i) is bool else i) for i in suf) if suf else None

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
    aparser.add_argument('-s', action='store_true', help='activate show single selection')
    aparser.add_argument('-la', action='store_true', help='active show all selections')
    aparser.add_argument('-ls', default=None, help='active show all selections for given DUT')
    pargs = aparser.parse_args()

    z = DiaScans(pargs.sel, pargs.verbose and not (pargs.ls or pargs.la))
    if pargs.p:
        print(z.get_all_ana_strings(pargs.d, pargs.tc, pargs.r))
    if pargs.s:
        z.show_selection()
    if pargs.ls or pargs.la:
        z.show_selections(pargs.ls)
        critical('finish')
