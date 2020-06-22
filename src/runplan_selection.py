#! /usr/bin/env python
# --------------------------------------------------------
#       DIAMOND RATE SCANS
# created on June 24th 2016 by M. Reichmann
# --------------------------------------------------------

from __future__ import print_function
from analysis import *
from collections import Counter
from json import load, dump

from ROOT import TMultiGraph, TH2F, TH1F, TGraph2DErrors
from re import split as splitname

from pad_collection import AnalysisCollection, PadCollection
from run import Run
from run_selection import RunSelection
from selector import collection_selector
from functools import partial
from binning import Bins


class DiaScans(Analysis):
    def __init__(self, selection_name=None, verbose=False):
        Analysis.__init__(self, verbose=verbose)

        self.print_start(run=selection_name, tc=False, prnt=verbose)

        # Main
        self.Name = selection_name
        self.Selections = self.load_selections()
        self.Selection = self.load_selection(selection_name)

        # Config
        self.DUTParser = load_parser(join(self.Dir, 'config', 'DiamondAliases.ini'))

        # Info
        self.RS = RunSelection()  # dummy for information
        self.DUTName = self.load_dut_name()
        self.RunPlans = self.load_runplans()
        self.TestCampaigns = self.load_test_campaigns()
        self.RunSelections = self.load_run_selections()
        self.RunInfos = self.load_run_infos()
        self.Info = self.load_selection_info()
        self.NPlans = len(self.Info) if self.Info else None

        self.Colors = get_color_gradient(self.NPlans) if self.Info else None

        self.set_results_dir(join('Results', 'selections', '' if self.Name is None else self.Name))

        self.print_finished(prnt=verbose)

    # ----------------------------------------
    # region INIT
    def load_selections(self):
        with open(join(self.Dir, self.MainConfig.get('MISC', 'runplan selection file'))) as f:
            return load(f, object_pairs_hook=OrderedDict)

    def load_selection(self, name):
        if name is None:
            return
        if name not in self.Selections.keys():
            warning('"{} does not exist in:\n{} '.format(name, sorted(self.Selections.iterkeys())))
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
        return sorted(self.RunPlans.iterkeys() if self.Name is None else list(set(self.Selection.iterkeys())))

    def load_run_selections(self):
        return OrderedDict((tc, RunSelection(tc)) for tc in self.TestCampaigns)

    def load_run_infos(self):
        return OrderedDict((tc, rs.RunInfos) for tc, rs in self.RunSelections.iteritems())

    def load_selection_info(self):
        selections = []
        if self.Selection is None:
            return
        for tc, rps in self.Selection.iteritems():
            for rp, dut_nrs in rps.iteritems():
                for dut_nr in array([dut_nrs]).flatten():
                    selections.append(SelectionInfo(self.RunSelections[tc].select_runs_from_runplan(rp, dut_nr, unselect=True)))
        return selections
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_rp_diamonds(self, tc, rp):
        sel = self.RunSelections[tc].select_runs_from_runplan(rp, unselect=True)
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
        return [make_irr_string(sel.Irradiation) if string else float(sel.Irradiation) for sel in self.Info]

    def get_bias_voltages(self):
        return [sel.Bias for sel in self.Info]

    def get_rp_values(self, sel, f, pickle_info=None, redo=False, load_tree=True, *args, **kwargs):
        pickle_path = self.make_pickle_path(pickle_info.SubDir, pickle_info.Name, sel.RunPlan, sel.DUTNr, pickle_info.Suffix, sel.TCString) if pickle_info else ''
        if file_exists(pickle_path) and not redo:
            with open(pickle_path) as f:
                return pickle.load(f)
        self.info('Did not find {}'.format(pickle_path), prnt=pickle_path)
        ana = collection_selector(sel.RunPlan, sel.DUTNr, sel.TCString, load_tree)
        try:
            pf = partial(f, ana, redo=redo, *args, **kwargs)
            return do_pickle(pickle_path, pf, redo=redo) if pickle_info else pf()
        except TypeError:
            pf = partial(f, ana, *args, **kwargs)
            return do_pickle(pickle_path, pf, redo=redo) if pickle_info else pf()

    def get_values(self, f, pickle_info=None, redo=False, load_tree=True, *args, **kwargs):
        return [self.get_rp_values(sel, f, pickle_info, redo, load_tree, *args, **kwargs) for sel in self.Info]

    def get_pulse_heights(self, avrg=False, redo=False):
        return self.get_values(AnalysisCollection.get_pulse_heights, PickleInfo('Ph_fit', 'PhVals', '{}'.format(int(avrg))), redo=redo, avrg=avrg)

    def get_rate_dependcies(self, redo=False):
        return self.get_values(AnalysisCollection.get_rate_dependence, PickleInfo('Ph_fit', 'RD'), redo=redo)

    def print_rate_dependencies(self):
        for i, (s1, s2) in zip(self.Info, self.get_rate_dependcies()):
            print(i)
            print('  Rel STD:    {:2.1f}'.format(s1.n * 100))
            print('  Rel Spread: {:2.1f} \\pm {:0.1f}'.format(s2.n * 100, s2.s * 100))

    def get_rp_pulse_heights(self, sel, corr=True, redo=False):
        return self.get_rp_values(sel, AnalysisCollection.get_pulse_heights, PickleInfo('Ph_fit', 'PhVals', '10000_{}'.format(corr)), redo=redo, corr=corr)

    def get_pedestals(self, redo=False):
        return self.get_values(PadCollection.get_pedestals, PickleInfo('Pedestal', 'Values'), redo=redo)

    def get_rel_errors(self, flux=105, redo=False):
        return self.get_values(AnalysisCollection.get_repr_error, PickleInfo('Errors', 'Repr', flux), redo=redo, flux=flux)

    def get_mean_uniformities(self, redo=False, low=False, high=False):
        return self.get_values(AnalysisCollection.get_mean_uniformity, PickleInfo('Uniformity', '', '{}{}'.format(int(low), int(high))), redo=redo, high_flux=high, low_flux=low)

    def get_uniformities(self, redo=False, low=False, high=False):
        return self.get_values(AnalysisCollection.get_uniformities, PickleInfo('Uniformity', 'SMSTD', '{}{}'.format(int(low), int(high))), redo=redo, high_flux=high, low_flux=low)

    def get_currents(self):
        return self.get_values(AnalysisCollection.get_currents, PickleInfo('Currents', 'Vals'))

    def get_fluxes(self, avrg=False):
        return self.get_values(AnalysisCollection.get_fluxes, PickleInfo('Flux', 'Vals', suf='{}'.format(int(avrg)) if avrg else ''), avrg=avrg)

    def get_all_infos(self):
        return [sel for tc in self.RunPlans.iterkeys() for sel in self.get_tc_infos(tc)]

    def get_dia_infos(self, dut_name):
        dut_name = self.RS.Run.translate_dia(dut_name)
        return [sel for tc in self.RunPlans.iterkeys() for sel in self.get_tc_infos(tc) if dut_name == sel.DUTName]

    def get_tc_infos(self, tc):
        rs = self.RunSelections[tc] if tc in self.RunSelections else RunSelection(tc)
        return [SelectionInfo(rs.select_runs_from_runplan(rp, dut + 1, unselect=True)) for rp in sorted(self.RunPlans[tc]) for dut in xrange(rs.get_n_duts(run_plan=rp))]

    def get_bias_str(self):
        return ' at {bias} V'.format(bias=self.Info[0].Bias) if len(set(self.get_bias_voltages())) == 1 else ''

    def get_all_ana_strings(self, dut=None, tc=None, redo=False):
        selections = self.get_all_infos() if dut is None and tc is None else self.get_dia_infos(dut) if tc is None else self.get_tc_infos(tc)
        selections = [sel for sel in selections if self.RS.Run.translate_dia(dut) == sel.DUTName] if dut is not None else selections
        redo = ' -rd' if redo else ''
        return '_'.join(['analyse {} {} -c -tc {} -d{}'.format(sel.RunPlan, sel.DUTNr, sel.TCString, redo) for sel in selections])
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW
    def show_selections(self):
        header = ['Name', 'Diamond', 'Campaigns']
        rows = []
        old_sel = deepcopy(self.Name)
        for name in self.Selections.iterkeys():
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
        for tc, runplans in sorted(self.RunPlans.iteritems()):
            if runplans:
                print_small_banner(tc, color='yellow')
                header = ['Runplan', 'Runs', 'Diamonds']
                rows = []
                for rp, dic in sorted(runplans.iteritems()):
                    runs = dic['runs']
                    rows.append([rp, '{:03d}-{:03d}'.format(runs[0], runs[-1]), ', '.join(self.get_rp_diamonds(tc, rp))])
                print_table(rows, header)
        self.set_selection_name(old_sel)

    def show_pulse_heights(self):
        rows = [[sel.TCString, sel.RunPlan] + ['{:2.2f}'.format(i) for i in mean_sigma([dic['ph'] for dic in phs.itervalues()])] for phs, sel in zip(self.get_pulse_heights(), self.Info)]
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
        self.RunSelections = self.load_run_selections()
        self.RunInfos = [rs.RunInfos for rs in self.RunSelections.itervalues()]
        self.Info = self.load_selection_info()
        self.set_results_dir(join('Results', 'selections', self.Name)) if name else do_nothing()

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
            log_warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def unselect_runplan(self, runplan, testcampaign=None):
        rp = make_runplan_string(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        try:
            self.Selection[tc].pop(rp)
        except KeyError:
            log_warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def add_selection(self):
        name = raw_input('Enter the name of the selection: ')
        self.clear_selection()
        self.Selections[name] = {} if name not in self.Selections else self.Selections[name]
        run_plan = raw_input('Enter test campaign, run plan number, channel: ')
        while run_plan:
            tc, rp, ch = [string.strip(' ') for string in run_plan.split(',')]
            self.select_runplan(rp, int(ch), tc)
            run_plan = raw_input('Enter test campaign, run plan number, channel (leave blank to finish): ')
        self.save_selection(name)

    def save_selection(self, name=None):
        name = raw_input('Enter the name of the selection: ') if name is None else name
        if not self.Selection:
            warning('Selection is empty!')
            return
        with open(join(self.Dir, self.MainConfig.get('MISC', 'runplan selection file')), 'w') as f:
            if name in self.Selections:
                query = raw_input('{} does already exist. Do you want to overwrite it? (y/n) '.format(name))
                if query.lower().strip() in ['no', 'n']:
                    return
            self.Selections[name] = self.Selection
            dump(self.Selection, f, indent=2, sort_keys=True)
            self.info('Saved {} to selections'.format(name))
    # endregion SELECTION
    # ----------------------------------------

    # ----------------------------------------
    # region DRAWING
    def draw_collimator_settings(self, show=True):
        h = TH2F('h_cs', 'Collimator Settings', 125, 50, 300, 75, 0, 150)
        for tc in self.TestCampaigns:
            for run, data in self.RunInfos[tc].iteritems():
                try:
                    h.Fill(data['fs11'], data['fs13'])
                except KeyError:
                    pass
        format_histo(h, x_tit='fs11', y_tit='fsh13', y_off=1.3, stats=0, z_off=1.1, z_tit='Number of Entries', z_range=[0, 80])
        self.save_histo(h, 'CollimatorSettings', show, draw_opt='colz', lm=.12, rm=.16)

    def draw_flux_vs_collimators(self, show=True):
        gr = TGraph2DErrors()
        gr.SetNameTitle('gr_fc', 'Flux Vs. Collimators')
        col_settings = Counter([(data['fs11'], data['fs13']) for tc in self.TestCampaigns for data in self.RunInfos[tc].itervalues() if 'fs11' in data and data['fs11'] > 0])
        i = 0
        for col, nr in col_settings.iteritems():
            if nr > 10:
                flux_fit = self.draw_flux_distribution(col[0], col[1], show=False)
                if flux_fit is not None:
                    gr.SetPoint(i, col[0], col[1], flux_fit.Parameter(1))
                    gr.SetPointError(i, 0, 0, flux_fit.Parameter(2))
                    i += 1
        self.draw_histo(gr, 'FluxVsCollimators', show, draw_opt='surf1', lm=.15, phi=17, theta=35)
        format_histo(gr, x_tit='fs11', x_off=1.3, y_tit='fsh13', y_off=1.9, stats=0, z_off=1.9, z_tit='Flux kHz/cm^{2}', markersize=2, y_range=[0, 130])
        self.save_plots('FluxVsCollimators', show=show, prnt=False)
        h = gr.Clone()
        h.Draw('samep')
        self.Objects.append(h)
        self.save_plots('FluxVsCollimators', show=show)

    def draw_flux_variations(self, show=True, rel_sigma=False):
        gr = self.make_tgrapherrors('gr_fd', 'Flux Deviations')
        col_settings = Counter([(data['fs11'], data['fs13']) for tc in self.TestCampaigns for data in self.RunInfos[tc].itervalues() if 'fs11' in data and data['fs11'] > 0])
        i = 0
        for col, nr in sorted(col_settings.iteritems()):
            if nr > 30:
                flux_fit = self.draw_flux_distribution(col[0], col[1], show=False)
                if flux_fit is not None:
                    gr.SetPoint(i, flux_fit.Parameter(1), flux_fit.Parameter(2) if not rel_sigma else flux_fit.Parameter(2) / flux_fit.Parameter(1))
                    yerr = flux_fit.ParError(2) + .5 * flux_fit.Parameter(2)
                    if rel_sigma:
                        yerr = flux_fit.Parameter(2) / flux_fit.Parameter(1) * sqrt(sum([((flux_fit.ParError(j) + .5 * flux_fit.Parameter(2) if rel_sigma else 0) / flux_fit.Parameter(j)) ** 2
                                                                                         for j in xrange(1, 3)]))
                    gr.SetPointError(i, flux_fit.ParError(1), yerr)
                    l1 = self.draw_tlatex(gr.GetX()[i] * 1.05, gr.GetY()[i], '{0}/{1}'.format(make_col_str(col[0]), make_col_str(col[1])), color=1, align=10, size=.03)
                    gr.GetListOfFunctions().Add(l1)
                    i += 1
        format_histo(gr, x_tit='Mean Flux [au]', y_tit='{0}Sigma [au]'.format('Relative ' if rel_sigma else ''), y_off=1.3)
        self.draw_histo(gr, 'FluxVariations', show, draw_opt='alp', logx=True, logy=not rel_sigma, lm=.12)
        gr.GetXaxis().SetLimits(gr.GetX()[0] / 2, gr.GetX()[gr.GetN() - 1] * 4)
        self.save_plots('FluxVariations{0}'.format('Rel' if rel_sigma else ''), show=show)

    def draw_flux_distribution(self, fs11, fsh13, tc=None, do_fit=True, show=True, run_thr=None):
        values = []
        for tc in self.TestCampaigns if tc is None else [tc]:
            for run, data in sorted(self.RunInfos[tc].iteritems()):
                info_run = Run(number=run, test_campaign=tc, tree=False)
                if run_thr is not None:
                    if run_thr > 0 and int(run) < run_thr:
                        continue
                    elif run_thr < 0 and int(run) > abs(run_thr):
                        continue
                try:
                    if data['fs11'] == fs11 and data['fs13'] == fsh13:
                        flux = info_run.Flux
                        # print tc, run, flux
                        values.append(flux) if flux > 1 else do_nothing()
                except KeyError:
                    pass
        if not values:
            return
        spread = max(values) - min(values)
        set_root_output(False)
        h = TH1F('h_fd', 'Flux Distribution for {0}/{1}'.format(fs11, fsh13), int(sqrt(len(values))) + 5, min(values) - .2 * spread, max(values) + .2 * spread)
        for val in values:
            h.Fill(val)
        self.format_statbox(only_fit=True, w=.25) if do_fit else do_nothing()
        format_histo(h, x_tit='Flux in kHz/cm^{2}', y_tit='Number of Entries', y_off=1.3, stats=0 if not do_fit else 1)
        self.draw_histo(h, '', show)
        fit = None
        if do_fit:
            h.SetName('Fit Results')
            set_root_output(show)
            fit = h.Fit('gaus', 'qs')
        self.save_plots('FluxDistribution{0}_{1}'.format(int(fs11), int(fsh13)), show=show)
        return fit

    def draw_dia_rate_scans(self, redo=False, irr=True, corr=True):
        mg = TMultiGraph('mg_ph', '{dia} Rate Scans{b};Flux [kHz/cm^{{2}}]; Pulse Height [mV]'.format(dia=self.DUTName, b=self.get_bias_str()))
        mgs = self.get_values(AnalysisCollection.draw_pulse_heights, PickleInfo('Ph_fit', 'MG', '10000_{}'.format(corr)), redo=redo, show=False, prnt=False)
        for i, (mgi, sel) in enumerate(zip(mgs, self.Info)):
            for g in mgi.GetListOfGraphs():
                format_histo(g, color=self.Colors[i], markersize=1.5, lw=2)
                if g.GetName() == 'gFirst':
                    format_histo(g, color=1, marker=26, markersize=2)
                elif g.GetName() == 'gLast':
                    format_histo(g, color=1, marker=23, markersize=2)
            mg.Add(mgi)
        legend = self.make_full_legend([mgi.GetListOfGraphs()[0] for mgi in mgs], irr)
        y = concatenate([get_graph_y(g) for g in mg.GetListOfGraphs()])
        format_histo(mg, draw_first=True, y_tit='Pulse Height [au]', y_range=[0, y.max().n * 1.1], tit_size=.05, lab_size=.05, y_off=.91, x_off=1.2, x_range=Bins().FluxRange)
        self.save_histo(mg, 'DiaScans{dia}'.format(dia=make_dia_str(self.DUTName)), draw_opt='a', logx=True, leg=legend, x=1.6, lm=.092, bm=.12, gridy=True)

    def draw_pedestals(self, rel=False, redo=False, show=True, irr=True):
        mg = TMultiGraph('mg_ph', '{dia} Pedestals{b};Flux [kHz/cm^{{2}}]; Pulse Height [mV]'.format(dia=self.DUTName, b=self.get_bias_str()))
        for i, (values, sel, fluxes) in enumerate(zip(self.get_pedestals(redo), self.Info, self.get_fluxes())):
            pedestals = array([make_ufloat(*tup) for tup in array(values).T])
            if rel:
                pedestals /= array([dic['ph'] for dic in self.get_rp_pulse_heights(sel, redo).itervalues()]) * .01
            g = self.make_tgrapherrors('gp{}'.format(i), '', x=fluxes.values(), y=pedestals, color=self.Colors[i])
            mg.Add(g, 'pl')
        legend = self.make_full_legend(mg.GetListOfGraphs(), irr)
        format_histo(mg, draw_first=True, y_tit='Pulse Height [au]', tit_size=.05, lab_size=.05, y_off=.91, x_off=1.2, x_range=Bins().FluxRange)
        self.save_histo(mg, '{}Pedestals'.format(self.Name), draw_opt='a', logx=True, leg=legend, x=1.6, lm=.092, bm=.12, gridy=True, show=show)

    def make_full_legend(self, graphs, irr=True):
        same_bias = len(set(self.get_bias_voltages())) == 1
        cols = 1 + (not same_bias) + irr
        legend = self.make_legend(.75, .4, w=.2 * cols, nentries=4, clean=True, cols=cols)
        for i, (g, sel) in enumerate(zip(graphs, self.Info)):
            legend.AddEntry(g, '{} - {}'.format(sel.RunPlan, make_tc_str(sel.TCString)), 'lp')
            legend.AddEntry(0, make_irr_string(sel.Irradiation), '') if irr else do_nothing()
            legend.AddEntry(0, get_bias_root_string(sel.Bias), '') if not same_bias else do_nothing()
        return legend

    def draw_title_pad(self, h, x0, lm, c_height):
        if self.Title:
            biases = list(set(self.get_bias_voltages()))
            bias_str = ' at {b}'.format(b=make_bias_str(biases[0])) if len(biases) == 1 else ''
            self.draw_tpad('p0', 'p0', pos=[x0, 1 - h / c_height, 1, 1], margins=[0, 0, 0, 0], transparent=True)
            self.draw_tpavetext('{dia} Rate Scans{b}'.format(dia=self.DUTName, b=bias_str), lm, 1, 0, 1, font=62, align=13, size=.5, margin=0)
            get_last_canvas().cd()

    def draw_scaled_rate_scans(self, irr=False, y_range=.07, pad_height=.18, scale=1, avrg=False):
        data = zip(self.get_pulse_heights(avrg=avrg), self.get_fluxes(avrg=avrg), get_color_gradient(self.NPlans))
        title_height = pad_height / 2 if self.Title else .03  # half of a pad for title
        c_height = (self.NPlans + .5) * pad_height + title_height  # half of a pad for the x-axis
        c_width = 1.3 * pad_height / .2  # keep aspect ratio for standard pad_height
        c = self.make_canvas(name='csrc', x=c_width, y=c_height, transp=True, logx=True, gridy=True)
        lm, rm, x0, size = .07, .02, .08, .22

        self.draw_title_pad(title_height, x0, lm, c_height)
        self.draw_tpad('p1', 'p1', pos=[0, 0, x0, 1], margins=[0, 0, 0, 0], transparent=True)           # info pad
        self.draw_tpavetext('Scaled Pulse Height', 0, 1, 0, 1, align=22, size=.5, angle=90, margin=0)   # y-axis title
        c.cd()

        for i, (ph, flux, color) in enumerate(data):
            c.cd()
            y0, y1 = [(c_height - title_height - pad_height * (i + j)) / c_height for j in [1, 0]]
            p = self.draw_tpad('p{i}'.format(i=i + 3), '', pos=[x0, y0, 1, y1], margins=[lm, rm, 0, 0], logx=True, gridy=True, gridx=True)
            g = self.make_tgrapherrors('gsph{}'.format(i), '', x=flux, y=ph)
            scale_graph(g, val=scale) if scale else do_nothing()
            format_histo(g, title=' ', color=color, x_range=Bins().FluxRange, y_range=[1 - y_range, 1 + y_range], marker=markers(i), lab_size=size, ndivy=505, markersize=1.5, x_ticks=.15)
            self.draw_histo(g, draw_opt='ap', canvas=p)
            self.draw_legend(i, g, irr, rm)
            c.cd()

        self.draw_tpad('p2', pos=[x0, 0, 1, pad_height / 2 / c_height], margins=[lm, rm, 0, 0], transparent=True)  # x-axis pad
        self.draw_x_axis(1, lm, 1 - rm, 'Flux [kHz/cm^{2}]', Bins().FluxRange, opt='', log=True, tick_size=0, lab_size=size * 2, tit_size=size * 2, off=1.1)
        self.save_plots('ScaledDiaScans{dia}'.format(dia=make_dia_str(self.DUTName)))

    def draw_scaled_distribution(self, excluded=None):
        values = concatenate(([vals / mean_sigma(vals)[0] for i, vals in enumerate(self.get_pulse_heights()) if i not in [excluded]]))
        h = TH1F('hsd', 'Scaled Pulse Height Distribution', 40, .9, 1.1)
        h.FillN(values.size, array([v.n for v in values], 'd'), full(values.size, 1, 'd'))
        self.format_statbox(all_stat=1)
        format_histo(h, x_tit='Scaled Pulse Height', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor)
        self.draw_histo(h, lm=.12)
        return values

    def make_plots(self, name, f, irr_pad=None, canvas=None, **kwargs):
        for sel in self.Info:
            self.info('Creating {} Plots for {}'.format(name, sel.TCString))
            self.get_rp_values(sel, f, **kwargs)
            self.draw_irradiation(make_irr_string(sel.Irradiation), irr_pad, left=False) if irr_pad is not None else do_nothing()
            self.save_plots('{}{}_{}_{}'.format(name, sel.TCString, sel.RunPlan, sel.DUTNr), canvas=get_object(canvas))

    def make_pulse_height_plots(self, y_range=None):
        self.make_plots('PH', AnalysisCollection.draw_pulse_heights, show=False, y_range=y_range, irr_pad=get_object('p1'))

    def make_current_plots(self, c_range=None):
        self.make_plots('Currents', AnalysisCollection.draw_currents, show=False, c_range=c_range, draw_opt='al')

    def make_current_flux_plots(self, c_range=None):
        self.make_plots('CF', AnalysisCollection.draw_currents, show=False, c_range=c_range, draw_opt='al', with_flux=True, canvas='cc')

    def draw_currents(self, align=False, show=True):
        mg = TMultiGraph('mgc', 'Leakage Current vs. Flux')
        legend = self.make_legend(nentries=len(self.RunSelections), w=.4, x2=.52)
        currents = [c.values() for c in self.get_currents()]
        fluxes = [f.values() for f in self.get_fluxes()]
        for i, (x, y) in enumerate(zip(fluxes, currents)):
            g = self.make_tgrapherrors('gf', 'gf', x=x, y=y)
            if align:
                fit = g.Fit('pol1', 'qs0')
                g = self.make_tgrapherrors('gc{}'.format(i), '', y=array(currents[i]) - fit.Parameter(0) + .1, x=fluxes[i])
            format_histo(g, color=self.get_color())
            legend.AddEntry(g, '{tc} - {hv}'.format(tc=self.Info[i].TCString, hv=self.get_rp_values(self.Info[i], AnalysisCollection.get_hv_name, load_tree=False)), 'pl')
            mg.Add(g, 'p')
        format_histo(mg, draw_first=True, y_tit='Current [nA]', x_tit='Flux [kHz/cm^{2}]', y_range=[.1, array(currents).max().n * 2], x_range=Bins().FluxRange)
        self.save_histo(mg, 'CurrentFlux{}'.format(self.Name), draw_opt='a', logx=True, logy=True, leg=legend, bm=.17, show=show)

    def get_titles(self, irr=False):
        if len(set(self.get_dut_names())) > 1:
            return self.get_dut_names()
        tits = self.get_irradiations() if irr else [make_tc_str(tc) for tc in self.TestCampaigns]
        if any(['rand' in word for word in self.get_run_types()]):
            for i, sel in enumerate(self.Info):
                tits[i] += ' (random)' if 'rand' in sel.Type.lower() else '         '
        return tits

    def draw_legend(self, ind, gr, irr, rm):
        add_bias = len(set(self.get_bias_voltages())) > 1
        tits = self.get_titles(irr)
        biases = [make_bias_str(bias) for bias in self.get_bias_voltages()] if add_bias else [''] * len(tits)
        x1 = 1 - max([(12 if irr else len(tit)) + len(bias) for tit, bias in zip(tits, biases)]) * .022
        legend = self.make_legend(x1, 1, x2=1 - rm, nentries=1.2, scale=5)
        legend.AddEntry(gr, tits[ind], 'pe')
        if add_bias:
            legend.SetNColumns(2)
            legend.AddEntry('', biases[ind], '')
        legend.Draw()

    def set_bin_labels(self, h):
        for i, sel in enumerate(self.Info):
            h.GetXaxis().SetBinLabel(h.GetXaxis().FindBin(i), '{} - {}'.format(make_tc_str(sel.TCString, 0), sel.RunPlan))

    def draw_mean_pedestals(self, redo=False, show=True):
        y = array([make_ufloat(mean(tc_values, axis=1)) for tc_values in self.get_pedestals(redo)])
        g = self.make_tgrapherrors('gmp', 'Mean Pedestals', x=arange(y.size), y=y)
        format_histo(g, x_tit='Run Plan', y_tit='Pulse Height [mV]', y_off=1.2, x_range=increased_range([0, y.size - 1], .3, .3), x_off=2.5)
        self.set_bin_labels(g)
        self.save_histo(g, 'PedestalMeans{}'.format(self.Name.title().replace('-', '').replace('_', '')), show, draw_opt='ap', bm=.2, x=1.5, y=.75, gridy=True)

    def draw_means(self, y_range=None, show=True):
        y = array([make_ufloat(mean_sigma([dic['ph'] for dic in ph.itervalues()])) for ph in self.get_pulse_heights()])
        g = self.make_tgrapherrors('gms', 'Pulse Height Evolution', x=arange(y.size), y=y)
        format_histo(g, x_tit='Run Plan', y_tit='Mean Pulse Height [mV]', y_off=1.2, x_range=increased_range([0, y.size - 1], .3, .3), x_off=2.5, y_range=y_range)
        self.set_bin_labels(g)
        self.save_histo(g, 'Means{}'.format(self.Name.title().replace('-', '').replace('_', '')), show, draw_opt='ap', bm=.2, x=1.5, y=.75, gridy=True)

    def draw_sigmas(self, y_range=None, show=True):
        mean_sigmas = [mean_sigma([dic['ph'] for dic in ph.itervalues()]) for ph in self.get_pulse_heights()]
        y = array([make_ufloat((0, 100 * s / m)) for m, s in mean_sigmas])
        g = self.make_tgrapherrors('gms', 'Pulse Height STD Evolution', x=arange(y.size), y=y)
        format_histo(g, x_tit='Run Plan', y_tit='Pulse Height Standard Deviation [%]', y_off=1.2, x_range=increased_range([0, y.size - 1], .3, .3), x_off=2.5, y_range=y_range)
        self.set_bin_labels(g)
        self.save_histo(g, 'Means{}'.format(self.Name.title().replace('-', '').replace('_', '')), show, draw_opt='ap', bm=.2, x=1.5, y=.75, gridy=True)

    def draw_uniformity(self, arg=2, redo=False, low=False, high=False):
        y_values = [v[arg] for v in self.get_mean_uniformities(redo=redo, low=low, high=high)]
        x_values = [make_ufloat((v, v * .2)) for v in self.get_irradiations(string=False)]
        g = self.make_tgrapherrors('gu', 'Uniformity', x=x_values, y=y_values)
        format_histo(g, x_tit='Irradiation [n/cm^{2}]', y_tit='FWHM/MPV', y_off=1.3)
        self.draw_histo(g, draw_opt='ap')

    def draw_peak_flux(self, show=True):
        leg = self.make_legend(x2=.45, w=.3)
        mg = TMultiGraph('mgph', 'Peak Flux vs. FAST-OR Flux')
        values_list = self.get_values(AnalysisCollection.get_peak_flux, PickleInfo('Peaks', 'Flux'))
        flux_list = self.get_fluxes()
        for sel, values, fluxes in zip(self.Info, values_list, flux_list):
            g = self.make_tgrapherrors('g{}'.format(sel.RunPlan), '', x=fluxes, y=values)
            format_histo(g, color=self.get_color())
            leg.AddEntry(g, '{} @ {:+1.0f}V'.format(sel.DUTName, sel.Bias), 'pl')
            mg.Add(g, 'p')
        x, y = concatenate(flux_list), concatenate(values_list)
        x_range, y_range = [.5 * min(x).n, 1.2 * max(x).n], [.5 * min(y).n, 1.2 * max(y).n]
        format_histo(mg, draw_first=True, y_tit='Peak Flux [kHz/cm^{2}] ', x_tit='FAST-OR Flux [kHz/cm^{2}]', x_range=x_range, y_off=1.8, y_range=y_range)
        self.save_histo(mg, 'PeakFluxes{}'.format(self.Name), draw_opt='a', leg=leg, show=show, lm=.13)
    # endregion DRAWING
    # ----------------------------------------


class SelectionInfo:
    def __init__(self, sel):
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

    def __str__(self):
        return 'Selection instance: {} {} {}'.format(self.TCString, self.RunPlan, self.DUTName)

    def __repr__(self):
        return self.__str__()

    def __call__(self):
        return [self.TCString, self.RunPlan, self.DUTName, self.DUTNr, '{:03d}-{:03d}'.format(self.Runs[0], self.Runs[-1]), '{:+4.0f}V'.format(self.Bias), self.Type, self.Irradiation]


class PickleInfo:
    def __init__(self, sub_dir, name, suf=None):
        self.SubDir = sub_dir
        self.Name = name
        self.Suffix = '' if suf is None else str(suf)


if __name__ == '__main__':

    aparser = ArgumentParser()
    aparser.add_argument('sel', nargs='?', default=None, help='name of the selection')
    aparser.add_argument('-d', nargs='?', default=None, help='dut')
    aparser.add_argument('-tc', nargs='?', default=None, help='test campaign')
    aparser.add_argument('-v', '--verbose', action='store_false')
    aparser.add_argument('-p', action='store_true', help='print analysis strings')
    aparser.add_argument('-r', action='store_true', help='redo')
    aparser.add_argument('-s', action='store_true', help='activate show single selection')
    aparser.add_argument('-sa', action='store_true', help='active show all selections')
    pargs = aparser.parse_args()

    z = DiaScans(pargs.sel, pargs.verbose)
    if pargs.p:
        print(z.get_all_ana_strings(pargs.d, pargs.tc, pargs.r))
    if pargs.s:
        z.show_selection()
    if pargs.sa:
        z.show_selections()
