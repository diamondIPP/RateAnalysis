# --------------------------------------------------------
#       DIAMOND RATE SCANS
# created on June 24th 2016 by M. Reichmann
# --------------------------------------------------------

from AnalysisCollection import AnalysisCollection
from Elementary import Elementary
from operator import itemgetter
from ConfigParser import ConfigParser, NoOptionError
from Utils import *
from argparse import ArgumentParser
from RunSelection import RunSelection
from Run import Run
from json import load, dump, loads
from collections import OrderedDict, Counter
from ROOT import TMultiGraph, TGraphErrors, kRed, kOrange, kBlue, kGreen, kCyan, kViolet, kPink, kYellow, gStyle, TF1, TH2F, TH1F, TGraph2DErrors
import pickle
from numpy import sqrt
from re import split
from uncertainties.unumpy import uarray


class DiaScans(Elementary):
    def __init__(self, selection, verbose=False):
        Elementary.__init__(self, verbose=verbose)
        self.Selection = []
        self.Name = None
        self.scale_factors = OrderedDict()

        self.Parser = self.load_diamond_parser()

        # information
        self.DiamondName = None
        self.TestCampaigns = self.load_testcampaigns(None)
        self.RunInfos = self.load_runinfos()
        self.AllRunPlans = self.load_all_runplans()
        self.RunPlans = self.find_diamond_runplans()
        self.Bias = None
        self.set_save_directory('Results/')
        self.save_dir = ''
        self.UsedColors = self.init_colors()

        # run plan selection
        self.Path = join(self.Dir, self.MainConfigParser.get('MISC', 'runplan_selection_file'))
        self.Selections = self.load_selections()
        self.Selection = self.load_selection()
        self.RunSelections = None

        # Save
        self.ROOTObjects = []

        self.set_selection(selection)

    @staticmethod
    def init_colors():
        return {'pulser': 0, 'neg': 0, 'pos': 0}

    def get_next_color(self, bias, pulser):
        colors = {'neg': [1, kBlue, kGreen + 2], 'pos': [kRed, kOrange, kPink], 'pulser': [kGreen, kViolet, kBlue, kYellow + 1, kCyan]}
        n_col = 3 if not pulser else 5

        def get_col(tag):
            color = colors[tag][self.UsedColors[tag] % n_col]
            self.UsedColors[tag] += 1
            return color

        if pulser:
            return get_col('pulser')
        elif bias < 0:
            return get_col('neg')
        else:
            return get_col('pos')

    # ==========================================================================
    # region INIT

    def load_selections(self):
        f = open(self.Path)
        selections = load(f, object_pairs_hook=OrderedDict)
        f.close()
        return selections

    def load_selection(self):
        return self.Selections[self.DiamondName] if self.DiamondName in self.Selections else {}

    def load_diamond_parser(self):
        parser = ConfigParser()
        parser.read('{0}/Configuration/DiamondAliases.cfg'.format(self.get_program_dir()))
        return parser

    def load_diamond(self, name):
        dia = name
        if 'all' in dia:
            return 'All'
        if 'test' in dia:
            return 'Test'
        if name.lower() not in self.Parser.options('ALIASES'):
            dia = split('[-_]', name)[-1]
        if dia.lower() not in self.Parser.options('ALIASES'):
            dia = '-'.join(split('[-_]', name)[:-1])
        if dia.lower() not in self.Parser.options('ALIASES'):
            dia = dia.split('-')[-1]
        try:
            return self.Parser.get('ALIASES', dia)
        except NoOptionError:
            log_warning('{0} is not a known diamond name! Please choose one from \n{1}'.format(dia, self.Parser.options('ALIASES')))

    def load_testcampaigns(self, tcs):
        if tcs is None:
            return ['201508', '201510']
        valid_tcs = self.find_test_campaigns()
        tcs = [tcs] if type(tcs) is not list else tcs
        if not all(tc in valid_tcs for tc in tcs):
            log_warning('You entered and invalid test campaign! Aborting!')
            exit()
        else:
            return tcs

    def load_all_runplans(self):
        runplan_path = join(self.Dir, self.MainConfigParser.get('MISC', 'runplan_file'))
        f = open(runplan_path, 'r')
        runplans = load(f)
        f.close()
        return runplans

    def load_runinfos(self):
        run_infos = {}
        for tc in self.TestCampaigns:
            self.set_test_campaign(tc)
            self.TCDir = self.generate_tc_directory()
            file_path = self.load_run_info_path()
            f = open(file_path)
            run_infos[tc] = load(f)
            f.close()
        return run_infos

    def find_diamond_runplans(self):
        runplans = {}
        for tc in self.TestCampaigns:
            runplans[tc] = {}
            for rp, dic in self.AllRunPlans[tc].iteritems():
                runs = dic['runs']
                for ch in [1, 2]:
                    if all(self.DiamondName == self.load_diamond(self.RunInfos[tc][str(run)]['dia{0}'.format(ch)]) for run in runs):
                        bias = self.RunInfos[tc][str(runs[0])]['dia{0}hv'.format(ch)]
                        if all(self.RunInfos[tc][str(run)]['dia{0}hv'.format(ch)] == bias for run in runs):
                            if bias not in runplans[tc]:
                                runplans[tc][bias] = {}
                            runplans[tc][bias][rp] = ch
        return runplans

    # endregion

    def get_diamond_names(self):
        return [sel.SelectedDiamond for sel in self.RunSelections]

    def get_run_types(self):
        return [sel.SelectedType.lower() for sel in self.RunSelections]

    def get_irradiations(self):
        return [make_irr_string(sel.get_irradiation()) for sel in self.RunSelections]

    def get_bias_voltages(self):
        return [sel.SelectedBias for sel in self.RunSelections]

    def set_selection(self, key=None):
        if key is None:
            key = self.DiamondName
        if key not in self.Selections.keys():
            log_warning('"{sel} does not exist in:'.format(sel=key))
            for sel in sorted(self.Selections.keys()):
                print sel
            return
        self.log_info('Set Selection {0}'.format(key))
        self.DiamondName = self.load_diamond(key)
        self.Selection = self.Selections[key]
        self.TestCampaigns = list(set(self.Selection.keys()))
        self.load_run_selections(redo=True)
        self.Name = key

    # ==========================================================================
    # region SHOW

    def show_runplans(self):
        for tc, vals in self.RunPlans.iteritems():
            print_small_banner(tc.rjust(15))
            for bias, val1s in vals.iteritems():
                print '{0} V:'.format(int(bias))
                for rp, ch in val1s.iteritems():
                    print ' ', rp.ljust(5), ch
                print

    def show_selections(self):
        print 'The following selections exists:'
        for key, sel in self.Selections.items():
            print '* "{key}": {sel}'.format(key=key, sel=sel)

    def show_selection_names(self):
        for key in self.Selections.iterkeys():
            print key

    def show_selection(self):
        """ Gives detailed information about the chosen selection """
        if self.Selection:
            header = ['Campaign', 'RunPlan', 'Diamond', 'Nr', 'Runs'.ljust(7), 'Voltage', 'Type'.ljust(11)]
            rows = []
            for sel in self.RunSelections:
                runs = sel.get_selected_runs()
                row = [sel.TCString.ljust(8)]                                               # Campaign
                row += [sel.SelectedRunplan.rjust(7)]                                       # Run Plan
                row += [sel.SelectedDiamond.rjust(7)]                                       # Diamond Name
                row += [str(sel.SelectedDiamondNr).rjust(2)]                                # Diamond Number
                row += ['{0}-{1}'.format(str(runs[0]).zfill(3), str(runs[-1]).zfill(3))]    # Selected Runs
                row += ['{0:+4.0f}V'.format(sel.SelectedBias).rjust(7)]                     # Bias Voltage
                row += [sel.SelectedType.ljust(11)]                                         # Run Plan Type
                rows.append(row)
            print_table(rows, header)
        else:
            log_warning('Selection is empty!')

    def get_runs(self, rp, tc):
        return self.AllRunPlans[tc][rp]['runs']

    def show_all_runplans(self):
        for tc in self.TestCampaigns:
            print_small_banner(tc)
            for rp, runs in sorted(self.AllRunPlans[tc]['rate_scan'].iteritems()):
                dias = [self.load_diamond(self.RunInfos[tc][str(runs[0])]['dia{0}'.format(ch)]) for ch in [1, 2]]
                print rp.ljust(5), '{0}-{1}'.format(str(runs[0]).zfill(3), str(runs[-1]).zfill(3)), dias[0].ljust(11), dias[1].ljust(11)

    # endregion

    # ==========================================================================
    # region SELECTION
    def select_runplan(self, runplan, ch=1, testcampaign=None):
        rp = self.make_runplan_string(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        if rp in self.AllRunPlans[tc]['rate_scan']:
            if tc not in self.Selection:
                self.Selection[tc] = {}
            self.Selection[tc][rp] = ch
        else:
            log_warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def unselect_runplan(self, runplan, testcampaign=None):
        rp = self.make_runplan_string(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        try:
            self.Selection[tc].pop(rp)
        except KeyError:
            log_warning('The runplan {0} does not exist in {1}!'.format(rp, tc))

    def select_runplans_by_bias(self, value):
        self.clear_selection()
        for tc, vals in self.RunPlans.iteritems():
            for bias, rps in vals.iteritems():
                if abs(bias) == value:
                    for rp, ch in rps.iteritems():
                        self.select_runplan(rp, ch, tc)

    def clear_selection(self):
        self.Selection = {}

    def save_selection(self, name):
        file_path = self.get_program_dir() + self.MainConfigParser.get('MISC', 'runplan_selection_file')
        f = open(file_path, 'r+')
        selections = load(f)
        if self.Selection:
            selections[name] = self.Selection
        else:
            log_warning('Selection is empty!')
        f.seek(0)
        dump(selections, f, indent=2, sort_keys=True)
        f.truncate()
        f.close()
        self.Selections = selections

    # endregion

    @staticmethod
    def make_runplan_string(nr):
        nr = str(nr)
        return nr.zfill(2) if len(nr) <= 2 else nr.zfill(4)

    @staticmethod
    def get_ph_below_flux(mg, flux=80, keys=None):
        keys = ['gFullError', 'data'] if keys is None else keys
        try:
            for g in mg.GetListOfGraphs():
                if g.GetName() in keys:
                    return DiaScans.get_ph_below_flux(g, flux)
            log_warning('cannot find correct data {0}'.format([g.GetName() for g in mg.GetListOfGraphs()]))
        except AttributeError:
            g = mg
            n, x, y, ex, ey = get_graph_data(g)
            ey = [0] * n if not hasattr(g, 'GetEY') else ey
            xy = zip(x, y)
            xy_filtered = filter(lambda x1: x1[0] < flux, xy)
            mean = calc_mean(map(itemgetter(1), xy_filtered))
            ymin = min(map(lambda val, err: val - err, y, ey))
            ymax = max(map(lambda val, err: val + err, y, ey))
            return mean + (ymin, ymax)

    def draw_hysteresis_graph(self, scans, limits=None):
        limits = [0, 30, 80, 250, 800, 2000, 4000, 1e10] if limits is None else limits

        def get_bucket(x_val, lim):
            for j in xrange(len(lim) - 1):
                if lim[j] < x_val < lim[j + 1]:
                    return j
            raise Exception()

        def calc_weighted_diff(x_val, y_val):
            x_val = list(x_val)
            y_val = list(y_val)
            retval = [0] * 4
            retval[0] = ((x_val[0][0] + y_val[0][0]) / 2)
            retval[1] = sqrt(x_val[0][1] ** 2 + y_val[0][1] ** 2) / 2
            retval[2] = ((x_val[1][0] - y_val[1][0]) / 2)
            retval[3] = sqrt(x_val[1][1] ** 2 + y_val[1][1] ** 2) / 2
            return retval

        n = len(limits) - 1
        keys = scans.keys()
        data = [{}] * 4
        values = None
        for key, d in scans.iteritems():
            values = [[]] * n
            for x, y, ex, ey in zip(d['x'], d['y'], d['ex'], d['ey']):
                i = get_bucket(x, limits)
                values[i].append([(x, ex), (y, ey)])
            for v in values:
                vv = zip(*v)
                xx = zip(*vv[0])
                yy = zip(*vv[1])
                x = calc_mean(xx[0])
                y = calc_weighted_mean(yy[0], yy[1])
                i = values.index(v)
                data[i][key] = (x, y)

        g = TGraphErrors(len(values))
        g.SetName('gHysteresis')
        tit = 'Flux_{{{k0}}} #minus Flux_{{{k1}}}'.format(k0=keys[0], k1=keys[1])
        g.SetTitle(tit + ';flux[kHz/cm^{2}];' + tit + ' [au]')
        for d in data:
            x = d[keys[0]]
            y = d[keys[1]]
            val = calc_weighted_diff(x, y)
            i = data.index(d)
            g.SetPoint(i, val[0], val[2])
            g.SetPointError(i, val[1], val[3])
        self.format_histo(g, y_off=1.4, x_off=1.3, x_tit='Flux [kHz/cm^{2}]', y_tit='hysteresis: ' + tit + ' [au]', draw_first=True)
        pname = 'Hysteresis_{key}'.format(key=self.Name)
        self.draw_histo(g, pname, lm=.14, draw_opt='ALP', logx=True, gridy=True)
        g.GetYaxis().SetNdivisions(509)
        fit = TF1('fit', 'pol0', 0, g.GetXaxis().GetXmax())
        gStyle.SetOptFit(11)
        g.Fit(fit, 'QS')
        self.save_plots(pname)

    def load_run_selections(self, redo=False):
        if self.RunSelections is not None and not redo:
            return self.RunSelections
        run_selections = []
        for tc, rps in self.Selection.iteritems():
            for rp, ch in rps.iteritems():
                sel = RunSelection(tc)
                sel.select_runs_from_runplan(rp, ch=ch)
                self.log_info('Loaded runplan {rp} of testcampaign {tc} and ch {ch} ({dia})'.format(rp=rp.rjust(4), tc=make_tc_str(tc), ch=ch, dia=sel.SelectedDiamond))
                run_selections.append(sel)
        self.RunSelections = run_selections
        return run_selections

    def draw_collimator_settings(self, show=True):
        h = TH2F('h_cs', 'Collimator Settings', 125, 50, 300, 75, 0, 150)
        for tc in self.TestCampaigns:
            for run, data in self.RunInfos[tc].iteritems():
                try:
                    h.Fill(data['fs11'], data['fs13'])
                except KeyError:
                    pass
        self.format_histo(h, x_tit='fs11', y_tit='fsh13', y_off=1.3, stats=0, z_off=1.1, z_tit='Number of Entries', z_range=[0, 80])
        self.save_histo(h, 'CollimatorSettings', show, draw_opt='colz', lm=.12, rm=.16)

    def draw_flux_vs_collimators(self, show=True):
        gr = TGraph2DErrors()
        gr.SetNameTitle('gr_fc', 'Flux Vs. Collimators')
        col_settings = Counter([(data['fs11'], data['fs13']) for tc in self.TestCampaigns for data in self.RunInfos[tc].itervalues() if 'fs11' in data and data['fs11'] > 0])
        i = 0
        for col, nr in col_settings.iteritems():
            if nr > 10:
                flux_fit = z.draw_flux_distribution(col[0], col[1], show=False)
                if flux_fit is not None:
                    gr.SetPoint(i, col[0], col[1], flux_fit.Parameter(1))
                    gr.SetPointError(i, 0, 0, flux_fit.Parameter(2))
                    i += 1
        self.draw_histo(gr, 'FluxVsCollimators', show, draw_opt='surf1', lm=.15, phi=17, theta=35)
        self.format_histo(gr, x_tit='fs11', x_off=1.3, y_tit='fsh13', y_off=1.9, stats=0, z_off=1.9, z_tit='Flux kHz/cm^{2}', markersize=2, y_range=[0, 130])
        self.save_plots('FluxVsCollimators', show=show, prnt=False)
        h = gr.Clone()
        h.Draw('samep')
        self.ROOTObjects.append(h)
        self.save_plots('FluxVsCollimators', show=show)

    def draw_flux_variations(self, show=True, rel_sigma=False):
        gr = self.make_tgrapherrors('gr_fd', 'Flux Deviations')
        col_settings = Counter([(data['fs11'], data['fs13']) for tc in self.TestCampaigns for data in self.RunInfos[tc].itervalues() if 'fs11' in data and data['fs11'] > 0])
        i = 0
        for col, nr in sorted(col_settings.iteritems()):
            if nr > 30:
                flux_fit = z.draw_flux_distribution(col[0], col[1], show=False)
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
        self.format_histo(gr, x_tit='Mean Flux [au]', y_tit='{0}Sigma [au]'.format('Relative ' if rel_sigma else ''), y_off=1.3)
        self.draw_histo(gr, 'FluxVariations', show, draw_opt='alp', logx=True, logy=not rel_sigma, lm=.12)
        gr.GetXaxis().SetLimits(gr.GetX()[0] / 2, gr.GetX()[gr.GetN() - 1] * 4)
        self.save_plots('FluxVariations{0}'.format('Rel' if rel_sigma else ''), show=show)

    def draw_flux_distribution(self, fs11, fsh13, tc=None, do_fit=True, show=True, run_thr=None):
        values = []
        for tc in self.TestCampaigns if tc is None else [tc]:
            for run, data in sorted(self.RunInfos[tc].iteritems()):
                info_run = Run(run_number=run, test_campaign=tc, tree=False)
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
        self.set_root_output(False)
        h = TH1F('h_fd', 'Flux Distribution for {0}/{1}'.format(fs11, fsh13), int(sqrt(len(values))) + 5, min(values) - .2 * spread, max(values) + .2 * spread)
        for val in values:
            h.Fill(val)
        set_statbox(only_fit=True, w=.25) if do_fit else do_nothing()
        self.format_histo(h, x_tit='Flux in kHz/cm^{2}', y_tit='Number of Entries', y_off=1.3, stats=0 if not do_fit else 1)
        self.draw_histo(h, '', show)
        fit = None
        if do_fit:
            h.SetName('Fit Results')
            self.set_root_output(show)
            fit = h.Fit('gaus', 'qs')
        self.save_plots('FluxDistribution{0}_{1}'.format(int(fs11), int(fsh13)), show=show)
        return fit

    def draw_dia_rate_scans(self):
        run_selections = self.load_run_selections()
        biases = self.get_bias_voltages()
        bias_str = ' at {bias} V'.format(bias=biases[0]) if len(biases) == 1 else ''
        mg = TMultiGraph('mg_ph', '{dia} Rate Scans{b};Flux [kHz/cm^{{2}}]; pulse height [au]'.format(dia=self.DiamondName, b=bias_str))
        legend = self.make_legend(.75, .4, nentries=4, felix=True)
        legend.SetNColumns(2) if len(biases) > 1 else do_nothing()
        colors = [4, 419, 2, 800, 3]
        # tits = [make_irr_string(v, p) for v, p in [(0, 0), (5, 14), (1.5, 15)]]
        tits = self.get_irradiations()
        for i, sel in enumerate(run_selections):
            path = self.make_pickle_path('Ph_fit', 'PulseHeights', sel.SelectedRunplan, self.DiamondName, 10000, sel.TESTCAMPAIGN)
            try:
                f = open(path, 'r')
                mg_ph_ana = pickle.load(f)
                f.close()
            except IOError:
                print 'Did not find', path
                Elementary(sel.generate_tc_str())
                t = load_root_files(sel, load=True)
                ana = AnalysisCollection(sel, threads=t)
                mg_ph_ana = ana.draw_pulse_heights(show=False)
                ana.close_files()
            for g in mg_ph_ana.GetListOfGraphs():
                self.format_histo(g, color=colors[i], markersize=1.5, lw=2)
                if g.GetName() == 'gFirst':
                    self.format_histo(g, color=1, marker=26, markersize=2)
                elif g.GetName() == 'gLast':
                    self.format_histo(g, color=1, marker=23, markersize=2)
            legend.AddEntry(mg_ph_ana.GetListOfGraphs()[0], tits[i], 'lp')
            legend.AddEntry(0, get_bias_root_string(sel.SelectedBias), '') if len(biases) > 1 else do_nothing()
            mg.Add(mg_ph_ana)
        x_vals = sorted([gr.GetX()[i] for gr in mg.GetListOfGraphs() for i in xrange(gr.GetN())])
        y_vals = sorted([gr.GetY()[i] for gr in mg.GetListOfGraphs() for i in xrange(gr.GetN())])
        self.format_histo(mg, draw_first=True, y_tit='Pulse Height [au]', y_range=[0, y_vals[-1] * 1.1], tit_size=.05, lab_size=.05, y_off=.91, x_off=1.2)
        mg.GetXaxis().SetLimits(x_vals[0] * 0.8, x_vals[-1] * 3)
        self.save_histo(mg, 'DiaScans{dia}'.format(dia=make_dia_str(self.DiamondName)), draw_opt='a', logx=True, l=legend, x_fac=1.6, lm=.092, bm=.12, gridy=True)

    def draw_scaled_rate_scans(self, irr=False, y_range=.15):
        run_selections = self.load_run_selections()
        biases = set(self.get_bias_voltages())
        bias_str = ' at {b}'.format(b=make_bias_str(biases.pop())) if len(biases) == 1 else ''
        colors = get_color_gradient(len(run_selections))
        graphs = self.get_pulse_height_graphs()
        x_vals = sorted([g.GetX()[i] for g in graphs for i in xrange(g.GetN())])
        y_range = [1 - y_range, 1 + y_range]

        c = self.make_canvas(x=1.6, y=.2 * (len(graphs) + 1), transp=True, logx=True, gridy=True)
        lm, rm = .03, .02
        lp = .04
        pad_height = .1 / (.2 * (len(graphs) + 1))
        self.draw_tpad('p0', 'p0', pos=[lp, 1 - pad_height, 1, 1], margins=[lm, rm, 0, .5], transparent=True)           # title pad
        self.draw_tlatex(lm, .5, '{dia} Rate Scans{b}'.format(dia=self.DiamondName, b=bias_str), size=.5, align=12)     # title
        c.cd()
        self.draw_tpad('p1', 'p1', pos=[0, 0, lp, 1], margins=[0, 0, 0, 0], transparent=True)                           # info pad
        self.draw_tlatex(.5, .5, '#font[42]{Scaled Pulse Height}', size=.7, align=22, angle=90)                         # title
        c.cd()

        for i, g in enumerate(graphs):
            last = i == len(self.RunSelections) - 1
            y0 = 0 if last else 1 - (pad_height * (2 * i + 3))
            y1 = 1 - (pad_height * (2 * i + 1))
            self.draw_tpad('p{i}'.format(i=i + 2), 'p{i}'.format(i=i + 2), pos=[lp, y0, 1, y1], margins=[lm, rm, .33 if last else 0, 0], logx=True, gridy=True, gridx=True)
            if last:
                self.format_histo(g, title=' ', color=colors[i], x_range=[x_vals[0] * 0.8, x_vals[-1] * 3], y_range=y_range, marker=markers(i), lab_size=.15 * 2 / 3, ndivy=505,
                                  markersize=1.5, x_tit='Flux [kHz/cm^{2}]', tit_size=.15 * 2 / 3, tick_size=.15 * 2 / 3)
            else:
                self.format_histo(g, title=' ', color=colors[i], x_range=[x_vals[0] * 0.8, x_vals[-1] * 3], y_range=y_range, marker=markers(i), lab_size=.15, ndivy=505,
                                  markersize=1.5, tick_size=.15)
            g.GetXaxis().SetLimits(x_vals[0] * 0.8, x_vals[-1] * 3)
            g.Draw('ap')
            self.draw_legend(i, g, irr, rm)
            c.cd()

        self.ROOTObjects.append(graphs)
        self.save_plots('ScaledDiaScans{dia}'.format(dia=make_dia_str(self.DiamondName)))

    def get_titles(self, irr=False):
        if len(set(self.get_diamond_names())) > 1:
            return self.get_diamond_names()
        if irr:
            tits = self.get_irradiations()
        else:
            tits = [make_tc_str(self.RunSelections[i].TCString) for i in xrange(len(self.RunSelections))]
        if 'rand' in self.get_run_types():
            for i, sel in enumerate(self.RunSelections):
                tits[i] += ' (random)' if 'rand' in sel.Type.lower() else '         '
        return tits

    def draw_legend(self, ind, gr, irr, rm):
        add_bias = len(set(self.get_bias_voltages())) > 1
        tits = self.get_titles(irr)
        biases = [make_bias_str(bias) for bias in self.get_bias_voltages()] if add_bias else [''] * len(tits)
        x1 = 1 - max([(12 if irr else len(tit)) + len(bias) for tit, bias in zip(tits, biases)]) * .016
        legend = self.make_legend(x1, 1, x2=1 - rm, nentries=1, scale=5 * (2 / 3. if ind == len(tits) - 1 else 1))
        legend.AddEntry(gr, tits[ind], 'pe')
        if add_bias:
            legend.SetNColumns(2)
            legend.AddEntry('', biases[ind], '')
        legend.Draw()
        self.ROOTObjects.append(legend)

    def get_pulse_height_graphs(self, scale=1, redo=False):
        run_selections = self.load_run_selections()
        graphs = []
        for i, sel in enumerate(run_selections):
            path = self.make_pickle_path('Ph_fit', 'PhVals', sel.SelectedRunplan, sel.SelectedDiamond, 10000, sel.TCString)
            try:
                if redo:
                    raise IOError
                f = open(path, 'r')
                phs = pickle.load(f)
                f.close()
            except IOError:
                print 'Did not find', path
                Elementary(sel.generate_tc_str())
                print
                t = load_root_files(sel, load=True)
                ana = AnalysisCollection(sel, t)
                phs = ana.get_pulse_heights(redo=redo)
            values, errors = self.scale_to(phs, scale)
            fluxes = [ph['flux'] for ph in phs.itervalues()]
            g = self.make_tgrapherrors('g{n}'.format(n=i), 'Rate Scans for {n}'.format(n=self.DiamondName))
            for j, (x, val, err) in enumerate(zip(fluxes, values, errors)):
                g.SetPoint(j, x, val)
                g.SetPointError(j, .1 * x, err)
            graphs.append(g)
        return graphs

    @staticmethod
    def scale_to(dic, scale=1):
        fluxes = [d['flux'] for d in dic.itervalues()]
        min_ind = fluxes.index(min(fluxes))
        err_scale = scale / dic.values()[min_ind]['ph'].Parameter(0)
        scale = scale / dic.values()[min_ind]['ph'].Parameter(0) if scale is not None else 1
        values = [d['ph'].Parameter(0) * scale for d in dic.itervalues()]
        errors = [d['ph'].ParError(0) * err_scale for d in dic.itervalues()]
        return values, errors

    def draw_beam_induced_currents(self, show=True):
        parser = ConfigParser()
        parser.read(join(self.Dir, 'Runinfos', 'beam_induced_currents.ini'))
        ymin, ymax, rate, tcs, min_err, max_err, phs = [array(loads(parser.get('MAIN', opt))) for opt in ['min', 'max', 'rate', 'tc', 'min_err', 'max_err', 'pulse_height']]
        ymin, ymax = uarray(ymin, min_err), uarray(ymax, max_err)
        rate = uarray(rate, rate * .1)
        phs = uarray(phs, phs * .02)
        yvals = (ymax - ymin) / rate / phs
        yvals = yvals / yvals[0].n
        irrads = {}
        sel = self.RunSelections[0]
        for tc in tcs:
            sel.set_test_campaign(tc)
            irrads[tc] = float(sel.get_irradiation('II6-B2')) / 1e14
        xvals = [irrads[tc] for tc in tcs]
        g = self.make_tgrapherrors('gbic', 'Beam Induced Currents')
        for i, (x, y) in enumerate(zip(xvals, yvals)):
            g.SetPoint(i, x, y.n)
            g.SetPointError(i, x * .1, y.s)
        self.format_histo(g, y_tit='Beam Induced Current / Pulse Height', x_tit='Irradiation [1e14 n/cm^{2}]')
        self.save_histo(g, 'BeamInducedErrors', draw_opt='alp', show=show)

    def add_selection(self):
        name = raw_input('Enter the name of the selection: ')
        self.Selections[name] = {} if name not in self.Selections else self.Selections[name]
        run_plan = raw_input('Enter test campaign, run plan number, channel: ')
        while run_plan:
            tc, rp, ch = [string.strip(' ') for string in run_plan.split(',')]
            rp = RunSelection.make_runplan_string(rp)
            if tc not in self.Selections[name]:
                self.Selections[name][tc] = {rp: int(ch)}
            else:
                self.Selections[name][tc][rp] = int(ch)
            run_plan = raw_input('Enter test campaign, run plan number, channel: ')
        self.save_selections()

    def save_selections(self):
        f = open(self.Path, 'r+')
        f.seek(0)
        dump(self.Selections, f, indent=2, sort_keys=True)
        f.truncate()
        f.close()


if __name__ == '__main__':
    main_parser = ArgumentParser()
    main_parser.add_argument('sel', nargs='?', default='test')
    main_parser.add_argument('-v', action='store_true')
    args = main_parser.parse_args()
    print_banner('STARTING DIAMOND RATE SCAN COLLECTION OF SELECTION {0}'.format(args.sel))

    Elementary(None, True, get_resolution())
    z = DiaScans(args.sel, verbose=args.v)
