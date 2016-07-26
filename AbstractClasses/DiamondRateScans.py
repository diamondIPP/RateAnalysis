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
from json import load, dump
from collections import OrderedDict
from ROOT import TMultiGraph, TGraphErrors, kRed, kOrange, kBlue, kGreen, kCyan, kViolet, kPink, kYellow, gStyle, TF1, TH2F, TH1F, TGraph2DErrors
import pickle
from numpy import mean, sqrt


class DiaScans(Elementary):
    def __init__(self, diamond=None, testcampaigns=None, verbose=False):
        Elementary.__init__(self, verbose=verbose)
        self.Selection = []
        self.Name = None
        self.scale_factors = OrderedDict()

        self.Parser = self.load_diamond_parser()

        # information
        self.DiamondName = self.load_diamond(diamond)
        self.TestCampaigns = self.load_tcs(testcampaigns)
        self.RunInfos = self.load_runinfos()
        self.AllRunPlans = self.load_all_runplans()
        self.RunPlans = self.find_diamond_runplans()
        self.Bias = None
        self.set_save_directory('Results/')
        self.save_dir = ''
        self.UsedColors = self.init_colors()

        # run plan selection
        self.Selections = self.load_selections()
        self.Selection = self.load_selection()
        self.RunSelections = None

        # Save
        self.ROOTObjects = []
        self.PickleDir = self.get_program_dir() + self.ana_config_parser.get('SAVE', 'pickle_dir')

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
        file_path = self.get_program_dir() + self.MainConfigParser.get('MISC', 'runplan_selection_file')
        f = open(file_path)
        selections = load(f)
        f.close()
        return selections

    def load_selection(self):
        return self.Selections[self.DiamondName] if self.DiamondName in self.Selections else {}

    def load_diamond_parser(self):
        parser = ConfigParser()
        parser.read('{0}/Configuration/DiamondAliases.cfg'.format(self.get_program_dir()))
        return parser

    def load_diamond(self, dia):
        try:
            return self.Parser.get('ALIASES', dia)
        except NoOptionError:
            if dia in [self.Parser.get('ALIASES', a) for a in self.Parser.options('ALIASES')]:
                return dia
            log_warning('{0} is not a known diamond name! Please choose one from \n{1}'.format(dia, self.Parser.options('ALIASES')))
            exit()

    def load_tcs(self, tcs):
        if tcs is None:
            return ['201505', '201508', '201510']
        valid_tcs = self.find_test_campaigns()
        tcs = [tcs] if type(tcs) is not list else tcs
        if not all(tc in valid_tcs for tc in tcs):
            log_warning('You entered and invalid test campaign! Aborting!')
            exit()
        else:
            return tcs

    def load_all_runplans(self):
        runplan_path = self.get_program_dir() + self.MainConfigParser.get('MISC', 'runplan_file')
        f = open(runplan_path, 'r')
        runplans = load(f)
        f.close()
        return runplans

    def load_runinfos(self):
        run_infos = {}
        for tc in self.TestCampaigns:
            self.TESTCAMPAIGN = tc
            parser = self.load_run_configs(None)
            file_path = parser.get('BASIC', 'runinfofile')
            f = open(file_path)
            run_infos[tc] = load(f)
            f.close()
        return run_infos

    def find_diamond_runplans(self):
        runplans = {}
        for tc in self.TestCampaigns:
            runplans[tc] = {}
            for rp, runs in self.AllRunPlans[tc]['rate_scan'].iteritems():
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
        names = []
        for sel, ch in self.RunSelections.iteritems():
            names.append(self.load_diamond(sel.Diamond))
        return list(set(names))

    def get_bias_voltages(self):
        voltages = []
        for sel, ch in self.RunSelections.iteritems():
            voltages.append(sel.SelectedBias)
        return list(set(voltages))

    def set_selection(self, key=None):
        if key is None:
            key = self.DiamondName
        if key not in self.Selections.keys():
            log_warning('"{selection} does not exist in {Selections}'.format(selection=key, Selections=self.Selections))
            return
        self.log_info('Set Selection {0}'.format(key))
        self.Selection = self.Selections[key]
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

    def show_selection(self):
        if self.Selection:
            for tc, rps in sorted(self.Selection.iteritems()):
                print_small_banner((tc + ':').ljust(15))
                for rp, ch in sorted(rps.iteritems()):
                    runs = self.get_runs(rp, tc)
                    print rp.ljust(5), '{0}-{1}'.format(str(runs[0]).zfill(3), str(runs[-1]).zfill(3))
        else:
            log_warning('Selection is empty!')

    def get_runs(self, rp, tc):
        return self.AllRunPlans[tc]['rate_scan'][rp]

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

    @staticmethod
    def scale_graph(mg, scale):
        try:
            mg.GetListOfGraphs()
            ymins, ymaxs = [], []
            for g in mg.GetListOfGraphs():
                g, ymin, ymax = DiaScans.scale_graph(g, scale)
                ymins.append(ymin)
                ymaxs.append(ymax)
            return mg, min(ymins), max(ymaxs)
        except AttributeError:
            g = mg
        n, x, y, ex, ey = get_graph_data(g)
        ey = [0] * n if g.GetName() in ['gLine', 'gFirst', 'gLast'] else ey
        y = map(lambda v: v * scale, y)
        ey = map(lambda v: v * scale, ey)
        for i in xrange(n):
            g.SetPoint(i, x[i], y[i])
            g.SetPointError(i, ex[i], ey[i])
        ymin = min(map(lambda val, err: val - err, y, ey))
        ymax = max(map(lambda val, err: val + err, y, ey))
        g.SetTitle('original' if scale == 1 else 'scaled by {0:3.2f}'.format(scale))
        return g, ymin, ymax

    def draw_rate_scans(self, show=False, pulser=False, include_zero=True, gridy=False, do_scale=True, draw_single_plots=False, scale_to=None, range_y=None):
        if not self.Selection:
            log_warning('Selection is empty!')
            return
        run_selections = self.load_run_selections()

        def scale_graphs(graphs, scale_to=None):
            print 'There are: {n} runs in the sample:'.format(n=len(graphs))
            i = 0
            biases = []
            for (sel, ch), (mg, mean_ph) in graphs.items():
                print ' ', i, sel.TESTCAMPAIGN, sel.SelectedBias, mean_ph
                biases.append(sel.SelectedBias)
                i += 1
            try:
                scale = float(scale_to)
                print 'scale', scale
                scale = scale, 1
            except:
                if len(graphs) == 1:
                    retval = 0
                if len(biases) == 2:
                    # print biases, biases[0], biases[1]
                    if biases[0] < 0 < biases[1]:
                        retval = 0
                    elif biases[0] > 0 > biases[1]:
                        retval = 1
                    else:
                        retval = -1
                while True:
                    try:
                        if retval < 0:
                            raise Exception()
                        i = int(retval)
                        scale = graphs.items()[i][1][1]
                        break
                    except:
                        retval = raw_input('Enter to which run the data should be scaled [0-{n}]: '.format(n=len(graphs) - 1))
                print 'data is scaled to run ', i, 'with ph', scale
                move_element(graphs, graphs.keys()[i], 0)
            ymins = []
            ymaxs = []
            for (sel, ch), (mg, mean_ph) in graphs.items():
                self.scale_factors[(sel, ch)] = scale[0] / mean_ph[0]
                mg, ymin, ymax = DiaScans.scale_graph(mg, scale[0] / mean_ph[0])
                mean_ph = mean_ph[:2] + (ymin, ymax)
                ymins.append(ymin)
                ymaxs.append(ymax)
            ymin = min(ymins)
            ymax = max(ymaxs)
            return graphs, ymin, ymax

        mg = TMultiGraph('mg_ph', '{dia} Rate Scans at {bias} V;Flux [kHz/cm^{{2}}]; pulse height [au]'.format(dia=self.DiamondName, bias=500))
        xndc = .47
        if do_scale:
            xndc = .2
        if not include_zero:
            nentries = len(run_selections)
            legend = self.make_legend(xndc, .2 + nentries * .05, nentries=nentries, scale=1., w=.7, x2=.94)
        else:
            legend = self.make_legend(xndc, .35, nentries=len(run_selections), scale=1., w=.7, x2=.94)
        legend.SetNColumns(4)
        # legend.SetX2NDC(.9)
        graphs = OrderedDict()
        legends = {}
        for sel, ch in run_selections.iteritems():
            # print 'get', ch, sel.Diamond, sel.TESTCAMPAIGN, sel.SelectedRunplan
            self.DiamondName = self.load_diamond(sel.Diamond)
            try:
                path = self.PickleDir + 'Ph_fit/PulseHeights_{tc}_{rp}_{dia}_{bin}.pickle'.format(tc=sel.TESTCAMPAIGN, rp=sel.SelectedRunplan, dia=self.DiamondName, bin=20000)
                if pulser:
                    path = self.PickleDir + 'Pulser/PulseHeights_{tc}_{rp}_{dia}.pickle'.format(tc=sel.TESTCAMPAIGN, rp=sel.SelectedRunplan, dia=self.DiamondName)
                f = open(path, 'r')
                mg_ph_ana = pickle.load(f)
                f.close()
            except Exception:
                Elementary(sel.TESTCAMPAIGN)
                ana = AnalysisCollection(sel, ch, self.verbose)
                mg_ph_ana = ana.draw_pulser_info(show=False, do_fit=False) if pulser else ana.draw_pulse_heights(show=False)
                ana.close_files()
            mean_ph = self.get_ph_below_flux(mg_ph_ana, flux=80)
            graphs[(sel, ch)] = (mg_ph_ana, mean_ph)
            l = self.make_legend(.65, .20, nentries=3)
            try:
                for g in mg_ph_ana.GetListOfGraphs():
                    if g.GetName() in ['gFullError', 'data']:
                        l.AddEntry(g, g.GetTitle(), 'pl')
                    elif g.GetName() == 'gStatError':
                        l.AddEntry(g, g.GetTitle(), 'l')
                    elif not g.GetName() == 'gLine':
                        l.AddEntry(g, g.GetTitle(), 'p')
            except:
                g = mg_ph_ana
                l.AddEntry(g, g.GetTitle(), 'pl')

            legends[(sel, ch)] = l
            self.ROOTObjects.append(l)

        if do_scale:
            graphs, ymin, ymax = scale_graphs(graphs, scale_to=scale_to)

        all_data = []
        for (sel, ch), (mg_ph_ana, scale) in graphs.iteritems():
            print_banner('CH: {0}\nScale: {1}'.format(ch, scale))
            resize_markers(mg_ph_ana, default_size=1, marker_sizes={'gFirst': 1.3, 'gLast': 1.3, 'gStatError': 0})
            if draw_single_plots and do_scale:
                y = ymin, ymax
                mg_y = y[0] * 1.3 - y[1] * .3
                delta_y = ymax - ymin
                y_range = increased_range(y, .2, .05)
                mg_y = y_range[0]
                self.format_histo(mg_ph_ana, y_range=y_range, y_off=1.75, x_off=1.3, x_tit='Flux [kHz/cm^{2}]', y_tit='pulse height [au]', do_marker=False)
                mg1 = mg_ph_ana.Clone()
                mg1.GetListOfGraphs()[0].SetLineColor(602)
                self.ROOTObjects.append(mg1)
                scale = self.scale_factors[(sel, ch)]
                hname = 'Combined{pre}{scale}PulseHeights_{dia}_{bias}_{tc}'.format(pre='Pulser' if pulser else '', dia=self.load_diamond(sel.Diamond), tc=sel.TESTCAMPAIGN,
                                                                                    scale='Scaled' if do_scale else'', bias=self.make_bias_string(sel.SelectedBias))
                hname = hname.replace('__', '_')
                hname = hname.replace('-', '_')

                runinfo = sel.get_runinfo(0 if ch == 1 else 3)
                if scale != 1:
                    legScale = self.make_legend(.70, .97, 1, scale=1.2)
                    legScale.AddEntry(None, 'Scaled by {:.2f}'.format(scale), '')
                    self.ROOTObjects.append(legScale)
                    draw_objects = [(legScale, '')]
                else:
                    draw_objects = None
                self.save_combined_pulse_heights(mg_ph_ana, mg1, legends[(sel, ch)], mg_y, show=True, name=hname, pulser_leg=None,
                                                 x_range=None, y_range=y_range, rel_y_range=None, run_info=runinfo, draw_objects=draw_objects)
                sel.run.scale_runinfo_legend(txt_size=.075, w=.35, h=0.1 / .268)
                self.save_plots(hname)
            bias = sel.SelectedBias
            color = self.get_next_color(bias, pulser)
            set_graph_color(mg_ph_ana, color, marker_size=0)
            resize_markers(mg_ph_ana, default_size=0, marker_sizes={'gFirst': 1.3, 'gLast': 1.3})
            self.ROOTObjects.append(mg_ph_ana)
            mg.Add(mg_ph_ana)
            tc = sel.print_testcampaign(pr=False)
            title = ''
            try:
                if pulser:
                    g = mg_ph_ana.GetListOfGraphs().FindObject('data')
                else:
                    g = mg_ph_ana.GetListOfGraphs().FindObject('gFullError')
                if g == None:
                    print [g.GetName() for g in mg_ph_ana.GetListOfGraphs()]
            except Exception as err:
                print log_warning(err)
                g = mg_ph_ana

            all_data.append(get_graph_data(g))
            # print all_data[-1]
            if g == None:
                print g, mg_ph_ana
                raise Exception
            if do_scale:
                title += g.GetTitle() + ', '
            else:
                title = ' '
            title3 = '{dia}'.format(dia=self.load_diamond(sel.Diamond))
            title4 = 'in {tc}, Bias'.format(tc=tc, )
            title2 = ' {bias:+5.0f} V'.format(tc=tc, bias=bias)
            title2 = title2.replace('-', '#minus ')
            title2 = title2.replace('+', '#plus ')
            legend.AddEntry(g, title, 'l')
            legend.AddEntry(None, title3, '')
            legend.AddEntry(None, title4, '')
            legend.AddEntry(None, title2, '')
        self.all_data = all_data
        tcs = '_'.join(self.Selection.keys())
        if include_zero:
            extension = 'Zero'
        else:
            extension = ''
        if do_scale:
            extension += 'Scaled'
        prefix = 'Pulser' if pulser else ''
        pname = 'Combined_{prefix}PH_{key}_RateScans{ext}_{tcs}'.format(tcs=tcs, ext=extension, prefix=prefix, key=self.Name).replace('-', '_')
        if self.DiamondName.replace('-', '_') not in pname: pname += '_' + self.DiamondName.replace('-', '_')
        pname = pname.replace('__', '_')
        self.format_histo(mg, draw_first=True, y_off=1.4, x_tit='Flux [kHz/cm^{2}]', y_tit='{0} Pulse Height [au]'.format('Pulser' if pulser else 'Signal'))
        mg.GetXaxis().SetLimits(3.5, 20e3)
        ymax = mg.GetYaxis().GetXmax()
        ymin = mg.GetYaxis().GetXmin()
        if range_y is None:
            if include_zero:
                range_y = [0, ymax * 1.2]
            else:
                range_y = [(ymin - .42 * (ymax - ymin)), (ymax + .00 * (ymax - ymin))]
        print 'Y Range:', range_y

        mg.SetMinimum(range_y[0])
        mg.SetMaximum(range_y[1])

        self.MainConfigParser.set('SAVE', 'short_name', 'True')
        legend.SetX2NDC(.94)

        draw_opt = 'A' if not pulser else 'APL'
        self.save_histo(mg, pname, lm=.14, draw_opt=draw_opt, l=legend, logx=True, gridy=gridy)

        return mg

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

    def create_up_down_scans(self, mg, pulser=False):
        scans = self.get_up_down_scans(mg)
        self.draw_hysteresis_graph(scans)
        graphs = {}
        mgg = TMultiGraph('mg_up_down', '{dia} Rate Scans at {bias} V;Flux [kHz/cm^{{2}}]; pulse height [au]'.format(dia=self.DiamondName, bias=500))
        xndc = .4

        legend = self.make_legend(xndc, .35, nentries=3, scale=1., w=.2)
        dia_names = self.get_diamond_names()
        biases = self.get_bias_voltages()
        bias_str = get_bias_root_string(biases)
        print 'diamond names', dia_names
        print 'biases', biases, bias_str
        legend.SetHeader('{names}, scaled, V: {bias}'.format(names=dia_names[0] if len(dia_names) == 1 else dia_names, bias=bias_str))
        for key, d in scans.iteritems():
            n = len(d['x'])
            gr = TGraphErrors(n)
            for i in range(n):
                gr.SetPoint(i, d['x'][i], d['y'][i])
                gr.SetPointError(i, d['ex'][i], d['ey'][i])
            gr.SetTitle('Flux ' + key)
            if key == 'up':
                gr.SetMarkerStyle(22)
                color = kRed
            elif key == 'down':
                gr.SetMarkerStyle(23)
                color = kBlue
            else:
                print 'cannot find', key
            gr.SetLineColor(color)
            gr.SetMarkerColor(color)
            graphs[key] = gr
            self.format_histo(gr, y_off=1.4, x_tit='Flux [kHz/cm^{2}]', y_tit='{0} Pulse Height [au]'.format('Pulser' if pulser else 'Signal'), color=color)
            mgg.Add(gr, 'p')
            legend.AddEntry(gr, gr.GetTitle(), 'lp')
        self.ROOTObjects.append(graphs)
        self.MainConfigParser.set('SAVE', 'short_name', 'True')
        pname = 'Combined_UpDown_PH_RateScans_{key}'.format(key=self.Name).replace('-', '_')
        dia = '_'.join(dia_names).replace('-', '_')
        if dia not in pname:
            pname += '_' + dia
        pname = pname.replace('__', '_')
        draw_opt = 'A'
        if pulser:
            draw_opt += 'PL'
        self.format_histo(mgg, y_off=1.4, x_off=1.3, x_tit='Flux [kHz/cm^{2}]', y_tit='{0} Pulse Height [au]'.format('Pulser' if pulser else 'Signal'), draw_first=True)
        self.draw_histo(mgg, pname, lm=.14, draw_opt=draw_opt, l=legend, logx=True)
        ymin = mgg.GetYaxis().GetXmin()
        ymax = mgg.GetYaxis().GetXmax()
        mgg.SetMinimum(ymin - .3 * (ymax - ymin))
        legend.SetX2NDC(.94)
        self.ROOTObjects.append(
            self.save_plots(pname, self.save_dir)
        )
        return mgg, scans

    @staticmethod
    def get_up_down_scans(mg, keys=None):
        keys = ['gFullError', 'data'] if keys is None else keys
        scans = {}
        for i in ['up', 'down']:
            scans[i] = {
                'x': [],
                'y': [],
                'ex': [],
                'ey': []
            }
        for g in mg.GetListOfGraphs():
            if g.GetName() not in keys:
                continue
            n, x, y, ex, ey = get_graph_data(g)
            if n < 2:
                continue
            if x[0] < x[1]:
                key = 'up'
            else:
                key = 'down'
            i = 0
            while True:
                scans[key]['x'].append(x[i])
                scans[key]['y'].append(y[i])
                scans[key]['ex'].append(ex[i])
                scans[key]['ey'].append(ey[i])
                if i == n - 1:
                    break
                if x[i] < x[i + 1] and key is not 'up':
                    key = 'up'
                    # print 'changed to up', x[i], x[i + 1]
                elif x[i] > x[i + 1] and key is not 'down':
                    key = 'down'
                    # print 'changed to down', x[i], x[i + 1]
                else:
                    i += 1
        return scans

    def load_run_selections(self, redo=False):
        if self.RunSelections is not None and not redo:
            return self.RunSelections
        run_selections = OrderedDict()
        for tc, rps in sorted(self.Selection.iteritems()):
            for rp, chs in sorted(rps.iteritems()):
                if type(chs) is not list:
                    chs = [chs]
                for ch in chs:
                    sel = RunSelection(tc)
                    sel.select_runs_from_runplan(rp, ch=ch)
                    self.log_info('Loaded runplan {rp} of testcampaign {tc} and ch {ch} ({dia})'.format(rp=rp.rjust(4), tc=datetime.strptime(tc, '%Y%m').strftime('%b %Y'), ch=ch, dia=sel.Diamond))
                    run_selections[sel] = ch
        self.RunSelections = run_selections
        return run_selections

    def create_combined_plots(self, flux_up_down=False, scaled=True):
        mg = self.draw_rate_scans(include_zero=False, do_scale=scaled, pulser=False, draw_single_plots=True)
        if flux_up_down and scaled:
            mg, scan = self.create_up_down_scans(mg)
            self.draw_hysteresis_graph(scan)
        self.draw_rate_scans(include_zero=False, do_scale=scaled, pulser=True, draw_single_plots=True)

    def plots_for_felix(self):
        pass
        for key in ['poly-B2_irradiated', 'poly-B2_unirradiated', 'poly-D']:
            self.set_selection(key)
            self.draw_rate_scans(include_zero=True, do_scale=False)

    def show_selections(self):
        print 'The following selections exists:'
        for key, sel in self.Selections.items():
            print '* "{key}": {sel}'.format(key=key, sel=sel)

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

    @staticmethod
    def calc_flux(info, tc):
        if 'for1' not in info or info['for1'] == 0:
            if 'measuredflux' in info:
                # return str('{0:5.0f}'.format(info['measuredflux'] * 2.48))
                return
        path = '/data/psi_{0}_{1}/masks/{mask}'.format(tc[:4], tc[-2:], mask=info['maskfile'])
        if file_exists(path):
            f = open(path, 'r')
        else:
            return
        data = []
        for line in f:
            if len(line) > 3:
                line = line.split()
                data.append([int(line[2])] + [int(line[3])])
        f.close()
        pixel_size = 0.01 * 0.015
        area = [(data[1][0] - data[0][0]) * (data[1][1] - data[0][1]) * pixel_size, (data[3][0] - data[2][0]) * (data[3][1] - data[2][1]) * pixel_size]
        # print area
        flux = [info['for{0}'.format(i + 1)] / area[i] / 1000. for i in xrange(2)]
        return mean(flux)


if __name__ == '__main__':
    main_parser = ArgumentParser()
    main_parser.add_argument('dia', nargs='?', default='S129')
    main_parser.add_argument('-tcs', nargs='?', default=None)
    args = main_parser.parse_args()
    print args
    print_banner('STARTING DIAMOND RATE SCAN COLLECTION OF DIAMOND {0}'.format(args.dia))

    z = DiaScans(args.dia, args.tcs, verbose=True)
    if False:
        for key_, s in z.Selections.items():
            print_banner('STARTING SELECTION {0}'.format(key_))
            z.set_selection(key_)
            z.create_combined_plots()
            if key_ in ['poly-D', 'poly-B2_irradiated']:
                z.create_combined_plots(flux_up_down=True)
    z.set_selection('poly-D')
    # if raw_input('press y for creating all plots').lower() == 'y':
    #     for key in ['poly-B2_irradiated',  # ,'poly-B2_unirradiated',
    #                 'poly-D'
    #                 ]:
    #         z.set_selection(key)
    #         mg = z.draw_rate_scans(include_zero=False, do_scale=True, pulser=False, draw_single_plots=True)
    #         mg, scan = z.create_up_down_scans(mg)
    #         # z.create_combined_plots(flux_up_down=True)
    #         # z.create_combined_plots(scaled=False)
    # z.set_selection('ext_neg_pulser_neg_pol')
    # z.create_combined_plots()
    # for key in [#'S129','S129_201510',
    #        'S129_n500','poly-B2']:
    #    z.set_selection(key)
    #    z.create_combined_plots()
    #    z.create_combined_plots(scaled=False)
    #    z.draw_rate_scans(include_zero=True,do_scale=False)
    # z.set_selection(u'ext_neg_pulser_neg_pol')
    # z.draw_rate_scans(pulser=True,do_scale=True,include_zero=False,reset_colors=False)
    # z.set_selection('ext_pul_pos_pol')
    #         # z.used_colors=[kGreen, kCyan, kViolet]
    #         # scale = raw_input('scale to:')
    #         # z.draw_rate_scans(pulser=True,do_scale=True,include_zero=False,reset_colors=False,scale_to=scale)
