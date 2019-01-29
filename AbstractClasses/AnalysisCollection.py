#! /usr/bin/env python
# ==============================================
# IMPORTS
# ==============================================
from argparse import ArgumentParser
from numpy import log, concatenate, zeros, sign

from ROOT import gROOT, TCanvas, TLegend, TExec, gStyle, TMultiGraph, THStack, TF1, TH1F, TH2F, TH2I, TProfile2D, TProfile, TCut

from CurrentInfo import Currents
from Elementary import Elementary
from Extrema import Extrema2D
from PadAnalysis import PadAnalysis
from RunSelection import RunSelection
from TelescopeAnalysis import Analysis
from VoltageScan import VoltageScan
from PulserCollection import PulserCollection
from InfoLegend import InfoLegend
from Run import Run
from Utils import *


# ==============================================
# MAIN CLASS
# ==============================================
class AnalysisCollection(Elementary):
    """
    An object of this class contains several analysis of runs.
    It gives the ability to compare the data from different runs.
    """
    current_run_number = -1

    def __init__(self, run_selection, threads=None):
        self.Runs = run_selection.get_selected_runs()
        Elementary.__init__(self, verbose=run_selection.verbose)
        # dict where all analysis objects are saved
        self.collection = OrderedDict()
        self.selection = run_selection

        self.Threads = threads
        self.LoadTrees = threads.values()[0].Tuple
        self.Type = run_selection.SelectedType
        # self.diamonds = self.load_diamonds(diamonds, list_of_runs)
        self.min_max_rate_runs = self.get_high_low_rate_runs()
        if self.LoadTrees:
            self.generate_slope_pickle()
            self.generate_threshold_pickle()
        self.add_analyses()
        self.FirstAnalysis = self.get_first_analysis()
        self.LastAnalysis = self.get_last_analysis()
        if self.LoadTrees:
            self.Plots = self.FirstAnalysis.Plots
        self.channel = self.load_channel()

        self.signalValues = None

        # important information
        self.NRuns = len(self.collection)
        self.DiamondName = self.FirstAnalysis.DiamondName
        self.DiamondNumber = self.FirstAnalysis.DiamondNumber
        self.Bias = self.FirstAnalysis.Bias

        # root stuff
        self.RunPlan = run_selection.SelectedRunplan
        self.save_dir = '{dia}/runplan{plan}'.format(tc=self.TESTCAMPAIGN[2:], plan=self.RunPlan, dia=self.DiamondName)
        self.RootObjects = []

        # sub classes
        self.StartTime = time_stamp(self.FirstAnalysis.Run.LogStart)
        self.VoltageScan = VoltageScan(self)
        self.InfoLegend = InfoLegend(self)
        self.Currents = Currents(self)
        self.Pulser = PulserCollection(self)

    def __del__(self):
        print '\ndeleting AnalysisCollection...'
        for nr, ana in self.collection.iteritems():
            ana.__del__()

    def close_files(self):
        for ana in self.collection.itervalues():
            ana.Run.tree.Delete()
            ana.Run.RootFile.Close()

    # ============================================
    # region INIT

    def delete_trees(self):
        for ana in self.collection.itervalues():
            ana.tree.Delete()

    def add_analyses(self):
        """ Creates and adds Analysis objects with run numbers in runs. """
        for run in self.Runs:
            thread = self.Threads[run]
            run_class = Run(run, tree=thread.Tuple, t_vec=thread.Time, verbose=self.verbose)
            analysis = PadAnalysis(run_class, self.selection.SelectedDiamondNr, self.min_max_rate_runs)
            self.collection[analysis.Run.RunNumber] = analysis
            self.current_run_number = analysis.Run.RunNumber

    def load_channel(self):
        binary = self.FirstAnalysis.run_config_parser.getint('ROOTFILE_GENERATION', 'active_regions')
        dia_nr = self.selection.SelectedDiamondNr
        return [i for i in xrange(16) if has_bit(binary, i)][dia_nr - 1]

    def get_high_low_rate_runs(self):
        runs = [Run(run_number=run, tree=False) for run in self.Runs]
        fluxes = OrderedDict()
        self.log_info('RUN FLUX [kHz/cm2]')
        for run in runs:
            fluxes[run.Flux] = run.RunNumber
            self.log_info('{run:3d} {flux:14.2f}'.format(run=run.RunNumber, flux=run.Flux.n))
        print '\n'
        return {'min': fluxes[min(fluxes)], 'max': fluxes[max(fluxes)]}

    def generate_slope_pickle(self):
        picklepath = self.make_pickle_path('TrackAngle', 'x', run=self.min_max_rate_runs['min'])
        if file_exists(picklepath):
            return
        Analysis(Run(self.min_max_rate_runs['min'], verbose=self.verbose))

    def generate_threshold_pickle(self):
        picklepath = self.make_pickle_path('Cuts', 'SignalThreshold', run=self.min_max_rate_runs['max'], ch=self.selection.SelectedDiamondNr)
        if file_exists(picklepath) or self.Type != 'pad':
            return
        PadAnalysis(self.min_max_rate_runs['max'], dia=self.selection.SelectedDiamondNr)

    def get_hv_name(self):
        return self.Currents.Name

    # endregion

    def create_all_single_run_plots(self):
        old_verbose = self.FirstAnalysis.verbose
        self.set_verbose(False)
        for key, ana in self.collection.iteritems():
            print 'Create Plots for Run ', key
            ana.compare_consecutive_cuts(scale=False)
            ana.compare_consecutive_cuts(scale=True)
            ana.show_cut_contributions(show=False)
            ana.draw_bucket_pedestal(show=False)
        self.set_verbose(old_verbose)

    # ============================================
    # region SIGNAL/PEDESTAL
    def print_all_off_results(self):
        string = 'Signal\tPedest.\tPulser\n'
        for ana in self.collection.itervalues():
            string += '{0}\n'.format(ana.print_off_results(prnt=False))
        print string

    def draw_all(self):

        old_verbose = self.FirstAnalysis.verbose
        self.set_verbose(False)
        if self.Type == 'voltage scan':
            self.VoltageScan.draw_all()
        else:
            self.draw_pulse_heights(show=False)
            self.draw_scaled_pulse_heights(show=False)
            self.draw_scaled_pulse_heights(show=False, vs_time=True)
            self.Pulser.draw_pulse_heights(do_fit=False, show=False)
            self.draw_pedestals(show=False)
            self.draw_noise(show=False)
            self.draw_pulser_pedestals(show=False)
        self.draw_ph_with_currents(show=False)
        self.draw_signal_distributions(show=False)
        self.save_signal_maps()
        self.save_signal_maps(hitmap=True)
        self.set_verbose(old_verbose)
        print
        self.print_all_off_results()

    def draw_little_all(self, redo=False):
        t0 = self.log_info('Generate all plots ... ')
        old_verbose = self.FirstAnalysis.verbose
        self.set_verbose(False)
        if self.Type == 'voltage scan':
            self.VoltageScan.draw_all(redo)
        else:
            self.draw_pulse_heights(show=False, redo=redo)
            self.draw_scaled_pulse_heights(show=False)
            self.draw_scaled_pulse_heights(show=False, vs_time=True)
            self.Pulser.draw_pulse_heights(show=False, redo=redo)
            self.draw_pedestals(show=False, redo=redo)
            self.draw_pedestals(show=False, sigma=True, redo=redo)
            self.draw_pulser_pedestals(show=False, redo=redo)
        self.draw_ph_with_currents(show=False)
        self.draw_signal_distributions(show=False, redo=redo)
        self.save_signal_maps(redo=redo)
        self.save_signal_maps(hitmap=True, redo=redo)
        self.draw_run_currents()
        self.draw_chi2s()
        self.set_verbose(old_verbose)
        print_elapsed_time(t0)
    
    def scale_current_gr(self, gr):
        vals = [gr.GetY()[i] for i in xrange(gr.GetN())]
        mi, ma, mn = min(vals), max(vals), mean(vals)
        y = [mi, mn + 3 * (mn - mi) if ma - mn > 50 else ma]
        self.format_histo(gr, name='cur', color=899)
        return increased_range(y, .3, .2)

    def draw_ph_with_currents(self, show=True, scale=1):
        ph = self.get_pulse_height_graph(vs_time=True, first_last=False, binning=10000, legend=False)
        self.Currents.set_graphs()
        cur = self.Currents.CurrentGraph.Clone()
        cur_range = self.scale_current_gr(cur)
        pul = self.Pulser.get_pulse_height_graph(vs_time=True, legend=False, show_flux=False)
        pul_fac = pul.GetListOfGraphs()[0].GetY()[0] if scale else ph.GetListOfGraphs()[0].GetY()[0] / pul.GetListOfGraphs()[0].GetY()[0]
        ph_fac = ph.GetListOfGraphs()[0].GetY()[0] if scale else sign(self.Bias)
        scale_multigraph(pul, ph.GetListOfGraphs()[0].GetY()[0] if scale is None else scale)
        scale_multigraph(ph, scale)

        # legends
        entries = [2, 3, 1]
        lm = .09
        positions = [[lm + .02, .95], [lm + .02, .95], [lm + .02, .95]]
        legends = [self.make_legend(*positions[i], nentries=entries[i], scale=2.5 * (27 / 36. if i == 2 else 1), x2=lm + .22) for i in xrange(3)]
        dummy_gr = self.make_tgrapherrors('g', 'g', width=2)
        legends[0].AddEntry(pul.GetListOfGraphs()[0], 'scaled pulser data', 'p')
        legends[0].AddEntry(0, 'factor {sgn}{fac:4.0f}'.format(fac=pul_fac, sgn='-' if self.FirstAnalysis.PulserPolarity < 0 else ''), '')
        legends[1].AddEntry(ph.GetListOfGraphs()[1], '{sgn:2.0f} * (signal - ped.) data'.format(sgn=ph_fac), 'p')
        legends[1].AddEntry(0, 'flux in kHz/cm^{2}', 'p')
        legends[1].AddEntry(dummy_gr, 'duration', 'l')
        legends[2].AddEntry(cur, 'current', 'l')

        # Drawing
        lab_size = .12
        for l in ph.GetListOfGraphs()[0].GetListOfFunctions():
            l.SetTextSize(.09)
            l.SetY(l.GetY() - .01)
        self.format_histo(ph, draw_first=True, lab_size=lab_size, tit_size=lab_size, y_off=.33)
        self.format_histo(cur, x_tit='Time [hh:mm]', lab_size=lab_size * 27 / 36., tit_size=lab_size * 27 / 36., y_off=.33 * 36 / 27.)
        self.format_histo(pul, color=859, draw_first=True, lab_size=lab_size, tit_size=lab_size, y_off=.33)
        self.format_histo(pul.GetListOfGraphs()[0], color=859)
        draw_opts = ['a', 'a', '']
        y_tits = ['Pulser ', 'Signal', 'Current [nA] ']
        y_ranges = [increased_range(scale_margins(pul, ph), .3), increased_range(scale_margins(ph, pul), .3), cur_range]
        y_ranges = [[.94, 1.06], [.94, 1.06], cur_range] if scale else y_ranges
        divs = [504, 504, None]
        ypos = [[1, .73], [.73, .46], [.46, .1]]
        margins = [[lm, .05, 0, 0], [lm, .05, 0, 0], [lm, .05, 9 / 36., 0]]
        x_range = increased_range([cur.GetX()[0], cur.GetX()[cur.GetN() - 1]], .1, .1)
        close_last_canvas()

        c = self.make_canvas('c_phc', 'c', 1.5, 1, transp=True, show=show)

        for i, gr in enumerate([pul, ph, cur]):
            self.format_histo(gr, title='', y_range=y_ranges[i], y_tit=y_tits[i], center_y=True, ndivy=divs[i], t_ax_off=self.Currents.Time[0])
            self.draw_tpad('p{}'.format(i), 'p{}'.format(i), pos=[0, ypos[i][0], 1, ypos[i][1]], gridx=True, gridy=True, margins=margins[i], transparent=True)
            gr.GetXaxis().SetLimits(*x_range)
            gr.Draw(draw_opts[i])
            legends[i].Draw()
            c.cd()

        self.save_plots('PhPulserCurrent', all_pads=False)
        self.RootObjects.append([ph, cur, pul, c, legends])

    def draw_slope_vs_voltage(self, show=True, gr=False):
        h = TH1F('hSV', 'PH Slope Distribution', 10, -1, 1) if not gr else self.make_tgrapherrors('gSV', 'PH Slope vs. Voltage')

        self.start_pbar(self.NRuns)
        for i, (key, ana) in enumerate(self.collection.iteritems()):
            ana.draw_pulse_height(corr=True, show=False)
            fit = ana.PulseHeight.Fit('pol1', 'qs')
            x = ana.Run.RunInfo['dia{nr}hv'.format(nr=self.FirstAnalysis.Run.channels.index(self.channel) + 1)]
            s, e = fit.Parameter(1), fit.ParError(1)
            if gr:
                h.SetPoint(i, x, s)
                h.SetPointError(i, 0, e)
            else:
                set_statbox(entries=4, only_fit=True)
                h.Fill(s)
                h.Fit('gaus', 'q')
            self.ProgressBar.update(i + 1)
        self.format_histo(h, x_tit='Slope [mV/min]', y_tit='Number of Entries', y_off=1.3)
        self.draw_histo(h, show=show, draw_opt='alp' if gr else '')

    def get_pulse_heights(self, binning=None, redo=False):
        pickle_path = self.make_pickle_path('Ph_fit', 'PhVals', self.RunPlan, ch=self.DiamondName, suf=self.FirstAnalysis.BinSize if binning is None else binning)

        def func():
            self.log_info('Getting pulse heights ... ')
            phs = OrderedDict()
            self.start_pbar(self.NRuns)
            for i, (key, ana) in enumerate(self.collection.iteritems()):
                ph = ana.draw_pulse_height(bin_size=binning, corr=True, show=False, redo=redo, prnt=False)[1]
                phs[key] = {'flux': ana.get_flux(), 'ph': make_ufloat(ph, par=0), 'time': ana.Run.get_time()}
                self.ProgressBar.update(i + 1)
            self.ProgressBar.finish()
            return phs

        pulse_heights = func() if redo else None
        return do_pickle(pickle_path, func, pulse_heights)

    def draw_log_flux(self, evts_per_bin=10000, rel_t=True, show=True):
        end_times = [t1[1] for t1 in self.get_start_end_times()]
        h = TH1F('hr1', 'Flux of Run Plan {r}'.format(r=self.RunPlan), *self.get_t_binning(evts_per_bin))
        ibin = 1
        for i, end_time in enumerate(end_times):
            ana = self.collection.values()[i]
            while h.GetBinCenter(ibin) <= end_time:
                h.SetBinContent(ibin, ana.Run.Flux)
                h.SetBinError(ibin, ana.Run.Flux * .1)
                ibin += 1
        self.format_histo(h, y_tit='Flux [kHz/cm^{2}]', x_tit='Time [hh:mm]', y_off=.8, fill_color=self.FillColor, stats=0)
        set_time_axis(h, off=self.FirstAnalysis.Run.StartTime if rel_t else 0)
        self.save_histo(h, 'FluxTime', show=show, draw_opt='hist', x_fac=1.5, y_fac=.75, lm=.065, gridy=True, logy=True)
        return h

    def draw_full_pulse_height(self, evts_per_bin=10000, show=True, rel_t=True, redo=False, with_flux=True):

        pickle_path = self.make_pickle_path('PulseHeight', 'FullPH', self.RunPlan, ch=self.DiamondNumber, suf=evts_per_bin)

        def f():
            histos = [ana.draw_pulse_height(evts_per_bin, corr=True, redo=redo, show=False)[0] for ana in self.collection.itervalues()]
            h1 = TH1F('hfph', 'Pulse Height for Run Plan {n}'.format(n=self.RunPlan), *self.get_t_binning(evts_per_bin))
            i_bin = 0  # the first bin is empty
            for hist in histos:
                for i in xrange(1, hist.GetNbinsX() + 1):
                    h1.SetBinContent(i_bin + 1, hist.GetBinContent(i))
                    h1.SetBinError(i_bin + 1, hist.GetBinError(i))
                    i_bin += 1
                i_bin += 1  # there is an empty bin after each run
            return h1

        histo = do_pickle(pickle_path, f, redo=redo)
        self.format_histo(histo, x_tit='Time [hh:mm]', y_tit='Mean Pulse Height [au]', y_off=.8, fill_color=self.FillColor, stats=0, y_range=[0, histo.GetMaximum() * 1.05])
        set_time_axis(histo, off=self.FirstAnalysis.Run.StartTime if rel_t else 0)
        c = self.draw_histo(histo, show=show, draw_opt='hist', x=1.5, y=.75, lm=.065, gridy=True, rm=.1 if with_flux else None)
        if with_flux:
            h = self.draw_fluxes(rel_time=rel_t, show=False)
            c.cd()
            self.draw_tpad('pr', margins=[.065, .1, .15, .1], transparent=True, logy=True)
            x_range = [histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax()]
            self.format_histo(h, title=' ', fill_color=2, fill_style=3002, lw=1, y_range=[1, h.GetMaximum() * 1.2], stats=0, y_off=1.05, x_range=x_range)
            h.Draw('histy+')
        self.save_plots('FullPulseHeight')

    def get_pulse_height_graph(self, binning=None, vs_time=False, first_last=True, redo=False, legend=True):

        marker_size = 1
        gStyle.SetEndErrorSize(4)
        ph = self.get_pulse_heights(binning, redo)
        y_values = [dic['ph'] for dic in ph.itervalues()]
        x_values = [dic['time' if vs_time else 'flux'] for dic in ph.itervalues()]
        g = self.make_tgrapherrors('g', 'stat. error', self.get_color(), marker_size=marker_size, x=x_values, y=y_values)
        rel_sys_error = self.get_repr_error(105, show=False)
        y_values = [make_ufloat((v.n, v.s + rel_sys_error * v.n)) for v in y_values]
        g_errors = self.make_tgrapherrors('gerr', 'full error', marker=0, color=602, marker_size=0, x=x_values, y=y_values)
        g_first = self.make_tgrapherrors('g1', 'first run', marker=22, color=2, marker_size=marker_size, x=[x_values[0].n], y=[y_values[0].n])
        g_last = self.make_tgrapherrors('g2', 'last run', marker=23, color=2, marker_size=marker_size, x=[x_values[-1].n], y=[y_values[-1].n])
        graphs = [g, g_errors]
        graphs += [g_first, g_last] if first_last else []
        l = self.make_legend(.75, .37, nentries=len(graphs))
        mg = TMultiGraph('mg_ph', 'Pulse Height vs {mod} - {dia}'.format(mod='Time' if vs_time else 'Flux', dia=self.DiamondName))
        for gr in graphs:
            l.AddEntry(gr, gr.GetTitle(), 'l' if gr.GetName() == 'gerr' else 'p')
            mg.Add(gr, 'pl')
        if legend:
            mg.GetListOfFunctions().Add(l)
        self.reset_colors()
        if vs_time:
            g = mg.GetListOfGraphs()[0]
            for i, (ana, x) in enumerate(zip(self.collection.itervalues(), x_values)):
                y, ey = g.GetY()[i], g.GetErrorY(i)
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(self.draw_tlatex(x.n, y + ey * 1.2, '{:1.0f}'.format(ana.get_flux().n), color=1, align=21, size=.02, show=0))
        return mg

    def draw_scaled_pulse_heights(self, scale=1, binning=None, vs_time=False, show=True, y_range=None, redo=False):

        mg = self.get_pulse_height_graph(binning, vs_time, first_last=not vs_time, redo=redo)
        scale_multigraph(mg, scale)
        xtit = 'Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]'
        y_range1 = [.95, 1.05] if y_range is None else y_range
        self.format_histo(mg, x_tit=xtit, y_tit='Scaled Pulse Height', y_off=1.75, x_off=1.3, draw_first=True, t_ax_off=0 if vs_time else None, y_range=y_range1, ndivx=503, center_y=1)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        move_legend(mg.GetListOfFunctions()[0], .75, .20)
        self.draw_histo(mg, '', show, lm=.14, draw_opt='a', logx=not vs_time, grid=vs_time, gridy=True, bm=.18)
        self.draw_irradiation(make_irr_string(self.selection.get_irradiation()))
        self.save_plots('ScaledPulseHeights{}'.format(xtit[:4]))
        return mg.GetListOfGraphs()[0]

    def draw_pulse_heights(self, binning=None, vs_time=False, show=True, show_first_last=True, save_comb=True, y_range=None, redo=False):

        mode = 'Time' if vs_time else 'Flux'
        mg = self.get_pulse_height_graph(binning, vs_time, show_first_last, redo)

        # small range
        self.format_histo(mg, color=None, x_tit='Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]', y_tit='Signal Pulse Height [au]', y_off=1.75, x_off=1.3, draw_first=True,
                          t_ax_off=0 if vs_time else None)
        ymin, ymax = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        yrange = increased_range([ymin, ymax], .5, .15) if y_range is None else y_range
        mg.GetYaxis().SetRangeUser(*yrange)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        self.save_histo(mg, 'PulseHeight{mod}'.format(mod=mode), show=False, lm=.14, draw_opt='A', logx=not vs_time, grid=vs_time)

        # no zero suppression
        mg1 = mg.Clone()
        mg1.GetListOfFunctions().Clear()
        mg1.SetName('mg1_ph')
        mg1.GetListOfGraphs()[0].SetLineColor(self.colors[0])
        mg1.GetYaxis().SetRangeUser(0, ymax * 1.1)
        self.save_histo(mg1, 'PulseHeightZero{mod}'.format(mod=mode), not save_comb, self.save_dir, lm=.14, draw_opt='Al', logx=not vs_time)
        self.reset_colors()

        if save_comb:
            y_min = increased_range([ymin, ymax], .5)[0] if y_range is None else y_range[0]
            # TODO fix vs time and comb plot
            self.save_combined_pulse_heights(mg, mg1, y_min, show=show, pulser_leg=self.__draw_signal_legend)

        return mg

    def get_pedestals(self, cut=None, pulser=False, redo=False):

        cut_name = '' if cut is None else TCut(cut).GetName()
        pickle_path = self.make_pickle_path('Pedestal', 'Values', self.RunPlan, self.DiamondNumber, suf='{}{}'.format('Pulser' if pulser else '', cut_name))

        def f():
            self.log_info('Getting {}pedestals ... '.format('pulser ' if pulser else ''))
            pedestals = OrderedDict()
            self.start_pbar(self.NRuns)
            for i, (key, ana) in enumerate(self.collection.iteritems()):
                ped = ana.Pulser.draw_pedestal(show=False, redo=redo, prnt=False) if pulser else ana.Pedestal.draw_disto_fit(cut=cut, show=False, redo=redo, prnt=False)
                pedestals[key] = {'flux': ana.get_flux(), 'ph': make_ufloat(ped, par=1), 'sigma': make_ufloat(ped, par=2), 'time': ana.Run.get_time()}
                self.ProgressBar.update(i + 1)
            self.ProgressBar.finish()
            return pedestals

        return do_pickle(pickle_path, f, redo=redo)

    def draw_pedestals(self, vs_time=False, sigma=False, cut=None, pulser=False, redo=False, show=True):

        mode = 'Time' if vs_time else 'Flux'
        pedestals = self.get_pedestals(cut, pulser, redo)
        y_values = [dic['sigma' if sigma else 'ph'] for dic in pedestals.itervalues()]
        x_values = [dic['time' if vs_time else 'flux'] for dic in pedestals.itervalues()]
        g = self.make_tgrapherrors('gps', '{}Pedestals'.format('Pulser ' if pulser else ''), x=x_values, y=y_values)
        self.format_histo(g, color=810, x_tit=self.make_x_tit(vs_time), y_tit='Mean Pedestal [au]', y_off=1.45, t_ax_off=0 if vs_time else None)
        cut_name = '' if cut is None else TCut(cut).GetName()
        save_name = '{p}Pedestal{s}{mod}{cut}'.format(mod=mode, cut=cut_name, s='Sigma' if sigma else 'Mean', p='Pulser' if pulser else '')
        self.save_histo(g, save_name=save_name, show=show, logx=False if vs_time else True, lm=.12, draw_opt='ap')
        return g

    def draw_noise(self, flux=True, show=True):
        return self.draw_pedestals(vs_time=flux, show=show, sigma=True)

    def draw_pulser_pedestals(self, show=True, redo=False):
        self.draw_pedestals(pulser=True, show=show, redo=redo)

    def draw_signal_distributions(self, show=True, off=3, redo=False):

        stack = THStack('hsd', 'Pulse Height Distributions')
        legend = self.make_legend(.67, .88, nentries=self.get_number_of_analyses())
        log_message('Generating signal distributions!')
        histos = []
        self.start_pbar(self.NRuns)
        for i, ana in enumerate(self.collection.itervalues(), 1):
            h = ana.draw_signal_distribution(show=False, redo=redo, prnt=False)
            self.format_histo(h, fill_color=0, fill_style=4000)
            histos.append(h)
            self.ProgressBar.update(i)
        self.ProgressBar.finish()
        for i, h in enumerate(histos):
            self.format_histo(h, lw=2, color=self.get_color())
            h.Scale(1 / h.GetMaximum())
            stack.Add(h)
            legend.AddEntry(h, '{0:06.1f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux().n), 'l')
            # legend.AddEntry(h, '{0:2d} V'.format(self.collection.values()[i].Bias), 'l')
        self.format_histo(stack, y_off=1.55, draw_first=True, x_tit='Pulse Height [au]', y_tit='Number of Entries')
        log_stack = stack.Clone()
        log_stack.SetMaximum(off)
        log_stack.SetNameTitle('hsdl', 'Signal Distribution LogY')
        self.save_histo(stack, 'SignalDistributions', False, self.save_dir, lm=.13, draw_opt='nostack', l=legend)
        self.save_histo(log_stack, 'SignalDistributionsLogY', False, self.save_dir, lm=.13, draw_opt='nostack', logy=True, l=legend)
        if show:
            c = TCanvas('c_sd1', 'Signal Distributions', 1500, 750)
            c.Divide(2)
            legends = [legend, legend.Clone()]
            for i, s in enumerate([stack, log_stack], 1):
                pad = c.cd(i)
                pad.SetLeftMargin(.14)
                pad.SetLogy() if i == 2 else do_nothing()
                s.Draw('nostack')
                legends[i - 1].Draw()
            self.RootObjects.append([c, legends])
        self.reset_colors()

    def draw_snrs(self, flux=True, draw=True):
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'SNR vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.collection.iteritems():
            snr = ana.calc_snr()
            x = ana.Run.Flux if flux else key
            gr.SetPoint(i, x, snr[0])
            gr.SetPointError(i, 0, snr[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        gr.Draw('ap')
        self.save_plots('AllSNRs', canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        gROOT.SetBatch(0)

    def draw_all_ph_distributions(self, binning=5000, show=False):

        pickle_path = self.FirstAnalysis.PickleDir + 'Ph_fit/Ph_distos_fits_{tc}_{rp}_{dia}_{bin}.pickle'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName, bin=binning)

        def func():
            collimator_settings = [(ana.Run.RunInfo['fs11'], ana.Run.RunInfo['fsh13']) for key, ana in self.collection.iteritems()]
            collimator_settings = set(collimator_settings)
            fits = {}
            for s in collimator_settings:
                retval = self.draw_ph_distributions(binning, show=show, fs11=s[0], fsh13=s[1])
                fits[s] = retval
            return fits

        res = func() if show else None
        return do_pickle(pickle_path, func, res)

    def get_repr_error(self, flux, show=True, redo=False):

        pickle_path = self.make_pickle_path('Errors', 'Repr', self.RunPlan, self.DiamondName, suf=flux)

        def f():
            runs = self.get_runs_below_flux(flux)
            if not runs:
                return 0
            vals = [make_ufloat(self.collection[run].draw_pulse_height(show=False, save=False)[1]) for run in runs]
            fluxes = [make_ufloat(self.get_fluxes()[run]) for run in runs]
            gr = self.make_tgrapherrors('gr_re', 'Pulse Heights Below {f} kHz/cm^{{2}}'.format(f=flux), x=fluxes, y=vals)
            set_statbox(entries=2, only_fit=True)
            gr.Fit('pol0', 'qs{s}'.format(s='' if show else '0'))
            self.format_histo(gr, x_tit='Flux [kHz/cm^{2}]', y_tit='Mean Pulse Height [au]', y_off=1.7)
            self.save_histo(gr, 'ReprErrors', show, draw_opt='ap', lm=.14, prnt=show)
            if len(vals) == 1:
                return .01  # take 1% if there is only one measurement below the given flux
            m, s = mean_sigma(vals)
            return s / m

        return do_pickle(pickle_path, f, redo=redo)

    def draw_ph_distributions(self, binning=5000, fsh13=.5, fs11=65, show=True):
        runs = self.get_runs_by_collimator(fsh13=fsh13, fs11=fs11)
        return self.draw_combined_ph_distributions(runs, binning, show)

    def draw_ph_distributions_below_flux(self, binning=None, flux=150, show=True, save_plot=True):
        binning = int(self.FirstAnalysis.Run.n_entries * .3 / 60 if binning is None else binning)
        pickle_path = self.make_pickle_path('Ph_fit', 'PhDistoBel', self.RunPlan, self.DiamondName, suf='{bin}_{flux}'.format(bin=binning, flux=flux))

        def func():
            log_message('Getting representative errors')
            runs = self.get_runs_below_flux(flux)
            return self.draw_combined_ph_distributions(runs, binning, show)

        err = func() if save_plot else None
        return do_pickle(pickle_path, func, err)

    def draw_combined_ph_distributions(self, runs, binning=200, show=True):
        stack = THStack('s_phd', 'Pulse Height Distributions')
        self.reset_colors()
        self.start_pbar(len(runs))
        for i, run in enumerate(runs, 1):
            ana = self.collection[run]
            set_root_output(False)
            h = ana.draw_ph_distribution(show=False, binning=binning, fit=False, save=False)
            set_root_output(True)
            self.format_histo(h, fill_color=4000)
            h.SetStats(False)
            h.SetLineColor(self.get_color())
            stack.Add(h)
            self.ProgressBar.update(i)
        self.ProgressBar.finish()
        self.draw_histo(stack, '', show, draw_opt='')
        summed_stack = stack.GetStack().Last().Clone()
        summed_stack.SetFillStyle(0)

        fit = TF1('fit', 'gaus', 0, 160)
        fit.SetNpx()
        fitptr = summed_stack.Fit(fit, 'SQ0', 'same')
        summed_stack.SetStats(1)
        set_statbox(only_fit=True)
        summed_stack.SetName('Fit Results')
        summed_stack.Draw('sames')
        stack.GetXaxis().SetRange(summed_stack.FindFirstBinAbove(0) - 2, summed_stack.FindLastBinAbove(0) + 2)
        fit.Draw('same')
        self.ROOTObjects.append(summed_stack)
        self.ROOTObjects.append(fit)
        self.save_plots('PulseHeightDistributions')
        return fitptr.Parameter(1), fitptr.Parameter(2), fitptr.Chi2() / fitptr.Ndf()

    def calc_pedestal_spread(self, fsh13=.5, fs11=65):
        runs = self.get_runs_by_collimator(fs11, fsh13)
        values = []
        for run in runs:
            ana = self.collection[run]
            fit = ana.show_pedestal_histo(save=False)
            values.append(fit.Parameter(1))
        return max(values) - min(values)

    def calc_all_pedestal_spreads(self, recalc=False):

        pickle_path = self.FirstAnalysis.PickleDir + 'Pedestal/PedSpread_{tc}_{rp}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName)

        def func():
            collimator_settings = [(ana.Run.RunInfo['fs11'], ana.Run.RunInfo['fsh13']) for key, ana in self.collection.iteritems()]
            collimator_settings = set(collimator_settings)
            mins = {}
            for s in collimator_settings:
                retval = self.calc_pedestal_spread(s[1], s[0])
                mins[s] = retval
            return mins

        res = func() if recalc else None
        return do_pickle(pickle_path, func, res)

    def calc_full_pedestal_spread(self):
        values = []
        for ana in self.collection.values():
            fit = ana.show_pedestal_histo(save=False)
            values.append(fit.Parameter(1))
        return max(values) - min(values)

    # endregion

    def __draw_signal_legend(self):
        sig = 'positive' if self.FirstAnalysis.Polarity > 0 else 'negative'
        l1 = self.make_legend(.17, .88, nentries=2, margin=.05, clean=True, x2=.5)
        l1.AddEntry(0, 'Signal Polarity:', '')
        l1.AddEntry(0, sig, '').SetTextAlign(12)
        l1.AddEntry(0, 'Pedestal Substraction:', '')
        l1.AddEntry(0, 'yes', '').SetTextAlign(12)
        l1.SetNColumns(2)
        l1.Draw()
        self.ROOTObjects.append(l1)

    def compare_pedestals(self, show=True):
        gr1 = self.draw_pedestals(show=False)
        gr2 = self.draw_pedestals(cut='pulser', show=False)
        graphs = [gr1, gr2]
        margins = find_graph_margins(graphs)
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Pedestal Comparison', 1000, 1000)
        c.SetLogx()
        legend = TLegend(.7, .78, .88, .88)
        names = ['Signal', 'Pulser', 'BeamOff']
        for i, gr in enumerate(graphs):
            gr.GetYaxis().SetRangeUser(*margins)
            self.format_histo(gr, color=self.get_color())
            gr.Draw('lp') if i else gr.Draw('alp')
            legend.AddEntry(gr, names[i], 'pl')
        legend.Draw()
        gROOT.SetBatch(0)
        self.save_plots('PulserPedestalComparison', sub_dir=self.save_dir)
        self.RootObjects.append([c, graphs, legend])

    def draw_beam_current(self, rel_t=True, show=True):
        graphs = [ana.draw_beam_current(show=False) for ana in self.collection.itervalues()]
        xvals = array([g.GetX()[i] for g in graphs for i in xrange(g.GetN())])
        xvals += (xvals[1] - xvals[0]) / 2
        h = TH1F('hbca', 'Beam Current for Run Plan {n}'.format(n=self.RunPlan), len(xvals) - 1, xvals)
        ibin = 1
        for g in graphs:
            for i in xrange(g.GetN()):
                h.SetBinContent(ibin, g.GetY()[i])
                ibin += 1
        self.format_histo(h, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', y_off=.85, fill_color=self.FillColor, stats=0, markersize=.3,
                          t_ax_off=self.FirstAnalysis.Run.StartTime if rel_t else 0)
        self.save_histo(h, 'AllBeamRate', show=show, draw_opt='hist', x_fac=1.5, y_fac=.75, lm=.065)

    # ============================================
    # region CUTS

    def draw_bucket_info(self, flux=True, show=True, mean_=True):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        gr1 = self.make_tgrapherrors('gr1', '', color=self.get_color())
        prefix = 'Number of Bucket Cut Events' if not mean_ else 'Mean Pulse Height with Different Bucket Cuts'
        gr2 = self.make_tgrapherrors('gr2', '{pref} vs {mod}'.format(pref=prefix, mod=mode), color=self.get_color())
        gr3 = self.make_tgrapherrors('gr3', '', color=self.get_color())
        i = 0
        for key, ana in self.collection.iteritems():
            x = ana.Run.Flux if flux else key
            if not mean_:
                n = ana.show_bucket_numbers(show=False)
                gr1.SetPoint(i, x, n['new'] / n['all'] * 100)
                gr2.SetPoint(i, x, n['old'] / n['all'] * 100)
            else:
                info = ana.show_bucket_means(show=False, plot_histos=False)
                gr1.SetPoint(i, x, info['new'][0])
                gr2.SetPoint(i, x, info['old'][0])
                gr3.SetPoint(i, x, info['no'][0])
                gr1.SetPointError(i, 0, info['new'][1])
                gr2.SetPointError(i, 0, info['old'][1])
                gr3.SetPointError(i, 0, info['no'][1])
            i += 1
        c = TCanvas('c', 'Bucket Numbers', 1000, 1000)
        c.SetLeftMargin(.13)
        if flux:
            c.SetLogx()
        self.format_histo(gr2, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Events [%]' if not mean_ else 'Mean [au]', y_off=1.7, color=None)
        gr2.Draw('apl')
        gr1.Draw('pl')
        if mean_:
            gr3.Draw('pl')
        leg = TLegend(.2, .8, .35, .9)
        leg.AddEntry(gr2, 'old cut', 'pl')
        leg.AddEntry(gr1, 'new cut', 'pl')
        if mean_:
            leg.AddEntry(gr3, 'no bucket', 'pl')
        leg.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('{mode}_' + mode, sub_dir=self.save_dir)
        self.RootObjects.append([c, gr1, gr2, gr3, leg])

    # endregion

    # ============================================
    # region PEAK VALUES
    def draw_signal_peaks(self, flux=True, draw=True, pulser=False):
        """
        Shows the means of the signal peak distribution.
        :param flux:
        :param draw:
        :param pulser:
        """
        mode = 'Flux' if flux else 'Run'
        signal = 'Pulser' if pulser else 'Signal'
        prefix = 'Mean of {pul} Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.Bias, pul='Pulser' if pulser else 'Signal')
        gr = self.make_tgrapherrors('gr', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_peak_values(draw=False, pulser=pulser)
            x = ana.Run.Flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(1))
            gr.SetPointError(i, 0, fit.ParError(1))
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'Mean of {0} Peaks'.format(signal), 1000, 1000)
        c.SetLogx()
        gr.Draw('alp')
        if not draw:
            c.Close()
        self.RootObjects.append([gr, c])

    def draw_signal_fwhm(self, flux=True, draw=True):
        """
        Shows the FWHM of the signal peak distribution.
        :param flux:
        :param draw:
        """
        mode = 'Flux' if flux else 'Run'
        prefix = 'FWHM of Signal Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.Bias)
        gr = self.make_tgrapherrors('gr1', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fwhm = ana.calc_peak_value_fwhm()
            x = ana.Run.Flux if flux else key
            gr.SetPoint(i, x, fwhm)
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'FWHM of Signal Peaks', 1000, 1000)
        gr.Draw('alp')
        if not draw:
            c.Close()
        self.RootObjects.append([gr, c])

    # endregion

    # ============================================
    # region 2D SIGNAL MAP
    def draw_mean_fwhm(self, saveplots=True, flux=True, draw=True):
        """
        Creates the FWHM Distribution of all selected MeanSignalHistograms
        :param saveplots: if True saves the plot
        :param flux: draw vs flux if True else vs run
        :param draw:
        """
        if not draw:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        prefix = 'FWHM of Mean Signal Histogram: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.Bias)
        gr = self.make_tgrapherrors('pedestal', prefix)
        conversion_factor = 2 * sqrt(2 * log(2))  # sigma to FWHM
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_sig_map_disto()
            x = ana.Run.Flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(2) * conversion_factor)
            gr.SetPointError(i, 0, fit.ParError(2) * conversion_factor)
            i += 1
        c = TCanvas('c', 'FWHM', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='Flux [kHz/cm2]', y_tit='FWHM [au]', y_off=1.1)
        gr.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        if saveplots:
            self.save_plots('Mean_FWHM_' + mode, canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        return gr

    def draw_log_fluxes(self):
        g1 = self.make_tgrapherrors('g_fp', 'Number of Peaks', color=self.get_color())
        pixel_fluxes = [ana.Run.Flux / 1000. for ana in self.collection.values()]
        g2 = self.make_tgrapherrors('g_ff', 'Pixel Fast-OR', x=self.Runs, y=pixel_fluxes, color=self.get_color())
        for i, (run, ana) in enumerate(self.collection.iteritems()):
            flux, err = ana.Peaks.get_flux()
            g1.SetPoint(i, run, flux / 1000.)
            g1.SetPointError(i, 0, err / 1000.)
        mg = TMultiGraph('mg_ff', 'Flux Comparison')
        mg.Add(g1, 'pl')
        mg.Add(g2, 'pl')
        l = self.make_legend(nentries=3, x2=.4)
        l.AddEntry(g1, g1.GetTitle(), 'pl')
        l.AddEntry(g2, g2.GetTitle(), 'pl')
        names = ['1x1', '2x1', '2x2', '4x2', '4x4']
        self.format_histo(mg, x_tit='Pattern', y_tit='Flux [MHz/cm^{2}]', y_off=1.4, draw_first=True)
        for i, run in enumerate(self.Runs):
            bin_x = mg.GetXaxis().FindBin(run)
            mg.GetXaxis().SetBinLabel(bin_x, names[i])
        self.save_histo(mg, 'FluxComparison', draw_opt='a', lm=.12, bm=.2, l=l)

    def draw_fluxes(self, bin_width=5, rel_time=False, show=True):

        pickle_path = self.make_pickle_path('Flux', 'FullFlux', self.RunPlan, self.channel, bin_width)

        def f():
            if not self.FirstAnalysis.has_branch('rate'):
                return self.draw_log_flux(rel_t=rel_time, show=show)
            else:
                histos = [ana.draw_flux(bin_width=bin_width, rel_t=False, show=False) for ana in self.collection.itervalues()]
                h1 = TH1F('hff', 'Flux Profile', *self.get_fixed_time_binning(bin_width))
                ibin = 1
                for h in histos:
                    for jbin in xrange(h.GetNbinsX()):
                        h1.SetBinContent(ibin, h.GetBinContent(jbin))
                        ibin += 1
                    ibin += 1
                return h1

        hist = do_pickle(pickle_path, f)
        self.format_histo(hist, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', t_ax_off=self.StartTime if rel_time else 0, fill_color=self.FillColor, y_range=[1, 20000], stats=0)
        self.save_histo(hist, 'FluxEvo', x_fac=1.5, y_fac=.75, show=show, logy=True)
        return hist

    def draw_flux_hist(self, bin_size=1, show=True):
        self.Currents.find_data()
        c_time = self.Currents.Time
        p = TProfile('hpr', 'Flux Evolution', int((c_time[-1] - c_time[0]) / bin_size), c_time[0], c_time[-1])
        g = self.draw_fluxes(show=False)
        for i in xrange(g.GetN()):
            p.Fill(g.GetX()[i], g.GetY()[i])
        # self.format_histo(p, x_tit='Time [hh:mm]', y_tit='Current [nA]', y_off=.8, fill_color=self.FillColor, markersize=.7, stats=0, t_ax_off=0)
        self.format_histo(p, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', y_off=.8, markersize=.7, stats=0, t_ax_off=0)
        self.draw_histo(p, '', show, lm=.08, draw_opt='p', x=1.5, y=.75)
        return p

    def draw_current_flux(self, c_range=None, fit=True, show=True):
        fluxes = self.get_fluxes().values()
        currents = self.get_currents().values()
        g = self.make_tgrapherrors('gcf', 'Leakage Current vs. Flux')
        for i, (flux, current) in enumerate(zip(fluxes, currents)):
            if current is not None:
                g.SetPoint(i, flux.n, current.n)
                g.SetPointError(i, flux.s + flux.n * .05, current.s)
        c_range = [.1, max([c.n for c in currents if c is not None]) * 2] if c_range is None else c_range
        self.format_histo(g, x_tit='Flux [kHz/cm^{2}', y_tit='Current [nA]', y_off=1.3, y_range=c_range, draw_first=True)
        g.GetXaxis().SetLimits(1, 20000)
        self.draw_histo(g, 'FluxCurrent', show, lm=.13, draw_opt='ap', logx=True, logy=True, bm=.17)
        if fit:
            set_statbox(only_fit=True, y=.33, entries=6, w=.22)
            f = TF1('fcf', 'pol1', .1, 1e5)
            f.SetParLimits(0, .1, 5)
            f.SetParLimits(1, 1e-5, 5e-3)
            g.Fit('fcf', 'q')
        self.draw_irradiation(make_irr_string(self.selection.get_irradiation()))
        self.save_plots('FluxCurrent', show=show)
        return g

    def save_signal_maps(self, hitmap=False, redo=False):

        name = 'signal' if not hitmap else 'hit'
        log_message('Generating {s} maps!'.format(s=name))
        self.start_pbar(self.NRuns)
        histos = []
        for i, ana in enumerate(self.collection.values(), 1):
            histos.append(ana.draw_signal_map(show=False, hitmap=hitmap, redo=redo, prnt=False))
            self.ProgressBar.update(i)
        self.ProgressBar.finish()

        # find min/max
        glob_max = (int(max([h.GetMaximum() for h in histos])) + 5) / 5 * 5
        glob_min = int(min([h.GetMinimum() for h in histos])) / 5 * 5
        for i, h in enumerate(histos):
            if not hitmap:
                self.format_histo(h, z_range=[glob_min, glob_max])
            self.save_histo(h, '{n}Map{nr}'.format(nr=str(i).zfill(2), n=name.title()), show=False, ind=i, draw_opt='colz', rm=.16, lm=.12, prnt=False)  # theta 55, phi 20

    def draw_cumulative_map(self, chi2=None, res=1.5, hitmap=False, redo=False, cut=None, show=True):

        self.start_pbar(self.NRuns)
        hitmaps, histos = [], []
        for i, ana in enumerate(self.collection.values(), 1):
            ana.Cut.set_chi2(chi2) if chi2 else do_nothing()
            hitmaps.append(ana.draw_dia_hitmap(show=False, cut='' if cut is None and hitmap else cut, redo=redo, prnt=False, res=res))
            histos.append(ana.draw_signal_map(show=False, redo=redo, prnt=False, res=res, cut=cut)) if not hitmap else do_nothing()
            self.ProgressBar.update(i)
        self.ProgressBar.finish()
        h_all = TH2I('hchm', 'Cumulative Diamond Hit Map', *self.Plots.get_global_bins(res)) if hitmap else TProfile2D('hcsm', 'Cumulative Signal Map', *self.Plots.get_global_bins(res))
        if hitmap:
            for h in hitmaps:
                h_all.Add(h)
        else:
            for i, h in enumerate(histos):
                for xbin in xrange(h.GetNbinsX()):
                    for ybin in xrange(h.GetNbinsY()):
                        # weight by number of hits
                        h_all.Fill(h.GetXaxis().GetBinCenter(xbin), h.GetYaxis().GetBinCenter(ybin), h.GetBinContent(xbin, ybin), hitmaps[i].GetBinContent(xbin, ybin))
        self.save_histo(h_all, 'Cumulative{s}Map'.format(s='Hit' if hitmap else 'Signal'), show, lm=.12, rm=.16, draw_opt='colz')

    def draw_signal_spreads(self, flux=True, draw=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'Relative Spread vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.collection.iteritems():
            rel_spread = ana.calc_signal_spread()
            x = ana.Run.Flux if flux else key
            gr.SetPoint(i, x, rel_spread[0])
            gr.SetPointError(i, 0, rel_spread[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Relative Spread [%]')
        gr.Draw('ap')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('RelativeSpread', canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        gROOT.SetBatch(0)

    def show_peak_distribution(self, show=True):
        """ Shows the positions of the peaks of the 2D map. """
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        # create an overall VotingHistogram
        ana = self.get_first_analysis()
        ana.draw_mean_signal_distribution(show=False)
        extrema = Extrema2D(ana.SignalMapHisto, ana.MeanSignalHisto)
        h = extrema.create_voting_histo()
        for run, ana in self.collection.iteritems():
            if ana.IsAligned:
                ana.find_2d_extrema(histo=h, show=False)
            else:
                print 'Run {run} is not aligned...'.format(run=run)
        c = TCanvas('c', 'Voting Histos', 1600, 800)
        c.Divide(2, 1)
        # new_pal = ar.array('i', [kYellow, kYellow, kOrange, kOrange - 3, kOrange + 7, kRed])
        ex = [TExec('ex1', 'gStyle->SetPalette(56);'), TExec('ex2', 'gStyle->SetPalette(51)')]
        if show:
            gROOT.SetBatch(0)
        for i, histo in enumerate(h.itervalues(), 1):
            c.cd(i)
            histo.Draw('col')
            ex[i - 1].Draw()
            histo.Draw('colz same')
        self.save_plots('PeakDistribution', sub_dir=self.save_dir)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.RootObjects.append([c, ex])

    def draw_signal_map(self, show=True, low=1e20, fid=True, factor=1):
        suffix = '{fid}{fac}'.format(fid='Fid' if fid else '', fac=factor)
        pickle_path = self.make_pickle_path('SignalMaps', 'SigMaps', self.RunPlan, self.channel, suffix)

        def func():
            histos = [ana.draw_signal_map(show=False, marg=False, fid=fid, res=factor) for ana in self.collection.itervalues() if ana.Run.Flux < low]
            h = histos[0]
            sig_map = TH2F('h_sms', 'Combined Signal Maps', h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax(), h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax())
            for h in histos:
                sig_map.Add(h)
            sig_map.Scale(1. / len(histos))
            self.format_histo(sig_map, x_tit='track_x [cm]', y_tit='track_y [cm]', y_off=1.4, z_off=1.3, stats=0, z_tit='Pulse Height [au]')
            self.__adjust_sig_map(sig_map)
            self.save_histo(sig_map, 'CombinedSignalMaps', show, lm=.12, rm=.16, draw_opt='colz')
            return sig_map

        hist = do_pickle(pickle_path, func)
        if not gROOT.FindObject('h_sms'):
            gStyle.SetPalette(53)
            self.__adjust_sig_map(hist)
            self.draw_histo(hist, '',  show, lm=.12, rm=.16, draw_opt='colz')
        return hist

    def __adjust_sig_map(self, h):
        h.GetZaxis().SetRangeUser(0, 500)
        h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinCenter(h.FindFirstBinAbove(0)), h.GetXaxis().GetBinCenter(h.FindLastBinAbove(0)))
        h.GetYaxis().SetRangeUser(h.GetYaxis().GetBinCenter(h.FindFirstBinAbove(0, 2)), h.GetYaxis().GetBinCenter(h.FindLastBinAbove(0, 2)))
        h1 = TH1F('h_av', 'hav', 100, 1, h.GetBinContent(h.GetMaximumBin()))
        for i in xrange(1, h.GetNbinsX() * h.GetNbinsY() + 1):
            h1.Fill(h.GetBinContent(i))
        hmax = h1.GetMaximum()
        ph_min, ph_max = h1.GetBinCenter(h1.FindFirstBinAbove(hmax * .08)), h1.GetBinCenter(h1.FindLastBinAbove(hmax * .02))
        h.GetZaxis().SetRangeUser(ph_min, ph_max)
        self.draw_histo(h1)

    def draw_hit_map(self, show=True):
        pass

    # endregion

    # ====================================================================================
    # region BEAM PROFILE
    def draw_beam_info(self, mean_=True, flux=True, show=True, direction='x', fit_margin=.6):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        title = 'Mean' if mean_ else 'Sigma'
        gr = self.make_tgrapherrors('gr', '{tit} of the Beam Profile in {dir}'.format(tit=title, dir=direction.title()))
        i = 0
        for key, ana in self.collection.iteritems():
            x = ana.Run.Flux if flux else key
            fit = ana.fit_beam_profile(direction, show=False, fit_margin=fit_margin)
            par = 1 if mean_ else 2
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Beam Profile', 1000, 1000)
        c.SetLeftMargin(.125)
        c.SetLogx() if flux else do_nothing()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='{tit} [cm]'.format(tit=title), y_off=1.8)
        gr.Draw('alp')
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('BeamProfile{tit}{dir}{mar}'.format(tit=title, dir=direction.title(), mar=fit_margin * 100), sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        return gr

    def draw_xy_profiles(self, flux=True, show=True, fitx=.4, fity=.7):
        gr1 = self.draw_beam_info(mean_=True, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr2 = self.draw_beam_info(mean_=True, flux=flux, show=False, direction='y', fit_margin=fity)
        gr3 = self.draw_beam_info(mean_=False, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr4 = self.draw_beam_info(mean_=False, flux=flux, show=False, direction='y', fit_margin=fity)
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulse Height Distribution', 1500, 1500)
        c.Divide(2, 2)
        for i, gr in enumerate([gr1, gr2, gr3, gr4], 1):
            self.format_histo(gr, y_off=1.3)
            pad = c.cd(i)
            pad.SetLogx() if flux else do_nothing()
            pad.SetBottomMargin(.15)
            gr.Draw('alp')
        gROOT.SetBatch(0)
        self.save_plots('BeamProfileOverview', sub_dir=self.save_dir)
        self.RootObjects.append([c, gr1, gr2, gr3, gr4])

    def draw_beam_profiles(self, show=True, direction='x'):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.draw_beam_profile(mode=direction, show=False, fit=False) for ana in self.collection.itervalues()]
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        legend = TLegend(.4, .6 - self.get_number_of_analyses() * 0.03, .6, .6)
        for i, h in enumerate(histos):
            h.SetStats(0)
            self.normalise_histo(h)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux().n) + '^{2}', 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.save_plots('AllBeamProfiles{mod}'.format(mod=direction.title()), sub_dir=self.save_dir)
        self.RootObjects.append([c, legend] + histos)

    # endregion

    # ====================================================================================
    # region TRACKS
    def show_chi2s(self, mode=None, show=True, disto=False):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.reset_colors()
        histos = [ana.draw_chi2(mode=mode, show=False, prnt=False) for ana in self.collection.itervalues()]
        yq = zeros(1)
        cuts = []
        for h in histos:
            h.GetQuantiles(1, yq, array([.9]))
            cuts.append(yq[0])
        cut = min(cuts)
        xmax = max(cuts) * 1.1
        legend = self.make_legend(y2=.95, nentries=self.get_number_of_analyses() + 1)
        stack = THStack('hx2', '#chi^{{2}}{mode}'.format(mode=' in ' + mode if mode is not None else ''))
        ymax = 0
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color(), lw=2)
            self.normalise_histo(h, to100=True)
            stack.Add(h)
            ymax = max(ymax, h.GetBinContent(h.GetMaximumBin()))
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux().n), 'l')
        self.format_histo(stack, x_tit='#chi^{2}', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True)
        stack.GetXaxis().SetRangeUser(0, xmax)
        mode = '' if mode is None else mode
        self.draw_histo(stack, '', show, self.save_dir, lm=.15, draw_opt='nostack', l=legend)
        l = self.draw_vertical_line(cut, -1e9, 1e9, color=2, style=2)
        legend.AddEntry(l, 'cut: 90% quantile', 'l')
        if disto:
            nominal_chi2 = TF1('f', '[1]*ROOT::Math::chisquared_pdf(x, {ndf})'.format(ndf=4 if not mode else 2), 0, xmax)
            histos[0].Fit(nominal_chi2, 'q0')
            nominal_chi2.SetNpx(1000)
            nominal_chi2.Draw('same')
            self.RootObjects.append(nominal_chi2)
        self.save_plots('AllChi2{mod}'.format(mod=mode.title()))

    def draw_all_chi2s(self, show=True):
        self.show_chi2s(show=show)
        self.show_chi2s('x', show)
        self.show_chi2s('y', show)

    def show_angles(self, mode='x', show=True):
        self.reset_colors()
        histos = [ana.draw_angle_distribution(mode=mode, show=False, print_msg=False) for ana in self.collection.itervalues()]

        legend = self.make_legend(nentries=self.get_number_of_analyses())
        stack = THStack('has', 'Track Angles in {mode}'.format(mode=mode.title()))
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color())
            self.normalise_histo(h, to100=True)
            stack.Add(h)
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux().n), 'l')
        self.format_histo(stack, x_tit='Angle [deg]', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True)
        stack.GetXaxis().SetRangeUser(-3, 4)
        self.RootObjects.append(self.save_histo(stack, 'AllTrackAngles{mod}'.format(mod=mode.title()), show, self.save_dir, lm=.15, draw_opt='nostack', l=legend))

    def draw_both_angles(self, show=True):
        self.show_angles('x', show)
        self.show_angles('y', show)

    def show_angle_peaks(self, mode='x', sigma=False, flux=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        tit = '{mod} of Track Angle Distribution in {dir}'.format(mod='Sigma' if sigma else 'Mean', dir=mode.title())
        gr = self.make_tgrapherrors('gr', tit)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.calc_angle_fit(mode, show=False)
            x = ana.Run.Flux if flux else key
            par = 2 if sigma else 1
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Angle Peaks', 1000, 1000)
        c.SetLeftMargin(.12)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod='Flux' if flux else 'Run', unit=' [kHz/cm2]' if flux else ''), y_tit='Sigma [deg]' if sigma else 'Mean [deg]', y_off=1.5)
        gr.Draw('ap')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('TrackAngle{0}s{mod}'.format('Sigma' if sigma else 'Mean', mod=mode.title()), canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])

    # endregion

    def draw_currents(self, v_range=None, rel_time=False, averaging=1, with_flux=False, c_range=None, f_range=None, draw_opt='ap', show=True):
        self.Currents.draw_indep_graphs(rel_time=rel_time, v_range=v_range, averaging=averaging, with_flux=with_flux, c_range=c_range, f_range=f_range, show=show, draw_opt=draw_opt)

    def draw_run_currents(self):
        log_message('Generating currents ...')
        self.start_pbar(self.NRuns)
        for i, ana in enumerate(self.collection.itervalues(), 1):
            ana.draw_current(relative_time=False, show=False)
            ana.Currents.get_current()
            self.ProgressBar.update(i)
        self.ProgressBar.finish()

    def draw_chi2s(self):
        self.start_pbar(self.NRuns)
        log_message('Generating chi2s ...')
        for i, ana in enumerate(self.collection.itervalues(), 1):
            ana.draw_all_chi2(show=False, prnt=False)
            self.ProgressBar.update(i)
        self.ProgressBar.finish()

    def get_hv_device(self):
        return self.FirstAnalysis.Currents.Name

    def get_currents(self):
        return OrderedDict((key, ana.Currents.get_current()) for key, ana in self.collection.iteritems())

    def make_flux_table(self):
        # for ana in self.collection.itervalues():
        #     print ana.Run.RunInfo
        fluxes = OrderedDict(sorted({cs: [] for cs in set([(ana.Run.RunInfo['fs11'], ana.Run.RunInfo['fs13']) for ana in self.collection.itervalues()])}.iteritems()))
        for ana in self.collection.itervalues():
            fluxes[(ana.Run.RunInfo['fs11'], ana.Run.RunInfo['fs13'])].append('{0:3.1f}'.format(ana.Run.calc_flux()))
        first_col = [[''] * max([len(col) for col in fluxes.itervalues()])]
        header = ['FS11/FSH13'] + ['{0}/{1}'.format(make_col_str(head[0]), make_col_str(head[1])) for head in fluxes.iterkeys()]
        print make_latex_table(header=header, cols=first_col + fluxes.values(), endline=True),
        print make_latex_table_row(['Mean'] + ['{0:3.1f}'.format(calc_mean(flux)[0]) for flux in fluxes.values()]),
        print make_latex_table_row(['Sigma'] + ['{0:3.1f}'.format(calc_mean(flux)[1]) for flux in fluxes.values()]),
        print '\\bottomrule'

    def make_signal_analysis(self, saveplots=True):
        """
        Run all available signal analyises together and plot them in an overview.
        :param saveplots:
        """
        start_time = time()
        c = TCanvas('c', 'overview', 800, 1000)
        c.Divide(1, 3)
        plots = [self.draw_pulse_heights(show=False), self.draw_pedestals(show=False), self.draw_mean_fwhm(draw=False)]
        for i, plot in enumerate(plots, 1):
            pad = c.cd(i)
            pad.SetLogx()
            plot.Draw()
        if saveplots:
            self.save_plots('Overview Plot', sub_dir=self.save_dir, canvas=c)
        self.RootObjects.append(c)

        print '\nThe preanalysis for this selection took', print_elapsed_time(start_time)

    def get_systematic_error_table(self, latex=False):
        f = open('PlotsFelix/table_{tc}_{rp}_{dia}.txt'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName), 'w')
        l = '&\t' if latex else ''
        print 'Pulser', 'SignalFlux', 'Signal < 150', 'Signal < 80', 'Full Pulser Spread', 'Single Pulser Spreads'
        pulser_err = self.Pulser.calc_all_errors().values()
        p_mean = calc_mean([x[0] for x in pulser_err])
        out = '{0:2.2f}\t{l}{1:2.2f}\t{l}{2:2.2f}\t{l}'.format(min(x[1] / x[0] for x in pulser_err) * 100, max(x[1] / x[0] for x in pulser_err) * 100, p_mean[1] / p_mean[0] * 100, l=l)
        scan_errors = self.draw_all_ph_distributions(show=False).values()
        out += '{0:2.2f}\t{l}{1:2.2f}\t{l}'.format(min(x[1] / x[0] for x in scan_errors) * 100, max(x[1] / x[0] for x in scan_errors) * 100, l=l)
        flux_errors = self.draw_ph_distributions_below_flux(flux=150, show=False)
        out += '{0:2.2f}\t{l}'.format(flux_errors[1] / flux_errors[0] * 100, l=l)
        flux_errors = self.draw_ph_distributions_below_flux(flux=80, show=False)
        out += '{0:2.2f}\t{l}'.format(flux_errors[1] / flux_errors[0] * 100, l=l)
        out += '{0:2.2f}\t{l}'.format(self.calc_full_pedestal_spread(), l=l)
        spreads = self.calc_all_pedestal_spreads().values()
        try:
            min_spread = min(x for x in spreads if x)
        except ValueError:
            min_spread = 0
        out += '{0:2.2f}\t{l}{1:2.2f}'.format(min_spread, max(x for x in spreads), l=l)
        print out
        f.write(out)

    def get_runs_by_collimator(self, fs11=65, fsh13=.5):
        return [key for key, ana in self.collection.iteritems() if ana.Run.RunInfo['fs11'] == fs11 and ana.Run.RunInfo['fs13'] == fsh13]

    def get_runs_below_flux(self, flux):
        return [key for key, ana in self.collection.iteritems() if ana.Run.Flux < flux]

    def select_runs_in_range(self, start, stop):
        new_collection = OrderedDict()
        for key, ana in self.collection.iteritems():
            if start <= key <= stop:
                new_collection[key] = ana
        if not new_collection:
            print 'You did not select any run! No changes were made!'
        else:
            self.collection = new_collection

    def get_binning(self, evts_per_bin, t_bins=True):
        binnings = [ana.Plots.get_binning(evts_per_bin, time_bins=t_bins) for ana in self.collection.itervalues()]
        return [sum(ibin[0] for ibin in binnings), concatenate([ibin[1] for ibin in binnings])]

    def get_t_binning(self, evts_per_bin):
        binnings = [ana.get_time_bins(evts_per_bin) for ana in self.collection.itervalues()]
        t_array = []
        off = 0
        for i, tup in enumerate(binnings):
            for t_sec in tup[1]:
                t_array.append(t_sec + off)
            # correct if last time of the
            if i < self.NRuns - 1 and t_array[-1] > binnings[i + 1][1][0]:
                off = self.get_break_time(i) + t_array[-1] - binnings[i + 1][1][0]
        return len(t_array) - 1, array(t_array)

    def get_fixed_time_binning(self, bin_width):
        times = [t1 for ana in self.collection.itervalues() for t1 in ana.Plots.get_time_binning(bin_width)[1]]
        return len(times) - 1, array(times)

    def get_break_time(self, ind):
        return (self.get_ana(ind + 1).Run.LogStart - self.get_ana(ind).Run.LogStart - self.get_ana(ind).Run.Duration).total_seconds()

    def get_start_end_times(self):
        ts = [[ana.Run.StartTime, ana.Run.EndTime] for ana in self.collection.itervalues()]
        for i in xrange(1, len(ts)):
            if ts[i][0] < ts[i - 1][1]:
                ts[i][1] += ts[i - 1][1] + self.get_break_time(i - 1) - ts[i][0]
                ts[i][0] = ts[i - 1][1] + self.get_break_time(i - 1)
        return ts

    def set_channel(self, ch):
        """
        Sets the channels to be analysed by the SignalAnalysisobjects.
        :param ch: int (0 || 3)
        """
        for ana in self.collection.values():
            ana.set_channel(ch)

    def get_fluxes(self):
        flux = OrderedDict()
        for key, ana in self.collection.iteritems():
            flux[key] = ana.get_flux()
        return flux

    @staticmethod
    def get_mode(vs_time):
        return 'Time' if vs_time else 'Flux'

    def get_first_analysis(self):
        return self.collection.values()[0]

    def get_last_analysis(self):
        return self.collection.values()[-1]

    def get_ana(self, ind):
        return self.collection.values()[ind]

    def get_run_numbers(self):
        """ :return: sorted list of run numbers in AnalysisCollection instance """
        return sorted(self.collection.keys())

    def get_number_of_analyses(self):
        """ :return: number of analyses that the analysis collection object contains """
        return len(self.collection)

    def show_information(self):
        print "ANALYSIS COLLECTION INFO:"
        self.get_first_analysis().print_info_header()
        for ana in self.collection.itervalues():
            ana.print_information(header=False)

    def print_loaded(self):
        print '\033[1A\rRuns {0}-{1} were successfully loaded!{2}\n'.format(self.Runs[0], self.Runs[-1], 20 * ' ')

    def set_verbose(self, status):
        self.verbose = status
        self.VoltageScan.verbose = status
        for ana in self.collection.itervalues():
            ana.verbose = status
            ana.Pulser.verbose = status
            ana.Pedestal.verbose = status
            ana.Cut.verbose = status

    def get_irradiation(self):
        return self.FirstAnalysis.get_irradiation()

    def get_attenuator(self):
        return self.FirstAnalysis.get_attenuator()

    @staticmethod
    def make_x_tit(vs_time):
        return '{mod}{unit}'.format(mod='Time' if vs_time else 'Flux', unit=' [hh:mm]' if vs_time else ' [kHz/cm^{2}]')


if __name__ == '__main__':
    st = time()
    main_parser = ArgumentParser()
    main_parser.add_argument('runplan', nargs='?', default=3)
    # noinspection PyTypeChecker
    main_parser.add_argument('dia', nargs='?', default=1, type=int)
    # noinspection PyTypeChecker
    main_parser.add_argument('dia2', nargs='?', default=1, type=int)
    main_parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
    main_parser.add_argument('-t', '--tree', action='store_false')
    main_parser.add_argument('-r', '--runs', action='store_true')
    main_parser.add_argument('-d', '--draw', action='store_true')
    main_parser.add_argument('-rd', '--redo', action='store_true')
    args = main_parser.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    run_plan = args.runplan
    diamond = args.dia
    a = Elementary(tc, True)
    sel = RunSelection(testcampaign=tc, verbose=True)
    sel.select_runs_from_runplan(run_plan, ch=diamond) if not args.runs else sel.select_runs([int(args.runplan), int(args.dia if args.dia else args.runplan)], args.dia2 if args.dia2 else 1)
    print_banner('STARTING PAD-ANALYSIS COLLECTION OF RUNPLAN {0}'.format(run_plan))
    a.print_testcampaign()
    t = load_root_files(sel, args.tree)
    z = AnalysisCollection(sel, threads=t)
    z.print_loaded()
    print_elapsed_time(st, 'Instantiation')
    if args.runs:
        z.Currents.draw_indep_graphs()
        raw_input('Press any button to exit')
    if args.draw:
        z.draw_little_all(args.redo)
