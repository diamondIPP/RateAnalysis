#!/usr/bin/env python
# --------------------------------------------------------
#       Class to draw the info legend for an analysis class
# created on Jan 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from os import chdir
from subprocess import check_output
from ROOT import gROOT
from helpers.draw import Draw, make_tc_str, timedelta, make_flux_string, make_irr_string, mean


class InfoLegend(object):
    def __init__(self, analysis=None):
        self.Ana = analysis
        self.ShowDate = analysis.MainConfig.getboolean('SAVE', 'date')
        self.ShowGit = analysis.MainConfig.getboolean('SAVE', 'git hash') and not self.ShowDate
        self.ShowInfo = analysis.MainConfig.getboolean('SAVE', 'info legend')

        self.Objects = []

    def __repr__(self):
        on = ['OFF', 'ON']
        return f'InfoLegend: legend {on[self.ShowInfo]}, git hash {on[self.ShowGit]}, date {on[self.ShowDate]}'
        
    def has_dut(self):
        return hasattr(self.Ana, 'DUT')
        
    def is_collection(self):
        return hasattr(self.Ana, 'Runs') or hasattr(self.Ana, 'IsCollection') and self.Ana.IsCollection

    def is_active(self):
        if self.is_collection():
            return self.Ana.Ana.LoadTree if hasattr(self.Ana, 'Ana') else self.Ana.LoadTree
        return hasattr(self.Ana, 'Run') and self.Ana.Run.Number is not None and bool(self.Ana.Tree.Hash() if hasattr(self.Ana, 'Tree') else self.Ana.Ana.Tree.Hash())

    def draw(self, canvas=None, all_pads=True, show=True):
        """
        Draws the run infos inside the canvas. If no canvas is given, it will be drawn into the active Pad.
        :param canvas: canvas or pad to draw the legend on (get last active of none is provided)
        :param all_pads: draw info legends in all subpads or just the provided canvas
        :param show: wether or not to draw the legend
        :return: [run info legend, git text]
        """
        if not self.is_active():
            return
        if canvas is not None:
            canvas.cd()
            if show and canvas.GetBottomMargin() < .105 and self.ShowInfo:
                canvas.SetBottomMargin(0.15)
        else:
            canvas = gROOT.GetSelectedPad()
            if not canvas:
                self.Ana.warning('Cannot access an active Pad')
                return

        run_str = self.get_run_string()
        info_str = self.get_info_string()
        # width = len(run_str) * min(canvas.GetWw(), canvas.GetWh()) * .01
        # if canvas.GetWw() > canvas.GetWh():
        width = float(max(len(run_str), len(info_str))) / canvas.GetWw() * Draw.Res / 1000 * 11
        legend = Draw.make_legend(x1=.005, y1=.003, w=width, nentries=3, clean=False, scale=.75, margin=.05, fix=True)

        legend.AddEntry(0, run_str, '')                # Run String
        legend.AddEntry(0, self.get_dia_string(), '')  # Detector and Test Campaign
        legend.AddEntry(0, info_str, '')               # Info String (bias, attenuator, irr)

        git_text = self.make_git_text()

        self.Objects.append([legend, git_text])

        pads = [i for i in canvas.GetListOfPrimitives() if i.IsA().GetName() == 'TPad'] if all_pads else [canvas]
        pads = [canvas] if not pads else pads
        for pad in pads:
            pad.cd()
            if self.ShowGit:
                git_text.Draw()
            if self.ShowDate:
                Draw.date(.995, .005, align=31, size=.02)
            if self.ShowInfo:
                legend.Draw()
            pad.Modified()
        canvas.Update()
        return legend, git_text

    def get(self, canvas=None):
        return self.draw(canvas, show=False)

    @staticmethod
    def make_git_text():
        chdir(Draw.Dir)
        txt = 'git hash: {ver}'.format(ver=check_output(['git', 'describe', '--always']).decode('utf-8').strip('\n'))
        return Draw.tlatex(.9, .02, txt, show=False, ndc=True, size=.02)

    def get_duration(self):
        dur = sum([ana.Run.Duration for ana in self.Ana.get_analyses()], timedelta()) if self.is_collection() else self.Ana.Run.Duration
        return dur - timedelta(microseconds=dur.microseconds)

    def get_run_string(self):
        run_str = 'Run{run}: {rate}, {dur}'.format(run=self.get_runnumber_str(), rate=self.get_rate_string(), dur=self.get_duration())
        run_str += '' if self.is_collection() else ' ({} evts)'.format(self.Ana.Run.NEvents)
        return run_str

    def get_runnumber_str(self):
        return 's {}-{}'.format(self.Ana.get_runs()[0], self.Ana.get_runs()[-1]) if self.is_collection() else ' {}'.format(self.Ana.Run.Number)

    def get_dia_string(self):
        dia_str = self.Ana.DUT.Name if self.has_dut() else ', '.join(dut.Name for dut in self.Ana.Run.DUTs)
        return 'Detector{b}: {d} ({tc})'.format(b='' if self.has_dut() else 's', tc=make_tc_str(self.Ana.TCString), d=dia_str)

    def get_rate_string(self):
        if self.is_collection():
            fluxes = [flux.n for flux in self.Ana.get_fluxes(pbar=False)]
            return '{} - {}'.format(make_flux_string(min(fluxes)), make_flux_string(max(fluxes))) if 'voltage' not in self.Ana.Type else make_flux_string(mean(fluxes), prec=0)
        else:
            return make_flux_string(self.Ana.Run.Flux.n)

    def get_info_string(self):
        voltage = '{0:+4.0f}V'.format(self.Ana.DUT.Bias) if self.has_dut() else '/'.join('{0:+4.0f}V'.format(dut.Bias) for dut in self.Ana.Run.DUTs)
        irradiation = make_irr_string(self.Ana.DUT.Irradiation[self.Ana.TCString]) if self.has_dut() else '/'.join(make_irr_string(dut.get_irradiation(self.Ana.TCString)) for dut in self.Ana.Run.DUTs)
        attenuator = 'Att: {}'.format(str(self.Ana.DUT.Attenuator)) if self.has_dut() and self.Ana.DUT.Attenuator else ''
        return 'Info: {v}, {i}{a}'.format(v=voltage, i=irradiation, a=', {}'.format(attenuator) if attenuator else '')
