#!/usr/bin/env python
# --------------------------------------------------------
#       Class to draw the info legend for an analysis class
# created on Jan 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gROOT, TLegend
from Utils import make_tc_str, timedelta, make_rate_str
from subprocess import check_output


class InfoLegend:
    def __init__(self, analysis):
        self.Analysis = analysis
        self.ShowGit = 1
        # self.ShowGit = analysis.MainConfigParser.getboolean('SAVE', 'git_hash')
        # self.ShowInfo = analysis.MainConfigParser.getboolean('SAVE', 'info_legend')
        self.ShowInfo = 1
        self.IsCollection = hasattr(analysis, 'collection')

        self.Objects = []

    def draw(self, canvas=None, all_pads=True, both_dias=False, show=True):
        """
        Draws the run infos inside the canvas. If no canvas is given, it will be drawn into the active Pad.
        :param all_pads: sets if the legens shall be drawn in all subpads or just the provided canvas
        :return: [run info legend, git text]
        """
        if canvas is not None:
            canvas.cd()
            if show and canvas.GetBottomMargin() < .105:
                canvas.SetBottomMargin(0.15)
        else:
            canvas = gROOT.GetSelectedPad()
            if not canvas:
                self.Analysis.log_warning('Cannot access an active Pad')
                return

        run_str = self.get_run_string()
        dia_str = self.get_diamond_str(both_dias)
        # width = len(run_str) * min(canvas.GetWw(), canvas.GetWh()) * .01
        # if canvas.GetWw() > canvas.GetWh():
        width = float(max(len(run_str), len(dia_str))) / canvas.GetWw() * 10.5
        legend = self.Analysis.make_legend(.005, .1, y1=.003, x2=width, nentries=3, clean=False, scale=.75, margin=.05)

        legend.AddEntry(0, 'Test Campaign: {tc}'.format(tc=make_tc_str(self.Analysis.generate_tc_str())), '')   # Test Campaign
        legend.AddEntry(0, run_str, '')                                                                         # Run String
        legend.AddEntry(0, dia_str, '')                                                                         # Diamond String

        git_text = self.make_git_text()

        self.Objects.append([legend, git_text])

        if show:
            pads = [i for i in canvas.GetListOfPrimitives() if i.IsA().GetName() == 'TPad'] if all_pads else [canvas]
            pads = [canvas] if not pads else pads
            for pad in pads:
                pad.cd()
                if self.ShowGit:
                    git_text.Draw()
                if self.ShowInfo:
                    legend.Draw()
                pad.Modified()
            canvas.Update()
        return legend, git_text

    def get(self, canvas=None):
        return self.draw(canvas, show=False)

    @staticmethod
    def make_git_text():
        git_text = TLegend(.85, 0, 1, .025)
        git_text.AddEntry(0, 'git hash: {ver}'.format(ver=check_output(['git', 'describe', '--always'])), '')
        git_text.SetLineColor(0)
        return git_text

    def get_duration(self):
        return str(sum([ana.run.Duration for ana in self.Analysis.collection.itervalues()], timedelta())) if self.IsCollection else self.Analysis.run.Duration

    def get_run_string(self):
        run_str = 'Run{run}: {rate}, {dur}'.format(run=self.get_runnumber_str(), rate=self.get_rate_string(), dur=self.get_duration())
        run_str += '' if self.IsCollection else ' ({} evts)'.format(self.Analysis.run.n_entries)
        return run_str

    def get_runnumber_str(self):
        return 's {}-{}'.format(self.Analysis.Runs[0], self.Analysis.Runs[-1]) if self.IsCollection else ' {}'.format(self.Analysis.RunNumber)

    def get_rate_string(self):
        if self.IsCollection:
            fluxes = self.Analysis.get_fluxes().values()
            return '{} - {}'.format(make_rate_str(min(fluxes)), make_rate_str(max(fluxes)))
        else:
            return make_rate_str(self.Analysis.run.Flux)

    def get_diamond_str(self, both_dias):
        if both_dias:
            dias = str(['{dia} @ {bias:+2.0f}V'.format(dia=dia, bias=bias) for dia, bias in zip(self.Analysis.run.DiamondNames, self.Analysis.run.Bias)])
            return 'Diamonds: {dias}'.format(dias=dias.strip('[]').replace('\'', ''))
        else:
            run_info = self.Analysis.FirstAnalysis.RunInfo if self.IsCollection else self.Analysis.RunInfo
            att = ', Attenuator: {a}'.format(a=run_info['att_dia{n}'.format(n=self.Analysis.DiamondNumber)]) if 'att_dia1' in run_info else ''
            return 'Diamond: {dia} @ {bias:+}V{a}'.format(dia=self.Analysis.DiamondName, bias=self.Analysis.Bias, a=att)
