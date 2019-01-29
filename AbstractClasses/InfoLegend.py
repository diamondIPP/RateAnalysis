#!/usr/bin/env python
# --------------------------------------------------------
#       Class to draw the info legend for an analysis class
# created on Jan 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gROOT, TLegend
from Utils import make_tc_str, timedelta, make_rate_str, make_irr_string
from subprocess import check_output


class InfoLegend:
    def __init__(self, analysis):
        self.Analysis = analysis
        self.ShowGit = analysis.MainConfigParser.getboolean('SAVE', 'git_hash')
        self.ShowInfo = analysis.MainConfigParser.getboolean('SAVE', 'info_legend')
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
        info_str = self.get_info_string(both_dias)
        # width = len(run_str) * min(canvas.GetWw(), canvas.GetWh()) * .01
        # if canvas.GetWw() > canvas.GetWh():
        width = float(max(len(run_str), len(info_str))) / canvas.GetWw() * self.Analysis.Res / 1000 * 10.5
        legend = self.Analysis.make_legend(.005, .1, y1=.003, x2=width, nentries=3, clean=False, scale=.75, margin=.05)

        legend.AddEntry(0, run_str, '')                         # Run String
        legend.AddEntry(0, self.get_dia_string(both_dias), '')  # Detector and Test Campaign
        legend.AddEntry(0, info_str, '')                        # Info String (bias, attenuator, irr)

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
        dur = sum([ana.Run.Duration for ana in self.Analysis.collection.itervalues()], timedelta()) if self.IsCollection else self.Analysis.Run.Duration
        return dur - timedelta(microseconds=dur.microseconds)

    def get_run_string(self):
        run_str = 'Run{run}: {rate}, {dur}'.format(run=self.get_runnumber_str(), rate=self.get_rate_string(), dur=self.get_duration())
        run_str += '' if self.IsCollection else ' ({} evts)'.format(self.Analysis.Run.n_entries)
        return run_str

    def get_runnumber_str(self):
        return 's {}-{}'.format(self.Analysis.Runs[0], self.Analysis.Runs[-1]) if self.IsCollection else ' {}'.format(self.Analysis.RunNumber)

    def get_dia_string(self, both_dias):
        dia_str = ', '.join(self.Analysis.Run.DiamondNames) if both_dias else self.Analysis.DiamondName
        return 'Detector{b}: {d} ({tc})'.format(b='s' if both_dias else '', tc=make_tc_str(self.Analysis.TCString), d=dia_str)

    def get_rate_string(self):
        if self.IsCollection:
            fluxes = [flux.n for flux in self.Analysis.get_fluxes().values()]
            return '{} - {}'.format(make_rate_str(min(fluxes)), make_rate_str(max(fluxes)))
        else:
            return make_rate_str(self.Analysis.Run.Flux.n)

    def get_info_string(self, both_dias):
        voltage = ', '.join('{0:+4d}V'.format(i) for i in self.Analysis.Run.Bias) if both_dias else '{0:+4d}V'.format(self.Analysis.Bias)
        irradiation = ', '.join(make_irr_string(irr) for irr in (self.Analysis.Run.get_irradiations() if both_dias else [self.Analysis.get_irradiation()]))
        attenuator = 'Att: {}'.format(', '.join(self.Analysis.Run.get_attenuators()) if both_dias else str(self.Analysis.get_attenuator()))
        return 'Info: {v}, {i}, {a}'.format(v=voltage, i=irradiation, a=attenuator)
