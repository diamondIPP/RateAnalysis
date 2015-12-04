__author__ = 'micha'

# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH2D, gStyle, TF1, TH1F, TCut, gROOT, TLegend
from newAnalysis import Analysis
# from collections import OrderedDict
from array import array
from math import sqrt
from argparse import ArgumentParser


# ==============================================
# MAIN CLASS
# ==============================================
class SignalAnalysis(Analysis):

    def __init__(self, run, channel, binning=5000):
        Analysis.__init__(self, run)

        # main
        self.channel = channel
        self.run_number = self.run.run_number
        self.ch_cut = self.cut[self.channel]
        self.save_dir = '{tc}_{run}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run_number, dia=self.run.diamondname[self.channel])

        # stuff
        self.draw_option = 'COLZ'
        self.bin_size = binning
        self.binning = self.__get_binning()
        self.time_binning = self.get_time_binning()
        self.n_bins = len(self.binning)
        self.polarity = self.polarities[self.channel]

        # names
        self.signal_name = self.signal_names[self.channel]
        self.pedestal_name = self.pedestal_names[self.channel]
        
        # projection
        self.signal_projections = {}
        # graphs
        self.pulse_height = TGraphErrors()
        self.signals = {}
        self.tmp_histos = {}
        self.canvas = None
        # histos
        self.histo = None
        self.signaltime = None

    def draw_pulse_height(self):
        titlePulseHeight = 'Run{run}: {dia} Signal Time Evolution'.format(run=self.run_number, dia=self.run.diamondname[self.channel])
        self.pulse_height.SetNameTitle('graph', titlePulseHeight)
        mode = 'mean'
        empty_bins = 0
        count = 0
        for i in xrange(self.n_bins):
            self.signal_projections[i] = self.signaltime.ProjectionY(str(self.run_number) + str(self.channel) + "signalprojection_bin_" + str(i).zfill(2), i + 1, i + 1)
            self.signal_projections[i].SetTitle("Run{run}Ch{channel} Signal Projection of Bin {bin}".format(run=self.run_number, channel=self.channel, bin=i))
            self.signal_projections[i].GetXaxis().SetTitle("Signal ({signal})".format(signal=self.signalname))
            if self.signal_projections[i].GetEntries() > 0:
                if mode in ["mean", "Mean"]:
                    self.pulse_height.SetPoint(count, (self.time_binning[i] - self.run.startTime) / 60e3, self.signal_projections[i].GetMean())
                    self.pulse_height.SetPointError(count, 0, self.signal_projections[i].GetRMS() / sqrt(self.signal_projections[i].GetEntries()))
                elif mode in ["fit", "Fit"]:
                    self.signal_projections[i].GetMaximum()
                    maxposition = self.signal_projections[i].GetBinCenter(self.signal_projections[i].GetMaximumBin())
                    self.signal_projections[i].Fit("landau", "Q", "", maxposition - 50, maxposition + 50)
                    fitfun = self.signal_projections[i].GetFunction("landau")
                    mpv = fitfun.GetParameter(1)
                    mpverr = fitfun.GetParError(1)
                    self.pulse_height.SetPoint(count, (i + 0.5) * self.run.totalMinutes / self.n_bins, mpv)
                    self.pulse_height.SetPointError(count, 0, mpverr)
                count += 1
            else:
                empty_bins += 1
        print 'Empty proj. bins:\t', str(empty_bins) + '/' + str(self.n_bins)
        self.__format_signal_graph('mean', None)

        self.pulse_height.Draw('ALP')

    def get_polarity(self):
        self.tree.GetEntry(0)
        return self.tree.polarities[self.channel]

    def show_signal_histo(self):
        canvas = TCanvas('bla', 'blub', 1000, 1000)
        self.histo = TH1F('signal b2', 'signal without cuts', 400, -100, 300)
        canvas.cd()
        self.tree.Draw("{name}>>signal b2".format(name=self.signal_name))
        self.SavePlots('signal_distribution', 'png', canvas=canvas, subDir=self.save_dir)

    def compare_cuts(self):
        canvas = TCanvas('bla', 'blub', 1000, 1000)
        c2 = TCanvas('all', 'all', 1000, 1000)
        cuts = self.ch_cut.cut_strings
        all_cuts = TCut('all_cuts', '')
        gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
        i_color = 2
        histos = []
        legend = TLegend(0.7, 0.3, 0.9, .7)
        for key, value in cuts.iteritems():
            if str(value):
                all_cuts += value
                print 'saving plot', key
                canvas.cd()
                name = 'signal_distribution_{cut}'.format(cut=key)
                histo = TH1F('signal b2', '', 400, -100, 300)
                hist_name = 'signal with cut ' + key
                histo.SetTitle(hist_name)
                # histo.GetYaxis().SetRangeUser(0, 7000)
                self.tree.Draw("{name}>>signal b2".format(name=self.signal_name), value)
                self.SavePlots(name, 'png', canvas=canvas, subDir=self.save_dir)
                c2.cd()
                c2.Update()
                histo.SetLineColor(i_color)
                if i_color == 2:
                    histo.SetTitle('signal distribution with different cuts')
                    histo.SetLineWidth(3)
                    histo.SetStats(0)
                    histo.Draw()
                else:
                    histo.SetLineWidth(1)
                    histo.Draw('same')
                histos.append(histo)
                if i_color == 4:
                   i_color += 2
                else:
                    i_color += 1
                legend.AddEntry(histo, key)
        # save all cuts
        print 'saving plot all cuts'
        canvas.cd()
        histo = TH1F('signal b2', '', 400, -100, 300)
        histo.SetTitle('signal with all cuts')
        # histo.GetYaxis().SetRangeUser(0, 7000)
        self.tree.Draw("{name}>>signal b2".format(name=self.signal_name), all_cuts)
        self.SavePlots('signal_distribution_all_cuts', 'png', canvas=canvas, subDir=self.save_dir)
        canvas.Close()
        c2.cd()
        histo.SetLineColor(1)
        histo.SetLineWidth(1)
        legend.AddEntry(histo, 'all cuts')
        histo.Draw('same')
        legend.Draw()
        self.SavePlots('all', 'png', canvas=c2, subDir=self.save_dir)
        self.SavePlots('all', 'root', canvas=c2, subDir=self.save_dir)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")

    def show_pedestal_histo(self):
        self.tmp_histos['1'] = TH1F('ped1', 'pedestal', 100, -20, 20)
        self.tmp_histos['2'] = TH1F('ped2', 'pedestal cut', 100, -20, 20)
        h1 = self.tmp_histos['1']
        h2 = self.tmp_histos['2']
        self.canvas = TCanvas('bla', 'blub', 1500, 1000)
        canvas = self.canvas
        canvas.Divide(2, 2)
        canvas.cd(1)
        self.tree.Draw("{name}>>ped1".format(name=self.pedestal_name), '', 'goff')
        h1.Scale(1 / h1.Integral(), 'width')
        h1.Draw()
        canvas.cd(2)
        self.tree.Draw("{name}>>ped2".format(name=self.pedestal_name), self.ch_cut.all_cut, 'goff')
        h2.Scale(1 / h2.Integral(), 'width')
        h2.Draw()
        canvas.cd(3)
        h2.SetLineColor(2)
        h1.Draw()
        h2.Draw('same')
        self.SavePlots('pedestal_cut', 'root', canvas=h2, subDir=self.save_dir)
        self.SavePlots('pedestal', 'root', canvas=h1, subDir=self.save_dir)
        canvas.Update()

    def __get_binning(self):
        jumps = self.ch_cut.jump_ranges
        n_jumps = len(jumps['start'])
        bins = [0, self.GetMinEventCut()]
        ind = 0
        for start, stop in zip(jumps['start'], jumps['stop']):
            gap = stop - start
            # continue if first start and stop outside minevent
            if stop < bins[-1]:
                ind += 1
                continue
            # if there is a jump from the start
            if start < bins[-1] < stop:
                bins.append(stop)
                ind += 1
                continue
            # add bins until hit interrupt
            while bins[-1] + self.bin_size < start:
                bins.append(bins[-1] + self.bin_size)
            # two jumps shortly after one another
            if ind < n_jumps - 2:
                next_start = jumps['start'][ind + 1]
                next_stop = jumps['stop'][ind + 1]
                if bins[-1] + self.bin_size + gap > next_start:
                    print 'double gap', bins[-1]
                    gap2 = next_stop - next_start
                    bins.append(bins[-1] + self.bin_size + gap + gap2)
                else:
                    bins.append(bins[-1] + self.bin_size + gap)
            else:
                bins.append(bins[-1] + self.bin_size + gap)
            ind += 1
        return bins

    def get_time_binning(self):
        time_bins = []
        for event in self.binning:
            time_bins.append(self.run.GetTimeAtEvent(event))
        return time_bins
    
    def __format_signal_graph(self, mode, setyscale):
        fit = TF1('fpol0', 'pol0')
        self.pulse_height.Fit(fit, 'Q')
        self.signals["signal"] = fit.GetParameter(0)
        print 'signal:\t\t\t', self.signals['signal']
        gStyle.SetOptFit(1)
        self.pulse_height.GetXaxis().SetTitleOffset(0.7)
        self.pulse_height.GetXaxis().SetTitle("time / min")
        self.pulse_height.GetXaxis().SetTitleSize(0.06)
        self.pulse_height.GetXaxis().SetLabelSize(0.06)
        # self.pulse_height.GetXaxis().SetRangeUser(0, self.analysis.run.totalMinutes)
        if mode in ["mean", "Mean"]:
            yTitlestr = "Mean Signal ({signalname})".format(signalname=(self.signaldefinition[self.channel]))
        else:
            yTitlestr = "MPV of Signal fit ({signalname})".format(signalname=(self.signaldefinition[self.channel]))
        self.pulse_height.GetYaxis().SetTitleOffset(0.9)
        self.pulse_height.GetYaxis().SetTitleSize(0.06)
        self.pulse_height.GetYaxis().SetLabelSize(0.06)
        self.pulse_height.GetYaxis().SetTitle(yTitlestr)
        if setyscale is not None:
            self.pulse_height.GetYaxis().SetRangeUser(setyscale[0], setyscale[1])
            self.pulse_height.Draw()
            # self.signalTimeCanvas.Update()

    def make_histos(self):
        # 2D Histos
        name = "signaltime_" + str(self.run_number)
        xbins = array('d', self.time_binning)
        # self.signaltime = TH2D(name, "signaltime", self.n_bins, 0, self.analysis.run.totalTime, 1000, -50, 300)
        self.signaltime = TH2D(name, "signaltime", len(xbins) - 1, xbins, 1000, -50, 300)
        self.tree.Draw("{name}:time>>{histo}".format(histo=name, name=self.signal_name), self.GetCut(self.channel), self.draw_option)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=392, type=int)
    args = parser.parse_args()
    test_run = args.run
    print '\nAnalysing run', test_run, '\n'
    z = SignalAnalysis(test_run, 0)
