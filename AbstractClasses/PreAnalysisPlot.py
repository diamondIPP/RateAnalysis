# ==============================================
# IMPORTS
# ==============================================
import ROOT
from AbstractClasses.Elementary import Elementary
from ROOT import gROOT, TGraphErrors, TCanvas, TH2D, gStyle, TF1
import types as t


# ==============================================
# MAIN CLASS
# ==============================================
class PreAnalysisPlot(Elementary):
    """
    Produces 3 plots inside one canvas:
     - mean signal vs time distribution
     - 2d signal vs time distribution
     - pedestal vs time distribution

    Cuts:
     - cuts from config file "AnalysisConfig.cfg"
     - first XXX events excluded (defined in same cfg file)

     to remove also beam interruption events, call the function
     Analysis.RemoveBeamInterruptions() first
    """

    def __init__(self, analysis, channel, canvas=None, binning=5000):
        Elementary.__init__(self)
        self.analysis = analysis
        if canvas is None:
            canvas = TCanvas("signalTimeCanvasRun" + str(self.analysis.run.run_number) + "Ch" + str(channel), "signalTimeCanvas" + "Ch" + str(channel), 650, 700)
        self.signalTimeCanvas = canvas
        self.channel = channel
        self.binning = binning
        self.padymargins = {}
        self.signals = {}
        self.nBins = int(self.analysis.get_event_at_time(9999999)) / int(self.binning)
        # graphs
        self.pulseHeight = TGraphErrors()
        self.pedgraph = TGraphErrors()
        # 2D Histos
        self.signaltime = TH2D("signaltime" + str(self.analysis.run.run_number), "signaltime", self.nBins, 0, self.analysis.run.totalTime, 200, -100, 500)
        self.pedestaltime = TH2D('pedestaltime' + str(self.analysis.run.run_number), 'pedestaltime', self.nBins, 0, self.analysis.run.totalTime, 200, -100, 500)
        self.signaltime2d = TH2D('signaltime2d{run}{channel}'.format(run=self.analysis.run.run_number, channel=self.channel), 'signaltime2d', self.nBins, self.analysis.run.startEvent / 1000,
                                 self.analysis.run.endEvent / 1000, 300, 0, 500)
        # projection
        self.signalProjection = {}
        # pads
        self.signalPad = None
        self.pedPad = None
        self.pad = None

        self.checkList = {
            "madeGraphs": {
                0: False,
                3: False}
        }

    def _Fill2DHistos(self, draw_option):
        # pedestal
        self.analysis.run.tree.Draw(
            self.analysis.pedestalname + "[{channel}]:(time-{starttime})>>pedestaltime{run}".format(channel=self.channel, starttime=self.analysis.run.startTime,
                                                                                                    run=self.analysis.run.run_number),
            self.analysis.get_cut(self.channel), draw_option)
        # signal
        nEntries = self.analysis.run.tree.Draw(
            (self.analysis.signaldefinition[self.channel] + ":(time-{starttime})>>signaltime{run}").format(channel=self.channel, starttime=self.analysis.run.startTime,
                                                                                                           run=self.analysis.run.run_number),
            self.analysis.get_cut(self.channel), draw_option)
        # printable signal
        self.analysis.run.tree.Draw(
            (self.analysis.signaldefinition[self.channel] + ":(event_number)/1000>>signaltime2d{run}{channel}").format(run=self.analysis.run.run_number, channel=self.channel),
            self.analysis.get_cut(self.channel), draw_option)
        return nEntries

    def _FillGraphs(self, mode):
        # set titles
        titlePulseHeight = 'Run ' + str(self.analysis.run.run_number) + ': ' + str(self.analysis.run.diamondname[self.channel]) + ' Signal Time Evolution'
        pedGraphTitle = 'Run ' + str(self.analysis.run.run_number) + ': ' + str(self.analysis.run.diamondname[self.channel]) + ' Pedestal Time Evolution'
        self.pulseHeight.SetNameTitle('graph', titlePulseHeight)
        self.pedgraph.SetNameTitle('ped_graph', pedGraphTitle)
        self.checkList['madeGraphs'][self.channel] = True

        count = 0
        empty_bins = 0
        runnumber = self.analysis.run.run_number
        for i in xrange(self.nBins):
            self.signalProjection[i] = self.signaltime.ProjectionY(str(runnumber) + str(self.channel) + "signalprojection_bin_" + str(i).zfill(2), i + 1, i + 1)
            self.signalProjection[i].SetTitle("Run{run}Ch{channel} Signal Projection of Bin {bin}".format(run=runnumber, channel=self.channel, bin=i))
            self.signalProjection[i].GetXaxis().SetTitle("Signal ({signal})".format(signal=self.analysis.signalname))
            binProjection_ped = self.pedestaltime.ProjectionY("proY_ped", i + 1, i + 1)
            if self.signalProjection[i].GetEntries() > 0:
                if mode in ["mean", "Mean"]:
                    self.pulseHeight.SetPoint(count, (i + 0.5) * self.analysis.run.totalMinutes / self.nBins, self.signalProjection[i].GetMean())
                    self.pulseHeight.SetPointError(count, 0, self.signalProjection[i].GetRMS() / ROOT.TMath.Sqrt(self.signalProjection[i].GetEntries()))
                elif mode in ["fit", "Fit"]:
                    self.signalProjection[i].GetMaximum()
                    maxposition = self.signalProjection[i].GetBinCenter(self.signalProjection[i].GetMaximumBin())
                    self.signalProjection[i].Fit("landau", "Q", "", maxposition - 50, maxposition + 50)
                    fitfun = self.signalProjection[i].GetFunction("landau")
                    mpv = fitfun.GetParameter(1)
                    mpverr = fitfun.GetParError(1)
                    self.pulseHeight.SetPoint(count, (i + 0.5) * self.analysis.run.totalMinutes / self.nBins, mpv)
                    self.pulseHeight.SetPointError(count, 0, mpverr)
                self.pedgraph.SetPoint(count, (i + 0.5) * self.analysis.run.totalMinutes / self.nBins, binProjection_ped.GetMean())
                self.pedgraph.SetPointError(count, 0, binProjection_ped.GetRMS() / ROOT.TMath.Sqrt(binProjection_ped.GetEntries()))
                count += 1
            else:
                empty_bins += 1
        print 'Empty proj. bins:\t', str(empty_bins) + '/' + str(self.nBins)

    def _FormatSignalGraph(self, mode, setyscale):
        fit = TF1('fpol0', 'pol0')
        self.pulseHeight.Fit(fit, 'Q')
        self.signals["signal"] = fit.GetParameter(0)
        print 'signal:\t\t\t', self.signals['signal']
        gStyle.SetOptFit(1)
        self.pulseHeight.GetXaxis().SetTitleOffset(0.7)
        self.pulseHeight.GetXaxis().SetTitle("time / min")
        self.pulseHeight.GetXaxis().SetTitleSize(0.06)
        self.pulseHeight.GetXaxis().SetLabelSize(0.06)
        self.pulseHeight.GetXaxis().SetRangeUser(0, self.analysis.run.totalMinutes)
        if mode in ["mean", "Mean"]:
            yTitlestr = "Mean Signal ({signalname})".format(signalname=(self.analysis.signaldefinition[self.channel]))
        else:
            yTitlestr = "MPV of Signal fit ({signalname})".format(signalname=(self.analysis.signaldefinition[self.channel]))
        self.pulseHeight.GetYaxis().SetTitleOffset(0.9)
        self.pulseHeight.GetYaxis().SetTitleSize(0.06)
        self.pulseHeight.GetYaxis().SetLabelSize(0.06)
        self.pulseHeight.GetYaxis().SetTitle(yTitlestr)
        if setyscale is not None:
            self.pulseHeight.GetYaxis().SetRangeUser(setyscale[0], setyscale[1])
            self.pulseHeight.Draw()
            self.signalTimeCanvas.Update()

    def _FormatPedGraph(self, setyscale):
        fit = TF1('fpol0', 'pol0')
        self.pedgraph.Fit(fit, 'Q')
        self.signals["pedestal"] = fit.GetParameter(0)
        print 'pedestal:\t\t', self.signals['pedestal']
        gStyle.SetOptFit(1)
        self.pedgraph.GetXaxis().SetTitleOffset(0.7)
        self.pedgraph.GetXaxis().SetTitle("time / min")
        self.pedgraph.GetXaxis().SetTitleSize(0.06)
        self.pedgraph.GetXaxis().SetLabelSize(0.06)
        self.pedgraph.GetXaxis().SetRangeUser(0, self.analysis.run.totalMinutes)
        yTitlestr = "Mean Pedestal ({pedestalname})".format(pedestalname=self.analysis.pedestalname + "[{channel}]".format(channel=self.channel))
        self.pedgraph.GetYaxis().SetTitleOffset(0.9)
        self.pedgraph.GetYaxis().SetTitleSize(0.06)
        self.pedgraph.GetYaxis().SetLabelSize(0.06)
        self.pedgraph.GetYaxis().SetTitle(yTitlestr)
        if setyscale is not None:
            self.pedgraph.GetYaxis().SetRangeUser(setyscale[0], setyscale[1])
            self.pedgraph.Draw()
            self.signalTimeCanvas.Update()

    def _Format2DHisto(self):
        gStyle.SetPalette(55)  # rainbow palette
        gStyle.SetNumberContours(200)
        # self.signaltime2d = gROOT.FindObject("signaltime2d{run}{channel}".format(run=self.analysis.run.run_number, channel=self.channel))
        self.signaltime2d.SetStats(0)
        self.signaltime2d.SetTitle("{signal} vs Event {cut}".format(signal=self.analysis.signaldefinition[self.channel],
                                                                    cut="{" + self.analysis.get_easy_cutstring(channel=self.channel) + "}"))
        self.signaltime2d.GetXaxis().SetLabelSize(0.06)
        self.signaltime2d.GetXaxis().SetTitle("event number / 1000")
        self.signaltime2d.GetXaxis().SetTitleSize(0.06)
        self.signaltime2d.GetXaxis().SetTitleOffset(0.7)
        self.signaltime2d.GetYaxis().SetTitleOffset(0.9)
        self.signaltime2d.GetYaxis().SetTitleSize(0.06)
        self.signaltime2d.GetYaxis().SetLabelSize(0.06)
        self.signaltime2d.GetYaxis().SetTitle(self.analysis.signaldefinition[self.channel])

    def _DrawFormattedGraphs(self, draw_option, setyscale_sig=None, setyscale_ped=None):
        if not bool(self.signalTimeCanvas):
            self.signalTimeCanvas = ROOT.TCanvas("signalTimeCanvas" + "Ch" + str(self.channel), "signalTimeCanvas" + "Ch" + str(self.channel), 650, 700)
            self.signalTimeCanvas.Divide(1, 3)

        self.signalTimeCanvas.cd(1)
        if setyscale_sig is not None:
            self.pulseHeight.GetYaxis().SetRangeUser(setyscale_sig[0], setyscale_sig[1])
            self.pulseHeight.Draw()
        else:
            self.pulseHeight.Draw()

        self.signalTimeCanvas.cd(2)
        self.signaltime2d.Draw(draw_option)

        self.signalTimeCanvas.cd(3)
        if setyscale_ped is not None:
            self.pedgraph.GetYaxis().SetRangeUser(setyscale_ped[0], setyscale_ped[1])
            self.pedgraph.Draw()
        else:
            self.pedgraph.Draw()

    def _SavePlots(self, savePlot):

        savename = "Run{run}_PreAnalysis_{diamond}".format(run=self.analysis.run.run_number, diamond=self.analysis.run.diamondname[self.channel])
        name_signal = 'signal_time'
        name_pedestal = 'pedestal_time'
        name_2D = '2D_distribution'
        subdir = '15080' + str(self.analysis.run.run_number) + '/' + str(self.analysis.run.diamondname[self.channel])
        # print "SAVENAME: ", savename
        if savePlot:
            print '\rSaving Plots for run', self.analysis.run.run_number,
            gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
            self.save_plots(savename, "png", canvas=self.signalTimeCanvas, sub_dir=self.analysis.run.diamondname[self.channel])
            self.save_plots(savename, "root", canvas=self.signalTimeCanvas, sub_dir="root")
            gROOT.SetBatch(1)
            saveCanvas = TCanvas('c1', 'c1', 1000, 600)
            saveCanvas.cd()
            self.pulseHeight.Draw("ALP")
            self.save_plots(name_signal, 'png', canvas=saveCanvas, sub_dir=subdir)
            self.pedgraph.Draw("ALP")
            self.save_plots(name_pedestal, 'png', canvas=saveCanvas, sub_dir=subdir)
            self.signaltime2d.Draw("colz")
            self.save_plots(name_2D, 'png', canvas=saveCanvas, sub_dir=subdir)
            saveCanvas.Close()
            gROOT.SetBatch(0)
        gROOT.ProcessLine("gErrorIgnoreLevel = 1;")

    def PrintInfos(self):
        print "Total Minutes:\t{tot}\nnbins:\t\t{nbins}".format(tot=self.analysis.run.totalMinutes, nbins=self.nBins)
        print "making PreAnalysis using\nSignal def:\t{signal}\nCut:\n\t{cut}".format(signal=self.analysis.signaldefinition[self.channel], cut=self.analysis.get_cut(self.channel))
        print "starttime:\t", self.analysis.run.startTime
        print "startevent:\t", self.analysis.run.startEvent
        print "endtime:\t", self.analysis.run.endTime
        print "endevent:\t", self.analysis.run.endEvent

    def Draw(self, mode="mean", savePlot=True, setyscale_sig=None, setyscale_ped=None):
        print mode
        assert (mode in ["mean", "Mean", "fit", "Fit"])
        assert (setyscale_sig is None or type(setyscale_sig) is t.ListType)

        drawOption2D = "COLZ"
        if not self.checkList['madeGraphs'][self.channel]:

            print '\nAnalysing run:', self.analysis.run.run_number

            # divide canvas
            self.signalTimeCanvas.Divide(1, 3)
            self.signalTimeCanvas.cd(1)

            nEntries = self._Fill2DHistos(drawOption2D)
            assert (int(nEntries) > 0), "Error: No signal event with current settings.. \nThe Cut is:\n\t" + self.analysis.get_cut(self.channel)

            self._FillGraphs(mode)

            # draw mean signal vs time
            self.signalPad = self.signalTimeCanvas.cd(1)
            self._FormatSignalGraph(mode, setyscale_sig)
            self.pulseHeight.Draw("ALP")
            self.signalPad.Update()
            # self.padymargins["signal"] = [self.signalPad.GetUymin(), self.signalPad.GetUymax()]

            # draw 2d distribution (high resolution)
            self.pad = self.signalTimeCanvas.cd(2)
            self.signalTimeCanvas.cd(2)
            self._Format2DHisto()
            self.signaltime2d.Draw(drawOption2D)
            self.analysis.draw_run_info(channel=self.channel, canvas=self.pad, infoid="preanalysis{run}{ch}".format(run=self.analysis.run.run_number, ch=self.channel))
            # self.pad.Update()

            # draw mean pedestal vs time
            self.pedPad = self.signalTimeCanvas.cd(3)
            self._FormatPedGraph(setyscale_ped)
            self.pedgraph.Draw("ALP")
            self.pedPad.Update()
            # self.padymargins["pedestal"] = [self.pedPad.GetUymin(), self.pedPad.GetUymax()]

            # self.signalTimeCanvas.Draw()

        else:
            self._DrawFormattedGraphs(drawOption2D, setyscale_sig, setyscale_ped)

        # update canvas
        self.signalTimeCanvas.Update()

        self._SavePlots(savePlot)

        # leave canvas open
        self.analysis.if_wait("showing MakePreAnalysis plots..")
