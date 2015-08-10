import ROOT
from ROOT import gROOT


class PreAnalysisPlot(object):
    '''
    Produces 3 plots inside one canvas:
     - mean signal vs time distribution
     - 2d signal vs time distribution
     - pedestal vs time distribution

    Cuts:
     - cuts from config file "AnalysisConfig.cfg"
     - first XXX events excluded (defined in same cfg file)

     to remove also beam interruption events, call the function
     Analysis.RemoveBeamInterruptions() first
    '''
    
    def __init__(self, analysis, channel=0 ,canvas=None):
        self.analysis = analysis
        if canvas == None:
            canvas = ROOT.TCanvas("signalTimeCanvas"+"Ch"+str(channel), "signalTimeCanvas"+"Ch"+str(channel), 650, 700)
        self.signalTimeCanvas = canvas
        self.channel = channel
        
    def Draw(self):
        drawOption2D = "COLZ"
        nbins = 30

        #define graphs:
        graph = ROOT.TGraphErrors()
        graphtitle = 'Run{runnumber}: {diamond} Signal Time Evolution'.format(runnumber=self.analysis.run.run_number, diamond=self.analysis.run.diamondname[self.channel])
        graph.SetNameTitle('graph', graphtitle)
        pedgraph = ROOT.TGraphErrors()
        pedgraphtitle = 'Run{runnumber}: {diamond} Pedestal Time Evolution'.format(runnumber=self.analysis.run.run_number, diamond=self.analysis.run.diamondname[self.channel])
        pedgraph.SetNameTitle('ped_graph', pedgraphtitle)

        #set Canvas
        self.signalTimeCanvas.Divide(1,3)
        self.signalTimeCanvas.cd(1)

        #fill graph
        self.analysis.run.tree.GetEvent(0)
        startevent = self.analysis.run.tree.event_number
        starttime = self.analysis.run.tree.time
        self.analysis.run.tree.GetEvent(self.analysis.run.tree.GetEntries()-1)
        endevent = self.analysis.run.tree.event_number
        endtime = self.analysis.run.tree.time
        totalMinutes = (endtime-starttime)/60000.
        print "Total Minutes: {tot} nbins={nbins}".format(tot=totalMinutes, nbins=nbins)
        signaltime = ROOT.TH2D("signaltime" ,"signaltime", nbins, 0, (endtime-starttime), 600, -100, 500)
        pedestaltime = ROOT.TH2D("pedestaltime" ,"pedestaltime", nbins, 0, (endtime-starttime), 600, -100, 500)
        test = self.analysis.run.tree.Draw((self.analysis.signaldefinition+":(time-{starttime})>>signaltime").format(channel=self.channel, starttime=starttime), self.analysis.cut.format(channel=self.channel), drawOption2D, 10000000000, self.analysis.excludefirst)
        self.analysis.run.tree.Draw(self.analysis.pedestalname+"[{channel}]:(time-{starttime})>>pedestaltime".format(channel=self.channel, starttime=starttime), self.analysis.cut.format(channel=self.channel), drawOption2D, 10000000000, self.analysis.excludefirst)

        print "starttime: ", starttime
        print "startevent:", startevent
        print "endtime:", endtime
        print "endevent:", endevent

        assert(int(test)>0), "Error: No signal event with current settings.. \nThe Cut is:\n\t"+self.analysis.cut.format(channel=self.channel)

        count = 0
        final_i = 0
        for i in xrange(nbins):
            binProjection = signaltime.ProjectionY("proY", i+1,i+1)
            binProjection_ped = pedestaltime.ProjectionY("proY_ped", i+1,i+1)
            if binProjection.GetEntries() > 0:
                graph.SetPoint(count, (i+0.5)*totalMinutes/nbins, binProjection.GetMean())
                graph.SetPointError(count, 0, binProjection.GetRMS()/ROOT.TMath.Sqrt(binProjection.GetEntries()))
                pedgraph.SetPoint(count, (i+0.5)*totalMinutes/nbins, binProjection_ped.GetMean())
                pedgraph.SetPointError(count, 0, binProjection_ped.GetRMS()/ROOT.TMath.Sqrt(binProjection_ped.GetEntries()))
                count += 1
                final_i = i
            else:
                print "bin", i, " EMPTY"

        #draw mean signal vs time
        self.signalTimeCanvas.cd(1)
        graph.Fit("pol0")
        ROOT.gStyle.SetOptFit(1)
        graph.GetXaxis().SetTitleOffset(0.7)
        graph.GetXaxis().SetTitle("time / min")
        graph.GetXaxis().SetTitleSize(0.06)
        graph.GetXaxis().SetLabelSize(0.06)
        graph.GetXaxis().SetRangeUser(0, totalMinutes)
        yTitlestr = "Mean Signal ({signalname})".format(signalname=(self.analysis.signaldefinition.format(channel=self.channel)) )
        # graph.GetYaxis().SetRangeUser(ymin, ymax)
        graph.GetYaxis().SetTitleOffset(0.9)
        graph.GetYaxis().SetTitleSize(0.06)
        graph.GetYaxis().SetLabelSize(0.06)
        graph.GetYaxis().SetTitle(yTitlestr)
        graph.Draw("ALP")
        savename= "Run{runnumber}_{diamondname}_SignalTimeEvolution.eps".format(runnumber=self.analysis.run.run_number, diamondname=self.analysis.run.diamondname[self.channel])
        self.analysis.SavePlots(savename)

        #2d distribution (high resolution)
        self.signalTimeCanvas.cd(2)
        self.analysis.run.tree.Draw((self.analysis.signaldefinition+":(event_number)/1000>>signaltime2d({bins}, {start}, {end}, 300, 0, 500)").format(bins=200, channel=self.channel, start=startevent/1000, end=endevent/1000), self.analysis.cut.format(channel=self.channel), drawOption2D, 10000000000, self.analysis.excludefirst)
        signaltime2d = gROOT.FindObject("signaltime2d")
        signaltime2d.SetStats(0)
        signaltime2d.GetXaxis().SetLabelSize(0.06)
        signaltime2d.GetYaxis().SetLabelSize(0.06)
        signaltime2d.GetXaxis().SetTitle("event number / 1000")
        signaltime2d.GetXaxis().SetTitleSize(0.06)
        signaltime2d.GetXaxis().SetTitleOffset(0.7)
        signaltime2d.Draw(drawOption2D)

        #draw mean pedestal vs time
        self.signalTimeCanvas.cd(3)
        pedgraph.Fit("pol0")
        ROOT.gStyle.SetOptFit(1)
        pedgraph.GetXaxis().SetTitleOffset(0.7)
        pedgraph.GetXaxis().SetTitle("time / min")
        pedgraph.GetXaxis().SetTitleSize(0.06)
        pedgraph.GetXaxis().SetLabelSize(0.06)
        pedgraph.GetXaxis().SetRangeUser(0, totalMinutes)
        yTitlestr = "Mean Pedestal ({pedestalname})".format(pedestalname= self.analysis.pedestalname+"[{channel}]".format(channel=self.channel))
        # pedgraph.GetYaxis().SetRangeUser(ymin, ymax)
        pedgraph.GetYaxis().SetTitleOffset(0.9)
        pedgraph.GetYaxis().SetTitleSize(0.06)
        pedgraph.GetYaxis().SetLabelSize(0.06)
        pedgraph.GetYaxis().SetTitle(yTitlestr)
        pedgraph.Draw("ALP")
        #savename= "Run{runnumber}_{diamondname}_PedestalTimeEvolution.eps".format(runnumber=self.analysis.run.run_number, diamondname=self.analysis.run.diamondname[self.channel])
        #self.analysis.SavePlots(savename)


        #update canvas
        self.signalTimeCanvas.Update()
        self.analysis.SavePlots("Run{run}_PreAnalysis_{diamond}.png".format(run=self.analysis.run.run_number, diamond=self.analysis.run.diamondname[self.channel]))

        self.analysis.IfWait("showing MakePreAnalysis plots..")

