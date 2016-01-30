# ====================================
# IMPORTS
# ====================================
from ROOT import TCanvas, TPad, TGaxis, TText, TGraph, TBox, TPaveText, TPaveLabel, kCanDelete
# import array
from numpy import array
# import functions
import os
import time
from CurrentInfo import KeithleyInfo

# ====================================
# CONSTANTS
# ====================================
marker_size = 0.5
pad_color = 0  # was 26
axis_title_size = 0.06
left_margin = 0.08

space = 0.015
title = 0.93
col_vol = 807
col_cur = 418


# ====================================
# HELPER FUNCTIONS
# ====================================


def write_run(run):
    if int(run) < 10:
        run = "00" + str(run)
    elif int(run) < 100:
        run = "0" + str(run)
    return str(run)


# ====================================
# ROOT CLASS
# ====================================
class RootGraphs:
    def __init__(self, currentinfo):
        self.infos = currentinfo
        # self.infos = KeithleyInfo('', '')

        # self.c = TCanvas("c", self.canvas_name(), 1500, 1000)
        # dividing canvas into pads
        # self.p5 = TPad("p5", "", space, space, 1 - space, title / 3 - space / 2)
        # self.p6 = TPad("p6", "", space, title / 3 + space / 2, 1 - space, title * 2 / 3 - space / 2)
        # self.p7 = TPad("p7", "", space, space, 1 - space, title - space / 2)
        # self.p8 = TPad("p8", "", 0, title, 1, 1)
        # self.draw_pads()
        # # title texts
        # self.t2 = TPaveText(0, 0, 1, 0.7, "NB")
        # self.t3 = TPaveText(0.5, 0, 1, 0.7, "NB")
        # self.draw_titles()
        # graphs
        self.CurrentGraph = None
        self.VoltageGraph = None
        self.make_graphs()
        self.Margins = self.get_margins()
        # graph pads
        self.VoltagePad = None
        self.CurrentPad = None
        self.TitlePad = None
        # self.p11 = {}
        # self.p22 = {}
        # self.p33 = {}
        # self.make_pads()
        # # second axis
        # self.a1 = {}
        # self.make_axis()
        # # box
        # self.b1 = self.make_box()
        # # pad titles
        # self.t1 = {}
        # self.make_pad_title()
        # # run lines
        # self.run_lines = []
        # self.pt = TPaveLabel(.01, .01, .1, .1, "", "NDC")
        # self.pt.ResetBit(kCanDelete)
        # self.pt.SetFillColor(0)
        # self.pt.SetShadowColor(0)
        # self.pt.SetBorderSize(0)
        self.histos = {}

    def draw_graphs(self):
        c = TCanvas('c', 'Keithley Currents for Run {0}'.format(self.infos.RunNumber), 1500, 1000)
        self.make_pads()
        self.draw_pads()

        self.VoltagePad.cd()
        self.draw_voltage_frame()
        self.VoltageGraph.Draw('p')
        a = self.make_voltage_axis()
        a.Draw()

        self.TitlePad.cd()
        t = self.make_pad_title()
        t.Draw()

        self.CurrentPad.cd()
        self.draw_current_frame()
        self.CurrentGraph.Draw('p')

        c.Update()

        self.histos[0] = [c, t, a]
        self.infos.save_plots('Currents', self.infos.analysis.save_dir)

    def make_pads(self):
        self.VoltagePad = self.make_tpad('p1', gridy=True, margins=[.08, .07, .15, .15])
        self.CurrentPad = self.make_tpad('p2', gridx=True, margins=[.08, .07, .15, .15], transparent=True)
        self.TitlePad = self.make_tpad('p3', transparent=True)

    def draw_pads(self):
        for p in [self.VoltagePad, self.TitlePad, self.CurrentPad]:
            p.Draw()

    # def init_loop(self):
    #     ind = 0
    #     for key in self.infos.keithleys:
    #         self.goto_pad(ind)
    #         self.p1[key].Draw()
    #         self.p3[key].Draw()
    #         self.p2[key].Draw()
    #         ind += 1

    # def main_loop(self):
    #     ind = 0
    #     self.goto_pad(ind)
    #     # first pad with voltage, frame and second axis
    #     # self.p1[key].Draw()
    #     self.p1[key].cd()
    #     self.draw_current_frame(key)
    #     self.g2[key].Draw("P")
    #     # self.a1[key].Draw()
    #     self.a1[key].DrawAxis(self.xmax[key], -1600, self.xmax[key], 1600, -1600, 1600, 510, "+L")
    #     # second pad with pad titles and box
    #     # self.p3[key].Draw()
    #     self.p3[key].cd()
    #     self.b1.Draw("L")
    #     self.t1[key].Draw()
    #     # third pad with current, frame and run lines
    #     # self.p2[key].Draw()
    #     self.p2[key].cd()
    #     self.draw_current_frame(key)
    #     # if not self.runmode and not self.infos.single_mode:
    #     #     self.make_lines(key)
    #     self.g1[key].Draw("P")
    #     self.pt.Draw()
    #     ind += 1
    #     self.c.Update()
    #     t = time.localtime()
    #     self.pt.SetLabel("%s" % time.strftime('%x - %H:%M:%S', t))
    #     self.c.Update()
    #
    # def update(self):
    #     self.c.Update()
    #
    # def goto_pad(self, index):
    #     if index == 0:
    #         self.p7.cd()
    #         return self.p7
    #     elif index == 1:
    #         self.p6.cd()
    #         return self.p6
    #     elif index == 2:
    #         self.p5.cd()
    #         return self.p5
    #
    # def make_lines(self, key):
    #     for i in range(int(self.infos.run_start), int(self.infos.run_stop) + 1):
    #         # run start
    #         ymin = self.ymin[key]
    #         ymax = self.ymax[key]
    #         start = functions.convert_time(self.infos.get_time(i, "start time"))
    #         a2 = TGaxis(start, ymin, start, ymax, ymin, ymax, 510, "+SU")
    #         tit = "run " + str(i) + "  "
    #         a2.SetLineColor(1)
    #         a2.SetTickSize(0)
    #         a2.SetLabelSize(0)
    #         a2.SetTitleSize(0.05)
    #         a2.SetTitleOffset(0.1)
    #         if self.infos.single_mode:
    #             a2.SetTitleSize(0.025)
    #             a2.SetTitleOffset(0.25)
    #             flux = str(self.infos.data[functions.convert_run(i)]["measured flux"])
    #             spaces = ""
    #             for j in range(5 - len(flux)):
    #                 spaces += " "
    #             tit = "flux " + flux + spaces
    #         a2.SetTitle(tit)
    #         # run stop
    #         stop = functions.convert_time(self.infos.get_time(i, "stop time"))
    #         a3 = TGaxis(stop, ymin, stop, ymax, ymin, ymax, 510, "-SU")
    #         a3.SetLineColor(13)
    #         a3.SetTickSize(0)
    #         a3.SetLabelSize(0)
    #         # draw only for runs longer than 4 minutes
    #         if stop - start > 20 * 60:
    #             a2.Draw()
    #             a3.Draw()
    #         self.c.Update()
    #         self.run_lines.append(a2)
    #         self.run_lines.append(a3)
    #
    def make_pad_title(self):
        text = 'Currents measured by {0}'.format(self.infos.Name)
        t1 = TText(0.08, 0.88, text)
        t1.SetTextSize(0.07)
        return t1

    # @staticmethod
    # def make_box():
    #     b1 = TBox(0, 0, 1, 1)
    #     b1.SetLineWidth(3)
    #     b1.SetFillStyle(4000)
    #     # b1.Draw("L")
    #     return b1
    #

    def get_margins(self):
        x = [min(self.infos.Time), max(self.infos.Time)]
        dx = .05 * (x[1] - x[0])
        y = [min(self.infos.Currents), max(self.infos.Currents)]
        dy = .01 * (y[1] - y[0])
        return {'x': [x[0] - dx, x[1] + dx], 'y': [y[0] - dy, y[1] + dy]}
    #
    # def set_canvas(self, color):
    #     self.c.SetFillColor(color)  # was 17
    #
    # def draw_pads(self):
    #     self.p5.Draw()
    #     self.p6.Draw()
    #     self.p7.Draw()
    #     self.p8.Draw()
    #
    # def draw_titles(self):
    #     self.p8.cd()
    #     sig_type = self.infos.type
    #     if sig_type != "signal" and sig_type != "pedestal" and sig_type != "pedastal":
    #         sig_type = "unknown"
    #     title_text1 = 'Overview of all runs of ' + self.infos.dia1 + ' & ' + self.infos.dia2
    #     title_text2 = ""
    #     if self.infos.single_mode:
    #         title_text1 = str(self.infos.start) + " - " + str(self.infos.stop)
    #         title_text2 = ""
    #     elif not self.infos.time_mode:
    #         title_text1 = "PSI" + self.infos.start_run + " - " + sig_type + " run"
    #         title_text2 = "Flux = " + str(self.infos.flux) + " kHz/cm^{2}"
    #         if self.infos.flux == -1:
    #             title_text2 = "Flux not measured"
    #         middle_run = str((int(self.infos.run_stop) + int(self.infos.run_start)) / 2 + 1)
    #         if not self.runmode:
    #             title_text1 = "PSI" + functions.convert_run(self.infos.run_start) + " - " + self.infos.run_stop
    #             title_text2 = 'all runs of ' + self.infos.get_info(middle_run, "diamond 1") + ' & ' + self.infos.get_info(middle_run, "diamond 2")
    #     self.t2.AddText(title_text1)
    #     self.t2.SetAllWith("", "Align", 11)
    #     self.t2.SetAllWith("", "Size", 0.5)
    #     self.t2.SetFillColor(0)
    #     self.t3.AddText(title_text2)
    #     self.t3.SetAllWith("", "Align", 31)
    #     self.t3.SetAllWith("", "Size", 0.5)
    #     self.t3.SetFillColor(0)
    #     self.t2.Draw()
    #     if not self.infos.time_mode:
    #         self.t3.Draw()

    def make_graphs(self):
        tit = ' measured by {0}'.format(self.infos.Name)
        x = array(self.infos.Time)
        # current
        y = array(self.infos.Currents)
        g1 = TGraph(len(x), x, y)
        self.infos.format_histo(g1, 'Current', 'Current' + tit, color=col_cur, markersize=marker_size)
        # voltage
        y = array(self.infos.Voltages)
        g2 = TGraph(len(x), x, y)
        self.infos.format_histo(g2, 'Voltage', 'Voltage' + tit, color=col_vol, markersize=marker_size)
        self.CurrentGraph = g1
        self.VoltageGraph = g2

    # def refresh_graphs(self):
    #     for key in self.infos.keithleys:
    #         x = array.array('d', self.infos.time_x[key])
    #         y = array.array('d', self.infos.current_y[key])
    #         pos = self.g1[key].GetN()
    #         for i in range(pos, len(x) - 1):
    #             self.g1[key].SetPoint(i, x[i + 1], y[i + 1])
    #         y = array.array('d', self.infos.voltage_y[key])
    #         pos = self.g2[key].GetN()
    #         for i in range(pos, len(x) - 1):
    #             self.g2[key].SetPoint(i, x[i + 1], y[i + 1])
    #

    @staticmethod
    def make_transparent(pad):
        pad.SetFillStyle(4000)
        pad.SetFillColor(0)
        pad.SetFrameFillStyle(4000)

    def make_tpad(self, name, tit='', fill_col=0, gridx=False, gridy=False, margins=None, transparent=False):
        margins = [.1, .1, .1, .1] if margins is None else margins
        p = TPad(name, tit, 0, 0, 1, 1)
        p.SetFillColor(fill_col)
        p.SetMargin(*margins)
        p.ResetBit(kCanDelete)
        if gridx:
            p.SetGridx()
        if gridy:
            p.SetGridy()
        if transparent:
            self.make_transparent(p)
        return p

    def make_voltage_axis(self):
        a1 = self.infos.make_tgaxis(self.Margins['x'][1], -1100, 1100, '#font[22]{Voltage [V]}', color=col_vol, offset=.6, tit_size=.05, line=False, opt='+L', width=2)
        a1.CenterTitle()
        a1.SetLabelSize(0.04)
        a1.SetLabelColor(col_vol)
        a1.SetLabelOffset(0.01)
        return a1

    def draw_current_frame(self):
        m = self.Margins
        h2 = self.CurrentPad.DrawFrame(m['x'][0], m['y'][0], m['x'][1], m['y'][1])
        # X-axis
        h2.GetXaxis().SetTitle("#font[22]{time [hh:mm]}")
        h2.GetXaxis().CenterTitle()
        h2.GetXaxis().SetTimeFormat("%H:%M")
        h2.GetXaxis().SetTimeOffset(-3600)
        h2.GetXaxis().SetTimeDisplay(1)
        h2.GetXaxis().SetLabelSize(.04)
        h2.GetXaxis().SetTitleSize(axis_title_size)
        h2.GetXaxis().SetTitleOffset(1.05)
        h2.GetXaxis().SetTitleSize(0.05)
        # Y-axis
        h2.GetYaxis().SetTitleOffset(0.6)
        h2.GetYaxis().SetTitleSize(0.05)
        h2.GetYaxis().SetTitle("#font[22]{Current [nA]}")
        h2.GetYaxis().SetTitleColor(col_cur)
        h2.GetYaxis().SetLabelColor(col_cur)
        h2.GetYaxis().SetAxisColor(col_cur)
        h2.GetYaxis().CenterTitle()
        h2.GetYaxis().SetLabelSize(.04)
        h2.GetYaxis().SetTitleSize(axis_title_size)
        h2.GetYaxis().SetTitleOffset(.6)

    def draw_voltage_frame(self):
        m = self.Margins
        h1 = self.VoltagePad.DrawFrame(m['x'][0], -1100, m['x'][1], 1100)
        h1.SetTitleSize(axis_title_size)
        h1.GetXaxis().SetTickLength(0)
        h1.GetYaxis().SetTickLength(0)
        h1.GetXaxis().SetLabelOffset(99)
        h1.GetYaxis().SetLabelOffset(99)
        h1.GetYaxis().SetAxisColor(0)
        h1.SetLineColor(0)
    #
    # # save the root to file
    # def save_as(self, formats):
    #     # check if dir exists
    #     dirs = os.path.dirname("runs/")
    #     try:
    #         os.stat(dirs)
    #     except OSError:
    #         os.mkdir(dirs)
    #     run = 'run' + write_run(self.infos.run_start) + '-' + write_run(self.infos.run_stop)
    #     if self.infos.time_mode:
    #         run = 'all_' + self.infos.dia1 + '_' + self.infos.dia2 + '(' + write_run(self.infos.run_start) + '-' + write_run(self.infos.run_stop) + ')'
    #     elif self.infos.run_stop == '-1':
    #         run = 'run' + write_run(self.infos.run_start)
    #     if formats == "all":
    #         ftypes = [".pdf", ".eps", ".root"]
    #         for ftype in ftypes:
    #             filename = "runs/" + str(run) + ftype
    #             self.c.SaveAs(filename)
    #     else:
    #         filename = "runs/" + str(run) + "." + str(formats)
    #         self.c.SaveAs(filename)
    #
    # def save_single(self, formats):
    #     dirs = os.path.dirname("singleruns/")
    #     try:
    #         os.stat(dirs)
    #     except OSError:
    #         os.mkdir(dirs)
    #     if self.infos.run_stop == '-1':
    #         run = 'run' + write_run(self.infos.run_start)
    #     else:
    #         return
    #     filename = {'Keithley1': "singleruns/" + str(run) + '_' + str(self.infos.dia1) + "." + str(formats),
    #                 'Keithley2': "singleruns/" + str(run) + '_' + str(self.infos.dia2) + "." + str(formats)}
    #
    #     # single canvas margins
    #     x = int((1 - 2 * space) * 1000)
    #     y = int((title - 2 * space) / 2 * 1000)
    #     # gROOT.SetBatch(kFALSE)
    #     save_canvas = [TCanvas('save', 'save', x, y), TCanvas('save1', 'save', x, y)]
    #     # save_canvas = TCanvas('save', 'save', x, y)
    #
    #     ind = 0
    #     # remake the pads
    #     self.make_pads('save')
    #     # draw everythin in the save canvas cause ROOT is just too supid and can't save sub pads...
    #     for key in self.infos.keithleys:
    #         save_canvas[ind].cd()
    #         # first pad with voltage, frame and second axis
    #         self.p11[key].Draw()
    #         self.p11[key].cd()
    #         self.draw_frame1(key, save=True)
    #         self.g2[key].Draw("P")
    #         # self.a1[key].Draw()
    #         self.a1[key].DrawAxis(self.xmax[key], -1600, self.xmax[key], 1600, -1600, 1600, 510, "+L")
    #         # second pad with pad titles and box
    #         self.p33[key].Draw()
    #         self.p33[key].cd()
    #         self.t1[key].Draw()
    #         # third pad with current, frame and run lines
    #         self.p22[key].Draw()
    #         self.p22[key].cd()
    #         self.draw_frame2(key, save=True)
    #         # if not self.runmode and not self.infos.single_mode:
    #         #     self.make_lines(key)
    #         self.g1[key].Draw("P")
    #         self.pt.Draw()
    #         save_canvas[ind].SaveAs(filename[key])
    #         ind += 1
