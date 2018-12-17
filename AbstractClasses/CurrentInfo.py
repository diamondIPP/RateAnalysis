# ====================================
# IMPORTS
# ====================================
from glob import glob
from json import load

from Elementary import Elementary
from RunSelection import RunSelection
from ROOT import TCanvas, TGraph, TProfile, TH1F, TH2F
from ConfigParser import ConfigParser
from argparse import ArgumentParser

from Utils import *

# ====================================
# CONSTANTS
# ====================================
axis_title_size = 0.06
label_size = .05
title_offset = 0.8
col_vol = 602  # 807
col_cur = 899  # 418
pad_margins = [.065, .09, .15, .1]
marker_size = .3


# ====================================
# CLASS FOR THE DATA
# ====================================
class Currents(Elementary):
    """reads in information from the keithley log file"""

    def __init__(self, analysis, averaging=False, points=10, start_run=None, verbose=False):
        self.Analysis = analysis
        self.IsCollection = hasattr(analysis, 'Runs')
        self.IsSelection = isinstance(analysis, RunSelection)
        Elementary.__init__(self, verbose=verbose)

        self.DoAveraging = averaging
        self.Points = points

        # config
        self.ConfigParser = self.load_parser()
        if self.ConfigParser is None:
            return

        # analysis/run info
        if self.IsCollection:
            self.RunPlan = analysis.RunPlan
            self.set_save_directory('Results/')
        self.RunLogs = self.load_runlogs()
        self.RunNumber = self.load_run_number()
        if not self.IsSelection:
            self.RunInfo = analysis.Run.RunInfo if not self.IsCollection else analysis.FirstAnalysis.RunInfo
            self.Channel = analysis.channel
            if 'dia1supply' not in self.RunInfo:
                return
        # todo: add a method to extract the currents for may
        self.DiamondName = self.load_dia_name()
        self.DiamondNumber = self.load_dia_number()
        self.Bias = self.load_bias()
        self.StartRun = start_run
        self.StartTime = self.load_start_time()
        self.StopTime = self.load_stop_time()

        # device info
        self.Number = self.get_device_nr()
        self.Channel = self.get_device_channel()
        self.Brand = self.ConfigParser.get('HV' + self.Number, 'name').split('-')[0].strip('0123456789')
        self.Model = self.ConfigParser.get('HV' + self.Number, 'model')
        self.Name = '{0} {1}'.format(self.Brand, self.Model)

        self.DataPath = self.find_data_path()
        self.OldDataPath = self.find_data_path(old=True)
        self.LogNames = None

        # data
        self.Currents = []
        self.IgnoreJumps = True
        self.Voltages = []
        self.Time = []
        self.MeanCurrent = 0
        self.MeanVoltage = 0
        self.NAveragedEvents = 0
        self.FoundStop = False

        # plotting
        self.CurrentGraph = None
        self.VoltageGraph = None
        self.Margins = None
        # graph pads
        self.VoltagePad = None
        self.CurrentPad = None
        self.TitlePad = None
        self.Histos = {}
        self.Stuff = []

    # ==========================================================================
    # region INIT
    def load_dia_name(self):
        return self.Analysis.DiamondName if not self.IsSelection else self.Analysis.SelectedDiamond

    def load_dia_number(self):
        if not self.IsSelection:
            return self.Analysis.DiamondNumber if not self.IsCollection else self.Analysis.FirstAnalysis.DiamondNumber
        return self.Analysis.SelectedDiamondNr

    def load_bias(self):
        if hasattr(self.Analysis, 'Type') and 'voltage' in self.Analysis.Type or self.IsSelection:
            return ''
        elif self.Analysis is not None:
            return self.Analysis.Bias

    def load_run_number(self):
        if not self.IsSelection:
            return self.Analysis.RunNumber if not self.IsCollection else self.Analysis.RunPlan

    def load_runlogs(self):
        filename = self.Analysis.Run.runinfofile if not self.IsCollection or self.IsSelection else self.Analysis.FirstAnalysis.Run.runinfofile
        try:
            f = open(filename)
            data = load(f)
            f.close()
        except IOError as err:
            log_warning('{err}\nCould not load default RunInfo!'.format(err=err))
            return None
        run_logs = OrderedDict(sorted(data.iteritems()))
        return run_logs

    def load_parser(self):
        parser = ConfigParser()
        if self.run_config_parser.has_option('BASIC', 'hvconfigfile'):
            file_path = self.run_config_parser.get('BASIC', 'hvconfigfile')
        else:
            file_path = join(self.DataDir, self.generate_tc_directory(), 'HV.cfg')
        if not file_exists(file_path):
            log_warning('HV info file "{f}" does not exist'.format(f=file_path))
            return None
        parser.read(file_path)
        self.log_info('HV Devices: {0}'.format([name for name in parser.sections() if name.startswith('HV')]))
        return parser

    def get_device_nr(self):
        if self.IsSelection:
            return str(self.Analysis.get_selection_runinfo().values()[0]['dia{}supply'.format(self.DiamondNumber)].split('-')[0])
        return str(self.RunInfo['dia{dia}supply'.format(dia=self.DiamondNumber)].split('-')[0])

    def get_device_channel(self):
        if self.IsSelection:
            full_str = self.Analysis.get_selection_runinfo().values()[0]['dia{}supply'.format(self.DiamondNumber)]
        else:
            full_str = self.RunInfo['dia{dia}supply'.format(dia=self.DiamondNumber)]
        return full_str.split('-')[1] if len(full_str) > 1 else '0'

    def find_data_path(self, old=False):
        if self.run_config_parser.has_option('BASIC', 'hvdatapath'):
            hv_datapath = self.run_config_parser.get('BASIC', 'hvdatapath')
        else:
            hv_datapath = join(self.DataDir, self.generate_tc_directory(), 'HVClient')
        if not dir_exists(hv_datapath):
            log_warning('HV data path "{p}" does not exist!'.format(p=hv_datapath))
        hv_datapath = join(hv_datapath, '{dev}_CH{ch}' if not old else '{dev}')
        return hv_datapath.format(dev=self.ConfigParser.get('HV' + self.Number, 'name'), ch=self.Channel)

    def load_start_time(self):
        if self.IsSelection:
            return
        ana = self.Analysis.FirstAnalysis if self.IsCollection else self.Analysis
        t = datetime.fromtimestamp(ana.Run.StartTime) if hasattr(ana.Run, 'StartTime') else ana.Run.LogStart
        return ana.Run.LogStart if t.year < 2000 or t.day != ana.Run.LogStart.day else t

    def load_stop_time(self):
        if self.IsSelection:
            return
        ana = self.Analysis.get_last_analysis() if self.IsCollection else self.Analysis
        t = datetime.fromtimestamp(ana.Run.EndTime) if hasattr(ana.Run, 'EndTime') else ana.Run.LogEnd
        return ana.Run.LogEnd if t.year < 2000 or t.day != ana.Run.LogEnd.day else t

    def set_start_stop(self, sta, sto=None):
        if not sta.isdigit():
            start_string = '{y}/{s}'.format(y=self.TESTCAMPAIGN[:4], s=sta)
            stop_string = '{y}/{e}'.format(y=self.TESTCAMPAIGN[:4], e=sto)
            self.StartTime = datetime.strptime(start_string, '%Y/%m/%d-%H:%M:%S')
            self.StopTime = datetime.strptime(stop_string, '%Y/%m/%d-%H:%M:%S')
        elif sto is None:
            self.StartRun = sta
            run = sta
            if not self.run_exists(run):
                return
            log = self.RunLogs[run]
            self.StartTime = datetime.strptime('{d} {t}'.format(d=log['begin date'], t=log['start time']), '%m/%d/%Y %H:%M:%S')
            self.StopTime = datetime.strptime('{d} {t}'.format(d=log['begin date'], t=log['stop time']), '%m/%d/%Y %H:%M:%S')
            if self.StartTime > self.StopTime:
                self.StopTime += timedelta(days=1)
        else:
            self.StartRun = sta
            log1 = self.RunLogs[sta]
            log2 = self.RunLogs[sto]
            try:
                self.StartTime = datetime.strptime('{d} {t}'.format(d=log1['begin date'], t=log1['start time']), '%m/%d/%Y %H:%M:%S')
                self.StopTime = datetime.strptime('{d} {t}'.format(d=log2['begin date'], t=log2['stop time']), '%m/%d/%Y %H:%M:%S')
            except KeyError:
                self.StartTime = datetime.strptime(log1['starttime0'], '%Y-%m-%dT%H:%M:%SZ') + timedelta(hours=1)
                self.StopTime = datetime.strptime(log2['endtime'], '%Y-%m-%dT%H:%M:%SZ') + timedelta(hours=1)

    def set_device(self):
        self.reset_data()
        self.Number = self.get_device_nr()
        self.Channel = self.get_device_channel()
        self.Brand = self.ConfigParser.get('HV' + self.Number, 'name').split('-')[0].strip('0123456789')
        self.Model = self.ConfigParser.get('HV' + self.Number, 'model')
        self.Name = '{0} {1}'.format(self.Brand, self.Model)
        self.DataPath = self.find_data_path()

    def reset_data(self):
        self.Currents = []
        self.Voltages = []
        self.Time = []
    # endregion

    # ==========================================================================
    # region DATA ACQUISITION
    def get_logs_from_start(self):
        log_names = sorted([name for name in glob(join(self.DataPath, '*'))] + [name for name in glob(join(self.OldDataPath, '*'))])
        start_log = None
        for i, name in enumerate(log_names):
            log_date = self.get_log_date(name)
            if log_date >= self.StartTime:
                break
            start_log = i
        self.log_info('Starting with log: {0}'.format(basename(log_names[start_log])))
        return log_names[start_log:]

    @staticmethod
    def get_log_date(name):
        log_date = basename(name).split('_')
        log_date = ''.join(log_date[-6:])
        return datetime.strptime(log_date, '%Y%m%d%H%M%S.log')

    def set_start(self, zero=False):
        self.Currents.append(self.Currents[-1] if not zero else 0)
        self.Voltages.append(self.Voltages[-1] if not zero else 0)
        self.Time.append(time_stamp(self.StartTime))

    def set_stop(self, zero=False):
        self.Currents.append(self.Currents[-1] if not zero else 0)
        self.Voltages.append(self.Voltages[-1] if not zero else 0)
        self.Time.append(time_stamp(self.StopTime))

    def find_data(self):
        if self.Currents:
            return
        self.LogNames = self.get_logs_from_start()
        for i, log_file_name in enumerate(self.LogNames):
            log_date = self.get_log_date(log_file_name)
            with open(log_file_name) as f:
                self.goto_start(f) if not i else do_nothing()  # jump to the correct line of the first file
                for line in f:
                    if not self.save_data(line, log_date):
                        break
            if self.FoundStop:
                break
        if self.Currents:
            self.set_stop()
        if not self.Currents:
            self.set_start(zero=True)
            self.set_stop(zero=True)
        if mean(self.Currents) < 0:
            self.Currents = [cur * -1 for cur in self.Currents]

    def save_data(self, line, log_date, shifting=False):
        info = line.split()
        if not isfloat(info[1]) or len(info) < 3:  # goto next line if there is device info in the line
            return True
        t = datetime.strptime('{}{}'.format(log_date.strftime('%Y%m%d'), info[0]), '%Y%m%d%H:%M:%S')
        current, voltage = float(info[2]) * 1e9, float(info[1])
        if t >= self.StopTime:
            self.FoundStop = True
            return False
        if self.StartTime < t < self.StopTime and current < 1e30:
            if self.DoAveraging:
                if not shifting:
                    self.MeanCurrent += current
                    self.MeanVoltage += voltage
                    self.NAveragedEvents += 1
                    if self.NAveragedEvents % self.Points == 0:
                        if mean(self.Currents) < 5 * self.MeanCurrent / self.Points:
                            self.Currents.append(self.MeanCurrent / self.Points)
                            self.Time.append(time_stamp(t))
                            self.Voltages.append(self.MeanVoltage / self.Points)
                            self.MeanCurrent = 0
                            self.MeanVoltage = 0
            else:
                if self.IgnoreJumps:
                    if len(self.Currents) > 100 and abs(self.Currents[-1] * 100) < abs(current):
                        if abs(self.Currents[-1]) > 0.01:
                            return True
                self.Currents.append(current)
                self.Time.append(time_stamp(t))
                self.Voltages.append(voltage)
        return True

    def goto_start(self, data):
        log_date = self.get_log_date(self.LogNames[0])
        lines = len(data.readlines())
        data.seek(0)
        if lines < 10000:
            return
        was_lines = 0
        for i in range(6):
            lines /= 2
            for j in xrange(lines):
                data.readline()
            while True:
                info = data.readline().split()
                if not info:
                    break
                if isfloat(info[1]):
                    now = datetime.strptime(log_date.strftime('%Y%m%d') + info[0], '%Y%m%d%H:%M:%S')
                    if now < self.StartTime:
                        was_lines += lines
                        break
                    else:
                        data.seek(0)
                        for k in xrange(was_lines):
                            data.readline()
                        break

    def convert_to_relative_time(self):
        zero = self.Time[0]
        for i in xrange(len(self.Time)):
            self.Time[i] = self.Time[i] - zero

    # endregion

    # ==========================================================================
    # region PLOTTING

    def draw_hist(self, bin_size=1, show=True):
        self.find_data()
        p = TProfile('hpr', 'Leakage Current', int((self.Time[-1] - self.Time[0]) / bin_size), self.Time[0], self.Time[-1])
        for t, c in zip(self.Time, self.Currents):
            p.Fill(t, c)
        self.format_histo(p, x_tit='Time [hh:mm]', y_tit='Current [nA]', y_off=.8, markersize=.7, stats=0, t_ax_off=0)
        self.draw_histo(p, '', show, lm=.08, draw_opt='p', x=1.5, y=.75)
        return p

    def draw_flux_correlation(self, show=True):
        p1 = self.draw_hist(show=False)
        p2 = self.Analysis.draw_flux_hist(show=False)
        ybins = log_bins(int(sqrt(p1.GetNbinsX())) * 3, .1, p1.GetMaximum() * 2)
        xbins = log_bins(int(sqrt(p1.GetNbinsX())), .1, p2.GetMaximum() * 2)
        h = TH2F('gfcc', 'Correlation of Flux and Current', *(xbins + ybins))
        for i in xrange(p1.GetNbinsX()):
            if p1.GetBinContent(i) and p2.GetBinContent(i):
                h.Fill(p2.GetBinContent(i), p1.GetBinContent(i))
        self.format_histo(h, y_tit='Current [nA]', x_tit='Flux [kHz/cm^{2}', stats=0, y_off=1.2)
        self.save_histo(h, 'FluxCurrent', draw_opt='colz', rm=.18, show=show, lm=.12, logy=True, logx=True)

    def set_graphs(self, averaging=1):
        self.find_data()
        sleep(.1)
        self.make_graphs(averaging)
        self.set_margins()

    def draw_distribution(self, show=True):
        self.find_data()
        m, s = mean_sigma(self.Currents)
        s = .1 if not s else s
        set_root_output(False)
        h = TH1F('hcd', 'Current Distribution', 5 * int(sqrt(len(self.Currents))), m - 2 * s, m + 2 * s)
        for current in self.Currents:
            h.Fill(current)
        self.format_histo(h, x_tit='Current [nA]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.draw_histo(h, '', show, lm=.13)
        return h

    def get_current(self):
        h = self.draw_distribution(show=False)
        if h.GetEntries() < 10:
            return None
        values = [h.GetBinCenter(i) for i in xrange(h.GetNbinsX())]
        weights = [h.GetBinContent(i) for i in xrange(h.GetNbinsX())]
        m, s = mean_sigma(values, weights)
        fit = h.Fit('gaus', 'sq0', '', m - s, m + s)
        fm, fs = fit.Parameter(1), fit.Parameter(2)
        if .8 * m < fit.Parameter(1) < 1.2 * m and s > 0 and fs < fm and fit.ParError(1) < m:  # only use gauss fit if its not deviating too much from the the mean
            return ufloat(fm, fs + .05 + .05 * fm)  # add .05 as uncertainty of the device and 5% systematic error
        return ufloat(h.GetMean(), h.GetMeanError() + .05 + .05 * h.GetMean())

    def draw_indep_graphs(self, rel_time=False, ignore_jumps=True, v_range=None, f_range=None, c_range=None, averaging=1, with_flux=False, draw_opt='ap', show=True):
        self.IgnoreJumps = ignore_jumps
        self.set_graphs(averaging)
        set_root_output(show)
        c = TCanvas('c', 'Keithley Currents for Run {0}'.format(self.RunNumber), int(self.Res * 1.5), int(self.Res * .75))
        self.draw_flux_pad(f_range, rel_time, draw_opt) if with_flux else self.draw_voltage_pad(v_range, draw_opt)
        self.draw_title_pad()
        self.draw_current_pad(rel_time, c_range, draw_opt)
        if self.IsCollection:
            self.draw_irradiation(make_irr_string(self.Analysis.selection.get_irradiation()))
        self.Stuff.append(c)
        run = self.Analysis.selection.SelectedRunplan if self.IsCollection else self.RunNumber
        save_name = 'Currents{}_{}_{}'.format(self.Analysis.TCString, run, self.DiamondNumber)
        self.save_canvas(c, name=save_name, sub_dir='Currents', show=show)

    def zoom_pads(self, low, high):
        self.VoltageGraph.GetXaxis().SetRangeUser(low, high)
        self.CurrentGraph.GetXaxis().SetRangeUser(low, high)

    def draw_current_pad(self, rel_t, c_range, draw_opt):
        self.draw_tpad('p3', gridx=True, margins=pad_margins, transparent=True)
        g = self.CurrentGraph
        self.format_histo(g, x_tit='#font[22]{Time [hh:mm]}', lab_size=label_size, x_off=1.05, tit_size=axis_title_size, t_ax_off=self.Time[0] if rel_t else 0, y_off=.55, yax_col=col_cur,
                          y_tit='#font[22]{Current [nA]}', center_y=True, x_range=[self.Time[0], self.Time[-1]], y_range=c_range, color=col_cur, markersize=marker_size)
        self.CurrentGraph.Draw(draw_opt)

    def draw_voltage_pad(self, v_range, draw_opt='ap'):
        self.draw_tpad('p1', gridy=True, margins=pad_margins, transparent=True)
        g = self.VoltageGraph
        v_range = [-1100, 1100] if v_range is None else v_range
        self.format_histo(g, y_range=v_range, y_tit='#font[22]{Voltage [V]}', x_range=[self.Time[0], self.Time[-1]], tit_size=axis_title_size, tick_size=0, x_off=99, l_off_x=99, center_y=True,
                          color=col_vol, y_off=title_offset, markersize=marker_size, yax_col=col_vol, lw=3, lab_size=label_size)
        g.Draw('{}y+'.format(draw_opt))

    def draw_flux_pad(self, f_range, rel_t=False, draw_opt='ap'):
        pad = self.draw_tpad('pr', margins=pad_margins, transparent=True, logy=True)
        h = self.Analysis.draw_fluxes(rel_time=rel_t, show=False)
        pad.cd()
        f_range = [1, 20000] if f_range is None else f_range
        self.format_histo(h, title=' ', y_tit='#font[22]{Flux [kHz/cm^{2}]}', fill_color=4000, fill_style=4000, lw=3, y_range=f_range, stats=0, x_off=99, l_off_x=99, tick_size=0,
                          center_y=True, tit_size=axis_title_size, y_off=.7)
        h.Draw('{}y+'.format(draw_opt) if 'TGraph' in h.Class_Name() else 'histy+')

    def draw_title_pad(self):
        self.draw_tpad('p2', transparent=True)
        bias_str = 'at {b} V'.format(b=self.Bias) if self.Bias else ''
        run_str = '{n}'.format(n=self.RunNumber) if not self.IsCollection else 'Plan {rp}'.format(rp=self.Analysis.RunPlan)
        text = 'Currents of {dia} {b} - Run {r} - {n}'.format(dia=self.DiamondName, b=bias_str, r=run_str, n=self.Name)
        self.draw_tlatex(pad_margins[0], 1.02 - pad_margins[-1], text, align=11, size=.06)

    def find_margins(self):
        x = [min(self.Time), max(self.Time)]
        dx = .05 * (x[1] - x[0])
        y = [min(self.Currents), max(self.Currents)]
        dy = .01 * (y[1] - y[0])
        return {'x': [x[0] - dx, x[1] + dx], 'y': [y[0] - dy, y[1] + dy]}

    def set_margins(self):
        self.Margins = self.find_margins()

    def make_graphs(self, averaging=1):
        xv = array(self.Time)
        xc = array(average_list(self.Time, averaging))
        # current
        y = array(average_list(self.Currents, averaging))
        g1 = TGraph(len(xc), xc, y)
        self.format_histo(g1, 'Current', '', color=col_cur, markersize=.5)
        g1.SetTitle('')
        # voltage
        y = array(self.Voltages)
        g2 = TGraph(len(xv), xv, y)
        self.format_histo(g2, 'Voltage', '', color=col_vol, markersize=.5)
        self.CurrentGraph = g1
        self.VoltageGraph = g2

    def draw_time_axis(self, y, opt=''):
        x = self.Margins['x']
        a1 = self.draw_x_axis(y, x[0], x[1], 'Time [hh:mm]    ', off=1.2, tit_size=.05, opt=opt, lab_size=.05, tick_size=.3, l_off=.01)
        a1.SetTimeFormat("%H:%M")
        a1.SetTimeOffset(-3600)

    # endregion

    def run_exists(self, run):
        if run in self.RunLogs:
            return True
        else:
            log_warning('Run {run} does not exist in {tc}!'.format(run=run, tc=self.print_testcampaign(pr=False)))
            return False

    def print_run_times(self, run):
        run = str(run)
        if not self.run_exists(run):
            return
        log = self.RunLogs[run]
        out = '{date}: {start}-{stop}'.format(date=log['begin date'], start=log['start time'], stop=log['stop time'])
        print out


if __name__ == "__main__":
    pars = ArgumentParser()
    pars.add_argument('start', nargs='?', default='01/01-00:00:00')
    pars.add_argument('stop', nargs='?', default=None)
    pars.add_argument('-d', '--dia', nargs='?', default='bcmprime-c1')
    pars.add_argument('-tc', '--testcampaign', nargs='?', default='')
    pars.add_argument('-v', '--verbose', action='store_false')
    pars.add_argument('-rp', '--runplan', nargs='?', default=None)
    args = pars.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    a = Elementary(tc)
    print_banner('STARTING CURRENT TOOL')
    a.print_testcampaign()
    start, end = args.start, args.stop
    sel = RunSelection(testcampaign=tc)
    if args.runplan is not None:
        sel.select_runs_from_runplan(args.runplan)
        start = str(sel.get_selected_runs()[0])
        end = str(sel.get_selected_runs()[-1])
    else:
        sel.select_diamond_runs(args.dia)
        if args.start.isdigit():
            sel.unselect_all_runs()
            sel.select_runs_in_range(int(args.start), int(args.stop))
    z = Currents(sel, start_run=start, verbose=args.verbose)
    z.set_start_stop(start, end)
