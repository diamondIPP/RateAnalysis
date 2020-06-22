from __future__ import print_function
from glob import glob

from analysis import Analysis, format_histo, join, basename
from run_selection import RunSelection
from ROOT import TGraph, TProfile, TH1F, TH2F
from numpy import genfromtxt, isnan, datetime64, vectorize, invert, zeros

from utils import *

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


class Currents(Analysis):
    """reads in information from the keithley log file"""

    def __init__(self, analysis=None, test_campaign=None, dut=None, begin=None, end=None, averaging=None, verbose=False):
        Analysis.__init__(self, test_campaign if analysis is None else analysis.TCString, verbose=verbose)

        # Settings
        self.Averaging = averaging
        self.TimeZone = timezone('Europe/Zurich')

        # Config
        self.Analysis = analysis
        self.IsCollection = hasattr(analysis, 'Runs')
        self.RunSelection = RunSelection(testcampaign=self.TCString)
        self.RunLogs = self.RunSelection.RunInfos
        self.Run = self.RunSelection.Run
        self.RunPlan = self.load_run_plan()  # required for plotting
        self.HVConfig = self.load_parser()
        self.set_save_directory('currents')
        self.Bias = self.load_bias()

        # Times
        self.Begin, self.End = self.load_times(begin, end)

        # DUT
        self.DUT = self.get_dut(dut)

        # HV Device Info
        self.Number = self.get_device_number()
        self.Channel = self.get_device_channel()
        self.Name = self.HVConfig.get('HV{}'.format(self.Number), 'name')
        self.Brand = remove_digits(self.Name.split('-')[0])
        self.Model = self.HVConfig.get('HV{}'.format(self.Number), 'model')
        self.Precision = .005 if '237' in self.Name else .05

        self.DataPath = self.find_data_path()

        # data
        self.LogNames = None
        self.IgnoreJumps = True
        self.Currents = []
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

    # ----------------------------------------
    # region INIT
    def load_bias(self):
        return self.Analysis.Bias if hasattr(self.Analysis, 'Bias') else None

    def load_run_number(self):
        return None if self.Analysis is None else self.Run.Number if not self.IsCollection else self.Analysis.RunPlan

    def load_run_plan(self):
        return self.RunSelection.SelectedRunplan if self.Analysis is None else self.Analysis.RunPlan if self.IsCollection else None

    def load_parser(self):
        parser = ConfigParser()
        file_path = self.MainConfig.get('BASIC', 'hvconfigfile') if self.MainConfig.has_option('BASIC', 'hvconfigfile') else join(self.Run.TCDir, 'HV.cfg')
        if not file_exists(file_path):
            critical('HV info file "{f}" does not exist'.format(f=file_path))
        parser.read(file_path)
        self.info('HV Devices: {}'.format(', '.join(name for name in parser.sections() if name.startswith('HV'))))
        return parser

    def load_times(self, begin, end):
        if self.Analysis is None:
            if str(begin).isdigit():  # run number or run plan is provided
                self.RunSelection.select_runs_in_range(begin, end if end is not None else begin) if end or end is None else self.RunSelection.select_runs_from_runplan(begin)
                return self.RunSelection.get_start_time(), self.RunSelection.get_end_time()
            else:  # actual time strings are provided
                return (self.TimeZone.localize(datetime.strptime('{}-{}'.format(self.TestCampaign.year, t), '%Y-%m/%d-%H:%M:%S')) for t in [begin, end])
        return self.load_ana_start_time(), self.load_ana_end_time()

    def get_dut(self, number):
        if self.Analysis is not None:
            return self.Analysis.DUT
        elif self.RunSelection.has_selected_runs():
            return self.RunSelection.SelectedDUT
        from dut_analysis import DUT
        return DUT(number, next(log for log in self.RunLogs.itervalues() if conv_log_time(log['starttime0']) > self.Begin))

    def get_device_str(self):
        if self.Analysis is not None:
            run_info = self.Analysis.FirstAnalysis.Run.RunInfo if self.IsCollection else self.Analysis.Run.RunInfo
        elif self.RunSelection.has_selected_runs():
            run_info = self.RunLogs[str(self.RunSelection.get_selected_runs()[0])]
        else:
            run_info = next(log for log in self.RunLogs.itervalues() if conv_log_time(log['starttime0']) > self.Begin)
        return str(run_info['dia{}supply'.format(self.DUT.Number)])

    def get_device_number(self):
        return self.get_device_str().split('-')[0]

    def get_device_channel(self):
        words = self.get_device_str().split('-')
        return words[1] if len(words) > 1 else '0'

    def find_data_path(self):
        hv_dir = self.MainConfig.get('BASIC', 'hvdatapath') if self.MainConfig.has_option('BASIC', 'hvdatapath') else join(self.Run.TCDir, 'HVClient')
        data_dir = join(hv_dir, '{}_CH{}'.format(self.Name, self.Channel))
        if not dir_exists(data_dir):
            critical('HV data path "{}" does not exist!'.format(data_dir))
        return data_dir

    def load_time(self, t, t_log):
        t = self.TimeZone.localize(datetime.fromtimestamp(t)) if t is not None else t_log
        return t_log if t.year < 2000 or t.day != t_log.day else t

    def load_ana_start_time(self):
        ana = self.Analysis if not self.IsCollection else self.Analysis.FirstAnalysis
        return self.load_time(ana.Run.StartTime if hasattr(ana.Run, 'StartTime') else None, ana.Run.LogStart)

    def load_ana_end_time(self):
        ana = self.Analysis if not self.IsCollection else self.Analysis.LastAnalysis
        return self.load_time(ana.Run.EndTime if hasattr(ana.Run, 'EndTime') else None, ana.Run.LogEnd)

    def reset_data(self):
        self.Currents = []
        self.Voltages = []
        self.Time = []
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region DATA ACQUISITION
    def get_log_names(self):
        log_names = sorted([name for name in glob(join(self.DataPath, '*'))])
        start_log = next((i for i, log_name in enumerate(log_names) if self.get_log_date(log_name) >= self.Begin), 1) - 1
        end_log = next((i for i, log_name in enumerate(log_names) if self.get_log_date(log_name) >= self.End), -1)
        return log_names[start_log:end_log]

    def get_log_date(self, name):
        log_date = ''.join(basename(name).split('_')[-6:])
        return self.TimeZone.localize(datetime.strptime(log_date, '%Y%m%d%H%M%S.log'))

    def find_data(self):
        if len(self.Currents) > 0:
            return
        f = vectorize(time_stamp)
        for i, log_file_name in enumerate(self.get_log_names()):
            data = genfromtxt(log_file_name, usecols=arange(3), dtype=[object, float, float])
            data = data[where(invert(isnan(data['f1'])))[0]]  # remove text entries
            log_date = self.get_log_date(log_file_name)
            data['f0'] = f((log_date.strftime('%Y-%m-%d ') + data['f0']).astype(datetime64).astype(datetime), log_date.utcoffset().seconds)
            data = data[where((data['f0'] >= time_stamp(self.Begin, off=True)) & (data['f0'] <= time_stamp(self.End, off=True)))]
            data = data[where(abs(data['f2']) < 1e-3)]  # filter our very high unphysical currents > 3mA
            if self.IgnoreJumps:  # filter out jumps
                data = data[where(abs(data['f2'][:-1]) * 100 > abs(data['f2'][1:]))[0] + 1]
            self.Time = concatenate([self.Time, data['f0'].astype('i4') + self.Begin.utcoffset().seconds])  # in local start time
            self.Voltages = concatenate([self.Voltages, data['f1']])
            self.Currents = concatenate([self.Currents, data['f2'] * 1e9])  # unit nA
        if not self.Currents.size:
            self.Time = array([time_stamp(self.Begin), time_stamp(self.End)])
            self.Currents = zeros(2)
            self.Voltages = zeros(2)
        if mean(self.Currents) < 0:
            self.Currents *= -1
        if self.Analysis is not None:
            self.Time -= self.Time[0] - self.Analysis.StartTime  # synchronise time vectors
    # endregion DATA ACQUISITION
    # ----------------------------------------

    # ----------------------------------------
    # region PLOTTING
    def get_x_bins(self, bw=5):
        bins = arange(self.Time[0], self.Time[-1], bw) if self.Analysis is None else self.Analysis.get_raw_time_bins(bw)[1] if self.IsCollection else self.Analysis.Bins.get_raw_time(bw)[1]
        return [bins.size - 1, bins]

    def draw_profile(self, bin_width=5, show=True):
        self.find_data()
        p = TProfile('hpr', 'Leakage Current', *self.get_x_bins(bin_width))
        for t, c in zip(self.Time, self.Currents):
            p.Fill(t, c)
        format_histo(p, x_tit='Time [hh:mm]', y_tit='Current [nA]', y_off=.8, markersize=.7, stats=0, t_ax_off=0)
        self.draw_histo(p, '', show, lm=.08, draw_opt='p', x=1.5, y=.75)
        return p

    def draw_flux_correlation(self, bin_fac=1, show=True):
        p1 = self.draw_profile(show=False)
        p2 = self.Analysis.draw_flux(show=False)
        ybins = log_bins(int(sqrt(p1.GetNbinsX()) * bin_fac) * 3, .1, p1.GetMaximum() * 2)
        xbins = log_bins(int(sqrt(p1.GetNbinsX()) * bin_fac), .1, p2.GetMaximum() * 2)
        h = TH2F('gfcc', 'Correlation of Flux and Current', *(xbins + ybins))
        for i in xrange(p1.GetNbinsX()):
            if p1.GetBinContent(i) and p2.GetBinContent(i):
                h.Fill(p2.GetBinContent(i), p1.GetBinContent(i))
        format_histo(h, y_tit='Current [nA]', x_tit='Flux [kHz/cm^{2}', stats=0, y_off=1.2)
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
        xmin, xmax = m - 3 * s, m + 3 * s
        max_bins = int((xmax - xmin) * 100 if '237' in self.Name else 10)
        h = TH1F('hcd', 'Current Distribution', min(5 * int(sqrt(len(self.Currents))), max_bins), xmin, xmax)
        for current in self.Currents:
            h.Fill(current)
        format_histo(h, x_tit='Current [nA]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.draw_histo(h, '', show, lm=.13)
        return h

    def get_current(self):
        if self.Analysis is not None and not self.Analysis.DUT.Bias:
            warning('Bias of run {} is 0!'.format(self.Run.Number))
            current = make_ufloat((0, 0))
        else:
            h = self.draw_distribution(show=False)
            if h.GetEntries() < 3:
                return None
            values = [h.GetBinCenter(i) for i in xrange(h.GetNbinsX())]
            weights = [h.GetBinContent(i) for i in xrange(h.GetNbinsX())]
            m, s = mean_sigma(values, weights)
            fit = h.Fit('gaus', 'sq0', '', m - s, m + s)
            fm, fs = fit.Parameter(1), fit.Parameter(2)
            if .8 * m < fit.Parameter(1) < 1.2 * m and s > 0 and fs < fm and fit.ParError(1) < m:  # only use gauss fit if its not deviating too much from the the mean
                current = ufloat(fm, fs + self.Precision + .03 * fm)  # add .05 as uncertainty of the device and 5% systematic error
            else:
                current = ufloat(h.GetMean(), h.GetMeanError() + .05 + .05 * h.GetMean())
        self.Analysis.server_pickle(self.make_pickle_path('Currents', run=self.Run.Number, ch=self.DUT.Number), current)
        return current

    def draw_iv(self, show=True):
        self.find_data()
        x = [v for v, c in zip(self.Voltages, self.Currents) if c != 590]
        y = [c for c in self.Currents if c != 590]
        g = self.make_tgrapherrors('giv', 'I-V Curve for {}'.format(self.Analysis.DUT.Name), x=x, y=y)
        format_histo(g, x_tit='Voltage [V]', y_tit='Current [nA]', y_off=1.4)
        self.draw_histo(g, draw_opt='ap', logy=True, lm=.12, show=show)
        return g

    def draw_indep_graphs(self, rel_time=False, ignore_jumps=True, v_range=None, f_range=None, c_range=None, averaging=1, with_flux=False, draw_opt='ap', show=True):
        self.IgnoreJumps = ignore_jumps
        self.set_graphs(averaging)
        c = self.make_canvas('cc', 'Keithley Currents for Run {0}'.format(self.Run.Number), x=1.5, y=.75, show=show)
        self.draw_flux_pad(f_range, rel_time, draw_opt) if with_flux else self.draw_voltage_pad(v_range, draw_opt)
        self.draw_title_pad()
        self.draw_current_pad(rel_time, c_range, draw_opt)
        if self.IsCollection:
            self.draw_irradiation(make_irr_string(self.Analysis.RunSelection.get_irradiation()))
        self.Stuff.append(c)
        run = self.Analysis.RunPlan if self.IsCollection else self.Run.Number
        save_name = 'Currents{}_{}_{}'.format(self.TCString, run, self.DUT.Number)
        self.save_canvas(c, name=save_name, sub_dir='currents', show=show, ftype='png')

    def zoom_pads(self, low, high):
        self.VoltageGraph.GetXaxis().SetRangeUser(low, high)
        self.CurrentGraph.GetXaxis().SetRangeUser(low, high)

    def draw_current_pad(self, rel_t, c_range, draw_opt):
        self.draw_tpad('p3', gridx=True, margins=pad_margins, transparent=True)
        g = self.CurrentGraph
        format_histo(g, x_tit='#font[22]{Time [hh:mm]}', lab_size=label_size, x_off=1.05, tit_size=axis_title_size, t_ax_off=self.Time[0] if rel_t else 0, y_off=.55, yax_col=col_cur,
                     y_tit='#font[22]{Current [nA]}', center_y=True, x_range=[self.Time[0], self.Time[-1]], y_range=c_range, color=col_cur, markersize=marker_size)
        self.CurrentGraph.Draw(draw_opt)

    def draw_voltage_pad(self, v_range, draw_opt='ap'):
        self.draw_tpad('p1', gridy=True, margins=pad_margins, transparent=True)
        g = self.VoltageGraph
        v_range = [-1100, 1100] if v_range is None else v_range
        format_histo(g, y_range=v_range, y_tit='#font[22]{Voltage [V]}', x_range=[self.Time[0], self.Time[-1]], tit_size=axis_title_size, tick_size=0, x_off=99, l_off_x=99, center_y=True,
                     color=col_vol, y_off=title_offset, markersize=marker_size, yax_col=col_vol, lw=3, lab_size=label_size)
        g.Draw('{}y+'.format(draw_opt))

    def draw_flux_pad(self, f_range, rel_t=False, draw_opt='ap'):
        pad = self.draw_tpad('p1', margins=pad_margins, transparent=True, logy=True)
        h = self.Analysis.draw_flux(rel_time=rel_t, show=False)
        pad.cd()
        f_range = [1, 20000] if f_range is None else f_range
        format_histo(h, title=' ', y_tit='#font[22]{Flux [kHz/cm^{2}]}', fill_color=4000, fill_style=4000, lw=3, y_range=f_range, stats=0, x_off=99, l_off_x=99, tick_size=0,
                     center_y=True, tit_size=axis_title_size, y_off=.7)
        h.Draw('{}y+'.format(draw_opt) if 'TGraph' in h.Class_Name() else 'histy+')

    def draw_title_pad(self):
        self.draw_tpad('p2', transparent=True)
        bias_str = 'at {b} V'.format(b=self.Bias) if self.Bias else ''
        run_str = '{n}'.format(n=self.Run.Number) if not self.IsCollection else 'Plan {rp}'.format(rp=self.Analysis.RunPlan)
        text = 'Currents of {dia} {b} - Run {r} - {n}'.format(dia=self.DUT.Name, b=bias_str, r=run_str, n=self.Name)
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
        format_histo(g1, 'Current', '', color=col_cur, markersize=.5)
        g1.SetTitle('')
        # voltage
        y = array(self.Voltages)
        g2 = TGraph(len(xv), xv, y)
        format_histo(g2, 'Voltage', '', color=col_vol, markersize=.5)
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
        out = '{date}: {start}-{stop}'.format(date=log['begin date'], start=log['begin time'], stop=log['stop time'])
        print(out)

    def get_time_from_log(self, t_str, year_str):
        return self.TimeZone.localize(datetime.strptime(year_str.strftime('%Y%m%d') + t_str, '%Y%m%d%H:%M:%S'))


if __name__ == '__main__':

    aparser = ArgumentParser()
    aparser.add_argument('dut', nargs='?', default=1, type=int, help='dut number [default: 1] (choose from 1,2,...)')
    aparser.add_argument('begin', nargs='?', default=12)
    aparser.add_argument('end', nargs='?', default=None)
    aparser.add_argument('-tc', '--testcampaign', nargs='?', default=None, help='YYYYMM beam test [default in main.ini]')
    aparser.add_argument('-v', '--verbose', action='store_false')
    aparser.add_argument('-c', '--collection', action='store_true', help='begin analysis collection')
    pargs = aparser.parse_args()

    z = Currents(test_campaign=pargs.testcampaign, dut=pargs.dut, begin=pargs.begin, end=pargs.end if not pargs.collection else False, verbose=pargs.verbose)
