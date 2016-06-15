# ====================================
# IMPORTS
# ====================================
from datetime import datetime, timedelta
from glob import glob
from time import sleep
from json import load
from collections import OrderedDict

from ROOT import TCanvas, TPad, TText, TGraph, kCanDelete
from ConfigParser import ConfigParser
from numpy import array, mean
from argparse import ArgumentParser

from Elementary import Elementary

# ====================================
# CONSTANTS
# ====================================
axis_title_size = 0.04
label_size = .03
col_vol = 602  # 807
col_cur = 418


# ====================================
# CLASS FOR THE DATA
# ====================================
class Currents(Elementary):
    """reads in information from the keithley log file"""

    def __init__(self, analysis, averaging=False, points=10, hv_nr=None, verbose=False):
        self.analysis = analysis
        self.IsCollection = hasattr(analysis, 'runs')
        Elementary.__init__(self, verbose=verbose)

        self.DoAveraging = averaging
        self.Points = points

        # config
        self.ConfigParser = self.load_parser()
        if self.ConfigParser is None:
            return

        # analysis/run info
        self.TimeOffset = self.run_config_parser.getint('BASIC', 'hvtimeoffset')
        self.RunNumber = self.load_run_number()
        self.RunLogs = self.load_runlogs()
        if analysis is not None:
            self.RunInfo = analysis.run.RunInfo if not self.IsCollection else analysis.get_first_analysis().RunInfo
            self.Channel = analysis.channel
            self.DiamondName = analysis.diamond_name
        self.StartTime = self.load_start_time()
        self.StopTime = self.load_stop_time()

        # device info
        self.Number = self.get_device_nr(hv_nr)
        self.Channel = self.get_device_channel()
        self.Brand = self.ConfigParser.get('HV' + self.Number, 'name').split('-')[0].strip('0123456789')
        self.Model = self.ConfigParser.get('HV' + self.Number, 'model')
        self.Name = '{0} {1}'.format(self.Brand, self.Model)

        self.DataPath = self.find_data_path()
        self.OldDataPath = self.find_data_path(old=True)
        self.LogNames = None

        # data
        self.Currents = []
        self.Voltages = []
        self.Time = []
        self.MeanCurrent = 0
        self.MeanVoltage = 0

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

    def __del__(self):
        print 'Cleaning up'
        for pad in [self.TitlePad, self.VoltagePad, self.CurrentPad]:
            try:
                print pad
            except Exception as err:
                self.log_warning(err)


    # ==========================================================================
    # region INIT
    def load_run_number(self):
        nr = None
        if self.analysis is not None:
            nr = self.analysis.run_number if not self.IsCollection else self.analysis.run_plan
        return nr

    def load_runlogs(self):
        try:
            f = open(self.run_config_parser.get('BASIC', 'runinfofile'))
            data = load(f)
            f.close()
        except IOError as err:
            self.log_warning('{err}\nCould not load default RunInfo!'.format(err))
            return None
        RunLogs = OrderedDict(sorted(data.iteritems()))
        return RunLogs

    def load_run_config(self):
        nr = None
        if self.analysis is not None:
            nr = self.analysis.run_number if not self.IsCollection else self.analysis.get_first_analysis().run_number
        return self.load_run_configs(nr)

    def load_parser(self):
        parser = ConfigParser()
        if not self.run_config_parser.has_option('BASIC', 'hvconfigfile'):
            self.log_warning('Missing hv info in RunConfig file')
            return None
        parser.read(self.run_config_parser.get('BASIC', 'hvconfigfile'))
        self.log_info('HV Devices: {0}'.format([name for name in parser.sections() if name.startswith('HV')]))
        return parser

    def get_device_nr(self, hv_nr):
        if self.analysis is None:
            return hv_nr
        full_str = self.RunInfo['dia{dia}supply'.format(dia=1 if not self.Channel else 2)]
        return str(full_str.split('-')[0])

    def get_device_channel(self):
        if self.analysis is None:
            return '0'
        full_str = self.RunInfo['dia{dia}supply'.format(dia=1 if not self.Channel else 2)]
        return full_str.split('-')[1] if len(full_str) > 1 else '0'

    def find_data_path(self, old=False):
        hv_datapath = self.run_config_parser.get('BASIC', 'hvdatapath')
        string = '{data}{dev}_CH{ch}/' if not old else '{data}{dev}/'
        return string.format(data=hv_datapath, dev=self.ConfigParser.get('HV' + self.Number, 'name'), ch=self.Channel)

    def load_start_time(self):
        if self.analysis is not None:
            return self.analysis.get_first_analysis().run.log_start + timedelta(hours=self.TimeOffset) if self.IsCollection else self.analysis.run.log_start + timedelta(hours=self.TimeOffset)
        else:
            return None

    def load_stop_time(self):
        if self.analysis is not None:
            return self.analysis.get_last_analysis().run.log_stop + timedelta(hours=self.TimeOffset) if self.IsCollection else self.analysis.run.log_stop + timedelta(hours=self.TimeOffset)
        else:
            return None

    def set_start_stop(self, start, stop):
        self.StartTime = datetime.strptime(start, '%Y/%m/%d-%H:%M:%S')
        self.StopTime = datetime.strptime(stop, '%Y/%m/%d-%H:%M:%S')
    # endregion

    # ==========================================================================
    # region ACQUIRE DATA
    def get_logs_from_start(self):
        log_names = sorted([name for name in glob(self.DataPath + '*')] + [name for name in glob(self.OldDataPath + '*')])
        start_log = None
        for i, name in enumerate(log_names):
            log_date = self.get_log_date(name)
            if log_date >= self.StartTime:
                break
            start_log = i
        self.log_info('Starting with log: {0}'.format(log_names[start_log].split('/')[-1]))
        return log_names[start_log:]

    @staticmethod
    def get_log_date(name):
        log_date = name.split('/')[-1].split('_')
        log_date = ''.join(log_date[-6:])
        return datetime.strptime(log_date, '%Y%m%d%H%M%S.log')

    def set_start(self, zero=False):
        self.Currents.append(self.Currents[-1] if not zero else 0)
        self.Voltages.append(self.Voltages[-1] if not zero else 0)
        self.Time.append((self.StartTime - datetime(self.StartTime.year, 1, 1)).total_seconds())

    def set_stop(self, zero=False):
        self.Currents.append(self.Currents[-1] if not zero else 0)
        self.Voltages.append(self.Voltages[-1] if not zero else 0)
        self.Time.append((self.StopTime - datetime(self.StopTime.year, 1, 1)).total_seconds())

    def find_data(self):
        if self.Currents:
            return
        stop = False
        self.LogNames = self.get_logs_from_start()
        for i, name in enumerate(self.LogNames):
            self.MeanCurrent = 0
            self.MeanVoltage = 0
            log_date = self.get_log_date(name)
            data = open(name, 'r')
            # jump to the correct line of the first file
            if not i:
                self.find_start(data, log_date)
            index = 0
            if index == 1:
                self.set_start()
            for line in data:
                # if index < 20:
                #     print line
                info = line.split()
                if self.isfloat(info[1]) and len(info) > 2:
                    now = datetime.strptime(log_date.strftime('%Y%m%d') + info[0], '%Y%m%d%H:%M:%S')
                    if self.StartTime < now < self.StopTime and float(info[2]) < 1e30:
                        self.save_data(now, info, index)
                        index += 1
                    if self.StopTime < now:
                        stop = True
                        break
            data.close()
            if stop:
                break
        if self.Currents:
            self.set_stop()
        if not self.Currents:
            self.set_start(zero=True)
            self.set_stop(zero=True)

    def save_data(self, now, info, index, shifting=False):
        total_seconds = (now - datetime(now.year, 1, 1)).total_seconds()
        if self.StartTime < now < self.StopTime and float(info[2]) < 1e30:
            index += 1
            if self.DoAveraging:
                if not shifting:
                    self.MeanCurrent += float(info[2]) * 1e9
                    self.MeanVoltage += float(info[1])
                    if index % self.Points == 0:
                        if mean(self.Currents) < 5 * self.MeanCurrent / self.Points:
                            self.Currents.append(self.MeanCurrent / self.Points)
                            self.Time.append(total_seconds)
                            self.Voltages.append(self.MeanVoltage / self.Points)
                            self.MeanCurrent = 0
                            self.MeanVoltage = 0
                # else:
                #     if index <= self.Points:
                #         self.mean_curr += float(info[2]) * 1e9
                #         dicts[1][key].append(self.mean_curr / index)
                #         if index == self.Points:
                #             self.mean_curr /= self.Points
                #     else:
                #         mean_curr = self.mean_curr * weight + (1 - weight) * float(info[2]) * 1e9
                #         dicts[1][key].append(mean_curr)
                #     dicts[0][key].append(convert_time(now))
                #     dicts[2][key].append(float(info[1]))
            else:
                if len(self.Currents) > 100 and abs(self.Currents[-1] * 100) < abs(float(info[2]) * 1e9):
                    if abs(self.Currents[-1]) > 0.01:
                        return
                self.Currents.append(float(info[2]) * 1e9)
                self.Time.append(total_seconds)
                self.Voltages.append(float(info[1]))

    def find_start(self, data, log_date):
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
                if self.isfloat(info[1]):
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
    def set_graphs(self, rel_time=True):
        self.find_data()
        self.convert_to_relative_time() if rel_time else self.do_nothing()
        sleep(.1)
        self.make_graphs()
        self.set_margins()

    def draw_graphs(self, relative_time=False):
        if not self.Currents:
            self.find_data()
        if self.Margins is None:
            self.set_margins()
        if relative_time:
            self.convert_to_relative_time()
            sleep(.1)
            self.make_graphs()
            self.set_margins()
        if self.CurrentGraph is None:
            self.make_graphs()
        c = TCanvas('c', 'Keithley Currents for Run {0}'.format(self.RunNumber), 1000, 1000)
        self.make_pads()
        self.draw_pads()

        self.VoltagePad.cd()
        self.draw_voltage_frame()
        self.VoltageGraph.Draw('p')
        a = self.draw_voltage_axis()
        a.Draw()

        self.TitlePad.cd()
        t = self.draw_title_pad()
        t.Draw()

        self.CurrentPad.cd()
        self.draw_current_frame()
        self.CurrentGraph.Draw('pl')

        c.Update()

        self.Histos[0] = [c, t, a]
        self.save_plots('Currents', sub_dir=self.analysis.save_dir)

    def draw_indep_graphs(self, rel_time=False):
        if not self.Currents:
            self.set_graphs(rel_time)
        c = TCanvas('c', 'Keithley Currents for Run {0}'.format(self.RunNumber), 1500, 750)
        self.make_pads()
        self.draw_pads()

        self.draw_voltage_pad()
        self.draw_title_pad()
        self.draw_current_pad()

        self.Stuff.append(c)
        self.save_plots('Currents', sub_dir='Currents')

    def draw_current_pad(self):
        self.CurrentPad.cd()
        self.draw_current_frame()
        self.CurrentGraph.Draw('pl')

    def draw_voltage_pad(self):
        self.VoltagePad.cd()
        self.draw_voltage_frame()
        self.VoltageGraph.Draw('p')
        self.draw_voltage_axis()

    def make_pads(self):
        self.VoltagePad = self.make_tpad('p1', gridy=True, margins=[.08, .07, .15, .15])
        self.CurrentPad = self.make_tpad('p2', gridx=True, margins=[.08, .07, .15, .15], transparent=True)
        self.TitlePad = self.make_tpad('p3', transparent=True)

    def draw_pads(self):
        for p in [self.VoltagePad, self.TitlePad, self.CurrentPad]:
            p.Draw()

    def draw_title_pad(self):
        self.TitlePad.cd()
        text = 'Currents measured by {0}'.format(self.Name)
        t1 = TText(0.08, 0.88, text)
        t1.SetTextSize(0.05)
        t1.Draw()
        self.Stuff.append(t1)

    def find_margins(self):
        x = [min(self.Time), max(self.Time)]
        dx = .05 * (x[1] - x[0])
        y = [min(self.Currents), max(self.Currents)]
        dy = .01 * (y[1] - y[0])
        return {'x': [x[0] - dx, x[1] + dx], 'y': [y[0] - dy, y[1] + dy]}

    def set_margins(self):
        self.Margins = self.find_margins()

    def make_graphs(self):
        tit = ' measured by {0}'.format(self.Name)
        x = array(self.Time)
        # current
        y = array(self.Currents)
        g1 = TGraph(len(x), x, y)
        self.format_histo(g1, 'Current', 'Current' + tit, color=col_cur, markersize=.5)
        # voltage
        y = array(self.Voltages)
        g2 = TGraph(len(x), x, y)
        self.format_histo(g2, 'Voltage', 'Voltage' + tit, color=col_vol, markersize=.5)
        self.CurrentGraph = g1
        self.VoltageGraph = g2

    def make_tpad(self, name, tit='', pos=None, fill_col=0, gridx=False, gridy=False, margins=None, transparent=False):
        margins = [.1, .1, .1, .1] if margins is None else margins
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(name, tit, *pos)
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

    def draw_voltage_axis(self):
        a1 = self.make_tgaxis(self.Margins['x'][1], -1100, 1100, '#font[22]{Voltage [V]}', color=col_vol, offset=.6, tit_size=.05, line=False, opt='+L', width=2)
        a1.CenterTitle()
        a1.SetLabelSize(label_size)
        a1.SetLabelColor(col_vol)
        a1.SetLabelOffset(0.01)
        a1.Draw()
        self.Stuff.append(a1)

    def draw_current_axis(self):
        y = self.Margins['y']
        yadd = (y[1] - y[0]) * .1
        diff = (y[1] - y[0])
        a1 = self.make_tgaxis(self.Margins['x'][1], y[0] - yadd, y[1] + yadd + diff, 'Current [nA]', offset=.7, tit_size=.04, line=False, opt='+L', width=2)
        a1.SetLabelSize(label_size)
        a1.SetLabelOffset(0.01)
        a1.Draw()
        self.Stuff.append(a1)

    def draw_pulser_axis(self, ymin, ymax):
        x = self.Margins['x'][0] - (self.Margins['x'][1] - self.Margins['x'][0]) * .07
        a1 = self.make_tgaxis(x, ymin, ymax, 'pulse height [au]', offset=.8, tit_size=.04, line=False, opt='-R', width=2)
        a1.SetLabelColor(859)
        a1.SetLineColor(859)
        a1.SetLabelSize(label_size)
        a1.SetLabelOffset(0.01)
        a1.Draw()
        self.Stuff.append(a1)

    def draw_current_frame(self):
        m = self.Margins
        h2 = self.CurrentPad.DrawFrame(m['x'][0], m['y'][0], m['x'][1], m['y'][1])
        # X-axis
        h2.GetXaxis().SetTitle("#font[22]{time [hh:mm]}")
        h2.GetXaxis().SetTimeFormat("%H:%M")
        h2.GetXaxis().SetTimeOffset(-3600)
        h2.GetXaxis().SetTimeDisplay(1)
        h2.GetXaxis().SetLabelSize(label_size)
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
        h2.GetYaxis().SetLabelSize(label_size)
        h2.GetYaxis().SetTitleSize(axis_title_size)
        h2.GetYaxis().SetTitleOffset(.6)
        self.Stuff.append(h2)

    def draw_ph_frame(self, pad, margins):
        pad.cd()
        m = self.Margins
        h2 = pad.DrawFrame(m['x'][0], margins[0], m['x'][1], margins[1])
        # X-axis
        h2.GetXaxis().SetTitle("time [hh:mm]")
        h2.GetXaxis().SetTimeFormat("%H:%M")
        h2.GetXaxis().SetTimeOffset(-3600)
        h2.GetXaxis().SetTimeDisplay(1)
        h2.GetXaxis().SetLabelSize(label_size)
        h2.GetXaxis().SetTitleSize(axis_title_size)
        h2.GetXaxis().SetTitleOffset(1.05)
        h2.GetXaxis().SetTitleSize(0.04)
        # Y-axis
        h2.GetYaxis().SetTitleOffset(0.6)
        h2.GetYaxis().SetTitleSize(0.04)
        # h2.GetYaxis().SetTitle("Pulse Height [au]")
        h2.GetYaxis().SetLabelSize(label_size)
        h2.GetYaxis().SetTitleSize(0.04)
        h2.GetYaxis().SetTitleOffset(0.75)

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

    def draw_curr_frame(self, pad):
        pad.cd()
        m = self.Margins
        y = m['y']
        yadd = (y[1] - y[0]) * .1
        diff = (m['y'][1] - m['y'][0])
        h1 = pad.DrawFrame(m['x'][0], m['y'][0] - yadd, m['x'][1], m['y'][1] + diff + yadd)
        h1.SetTitleSize(axis_title_size)
        h1.GetXaxis().SetTickLength(0)
        h1.GetYaxis().SetTickLength(0)
        h1.GetXaxis().SetLabelOffset(99)
        h1.GetYaxis().SetLabelOffset(99)
        h1.SetLineColor(0)
        
    def draw_pulser_frame(self, pad, ymin, ymax):
        pad.cd()
        m = self.Margins
        h1 = pad.DrawFrame(m['x'][0], ymin, m['x'][1], ymax)
        self.format_frame(h1)

    @staticmethod
    def format_frame(frame):
        fr = frame
        fr.GetYaxis().SetTitleSize(.06)
        fr.GetYaxis().SetTitleOffset(.6)
        fr.GetYaxis().SetLabelSize(.06)
        fr.SetTitleSize(axis_title_size)
        fr.GetXaxis().SetTickLength(0)
        # fr.GetYaxis().SetTickLength(0)
        fr.GetXaxis().SetLabelOffset(99)
        # fr.GetYaxis().SetLabelOffset(99)
        fr.SetLineColor(0)
        fr.GetXaxis().SetTimeDisplay(1)

    def draw_frame(self, pad, ymin, ymax, tit):
        pad.cd()
        x = self.Margins['x']
        fr = pad.DrawFrame(x[0], ymin, x[1], ymax)
        pad.Modified()
        fr.GetYaxis().SetTitle(tit)
        self.format_frame(fr)

    def draw_time_axis(self, y, opt=''):
        x = self.Margins['x']
        a1 = self.make_tgxaxis(x[0], x[1], y, 'Time [hh:mm]    ', offset=1.2, tit_size=.05, line=False, opt=opt)
        a1.SetLabelSize(.05)
        a1.SetTickSize(.3)
        a1.SetTitle("Time [hh:mm]")
        a1.SetTimeFormat("%H:%M")
        a1.SetTimeOffset(-3600)
        a1.SetLabelOffset(0.01)
        a1.Draw()
        self.Stuff.append(a1)
        
    # endregion

    @staticmethod
    def make_transparent(pad):
        pad.SetFillStyle(4000)
        pad.SetFillColor(0)
        pad.SetFrameFillStyle(4000)


if __name__ == "__main__":
    pars = ArgumentParser()
    pars.add_argument('start', nargs='?', default='')
    pars.add_argument('stop', nargs='?', default='2099')
    pars.add_argument('-hv', nargs='?', default='1')
    pars.add_argument('-tc', '--testcampaign', nargs='?', default='')
    pars.add_argument('-v', '--verbose', nargs='?', default=True, type=bool)
    args = pars.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    a = Elementary(tc)
    a.print_banner('STARTING CURRENT TOOL')
    a.print_testcampaign()
    z = Currents(None, hv_nr=args.hv, verbose=args.verbose)
    z.set_start_stop('{y}/{s}'.format(y=a.TESTCAMPAIGN[:4], s=args.start), '{y}/{e}'.format(y=a.TESTCAMPAIGN[:4], e=args.stop))
