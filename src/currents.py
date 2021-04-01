from glob import glob
from os.path import getsize
from ROOT import TH2F
from numpy import genfromtxt, isnan, datetime64, invert, uint32, char
from helpers.save_plots import *
from src.analysis import Analysis
from src.run_selection import RunSelector


class Currents(Analysis):
    """reads in information from the keithley log file"""

    VCol = 602  # 807
    CCol = 899  # 418
    MS = .3

    def __init__(self, analysis=None, test_campaign=None, dut=None, begin=None, end=None, averaging=None, verbose=False):
        Analysis.__init__(self, test_campaign if analysis is None else analysis.TCString, verbose=verbose, sub_dir='currents')

        # Settings
        self.Averaging = averaging
        self.TimeZone = timezone('Europe/Zurich')
        self.DataDir = join(self.TCDir, 'hv')

        # Config
        self.Ana = analysis
        self.IsCollection = hasattr(analysis, 'Runs')
        self.Type = self.Ana.Type if analysis is not None and self.IsCollection else 'None'
        self.RunSelection = RunSelector(testcampaign=self.TCString)
        self.RunLogs = self.RunSelection.RunInfos
        self.Run = self.RunSelection.Run if analysis is None else self.Ana.FirstAnalysis.Run if self.IsCollection else self.Ana.Run
        self.RunPlan = self.load_run_plan()  # required for plotting
        self.HVConfig = self.load_parser()
        self.Bias = self.Ana.Bias if hasattr(self.Ana, 'Bias') else None

        # Times
        self.Begin, self.End = self.load_times(begin, end, dut)

        # DUT
        self.DUT = self.init_dut(dut)

        # HV Device Info
        self.Number = self.load_device_number()
        self.Channel = self.load_device_channel()
        self.Name = self.HVConfig.get('HV{}'.format(self.Number), 'name')
        self.Brand = remove_digits(self.Name.split('-')[0])
        self.Model = self.HVConfig.get('HV{}'.format(self.Number), 'model')
        self.Precision = .005 if '237' in self.Name else .05

        # data
        self.IgnoreJumps = True
        self.Data = self.load_data()

    # ----------------------------------------
    # region INIT
    def load_data(self):
        data_file = join(self.DataDir, 'data.hdf5')
        if not file_exists(data_file):
            self.convert_data()
        data = h5py.File(data_file, 'r')['{}_CH{}'.format(self.Name, self.Channel)]
        data = data[(data['timestamps'] >= time_stamp(self.Begin)) & (data['timestamps'] <= time_stamp(self.End))]
        if not data.size:
            return None
        if self.IgnoreJumps:  # filter out jumps
            data = data[where(abs(data['currents'][:-1]) * 100 > abs(data['currents'][1:]))[0] + 1]  # take out the events that are 100 larger than the previous
        data['currents'] *= 1e9 * sign(mean(data['currents']))  # convert to nA and flip sign if current is negative
        if self.Ana is not None and data.size:
            data['timestamps'] -= uint32(data['timestamps'][0] - time_stamp(self.load_ana_start_time()))  # synchronise time vectors
        return data

    def reload_data(self, ignore_jumps):
        if ignore_jumps != self.IgnoreJumps:
            self.IgnoreJumps = ignore_jumps
            self.Data = self.load_data()

    def load_run_plan(self):
        return self.RunSelection.SelectedRunplan if self.Ana is None else self.Ana.RunPlan if self.IsCollection else None

    def load_parser(self):
        file_path = join(self.DataDir, 'config.ini')
        if not file_exists(file_path):
            critical('HV info file "{f}" does not exist'.format(f=file_path))
        return Config(file_path)

    def load_times(self, begin, end, dut=1):
        if self.Ana is None:
            if str(begin).isdigit():  # run number or run plan is provided
                self.RunSelection.select_runs_in_range(begin, end if end is not None else begin, dut) if end or end is None else self.RunSelection.select_runs_from_runplan(begin)
                return self.RunSelection.get_start_time(), self.RunSelection.get_end_time()
            else:  # actual time strings are provided
                return (self.TimeZone.localize(datetime.strptime('{}-{}'.format(self.TestCampaign.year, t), '%Y-%m/%d-%H:%M:%S')) for t in [begin, end])
        return self.load_ana_start_time(), self.load_ana_end_time()

    def load_time(self, t, t_log):
        t = self.TimeZone.localize(datetime.fromtimestamp(t)) if t is not None else t_log
        return t_log if t.year < 2000 or t.day != t_log.day else t

    def load_ana_start_time(self):
        ana = self.Ana if not self.IsCollection else self.Ana.FirstAnalysis
        return self.load_time(ana.Run.StartTime if hasattr(ana.Run, 'StartTime') else None, ana.Run.LogStart)

    def load_ana_end_time(self):
        ana = self.Ana if not self.IsCollection else self.Ana.LastAnalysis
        return self.load_time(ana.Run.EndTime if hasattr(ana.Run, 'EndTime') else None, ana.Run.LogEnd)

    def init_dut(self, number):
        if self.Ana is not None:
            return self.Ana.DUT
        elif self.RunSelection.has_selected_runs():
            return self.RunSelection.SelectedDUT
        from dut import DUT
        return DUT(number, next(log for log in self.RunLogs.values() if conv_log_time(log['starttime0']) > self.Begin))

    def load_device_str(self):
        if self.Ana is not None:
            run_info = self.Ana.FirstAnalysis.Run.Info if self.IsCollection else self.Ana.Run.Info
        elif self.RunSelection.has_selected_runs():
            run_info = self.RunLogs[self.RunSelection.get_selected_runs()[0]]
        else:
            run_info = next(log for log in self.RunLogs.values() if conv_log_time(log['starttime0']) > self.Begin)
        return str(run_info['dia{}supply'.format(self.DUT.Number)])

    def load_device_number(self):
        return self.load_device_str().split('-')[0]

    def load_device_channel(self):
        words = self.load_device_str().split('-')
        return words[1] if len(words) > 1 else '0'
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_runs(self):
        return self.Ana.get_runs()

    def get_fluxes(self, pbar=False):
        return self.Ana.get_fluxes(pbar=pbar)

    def get_analyses(self):
        return self.Ana.get_analyses()
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DATA ACQUISITION
    def get_log_date(self, name):
        log_date = ''.join(basename(name).split('_')[-6:])
        return self.TimeZone.localize(datetime.strptime(log_date, '%Y%m%d%H%M%S.log'))

    def convert_data(self):
        info('converting hv text files to hdf5 ...')
        self.PBar.start(len(glob(join(self.DataDir, '*', '*.log'))))
        f = h5py.File(join(self.DataDir, 'data.hdf5'), 'w')
        for d in glob(join(self.DataDir, '*_*')):
            arrays = []
            for file_name in sorted(glob(join(d, '*.log'))):
                if getsize(file_name) == 0:
                    remove_file(file_name)
                    self.PBar.update()
                    continue
                log_date = self.get_log_date(file_name)
                data = genfromtxt(file_name, usecols=arange(3), dtype=[('timestamps', 'U10'), ('voltages', 'f2'), ('currents', 'f4')])
                data = data[invert(isnan(data['voltages']))]  # remove text entries
                data['timestamps'] = array(log_date.strftime('%Y-%m-%d ') + char.array(data['timestamps'])).astype(datetime64) - datetime64('1970-01-01T00:00:00') - log_date.utcoffset().seconds
                data = data.astype([('timestamps', 'u4'), ('voltages', 'f2'), ('currents', 'f4')])
                arrays.append(data)
                self.PBar.update()
            if len(arrays):
                f.create_dataset(basename(d), data=concatenate(arrays))
    # endregion DATA ACQUISITION
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_current(self):
        if self.Ana is not None and not self.Ana.DUT.Bias:
            warning('Bias of run {} is 0!'.format(self.Run.Number))
            current = ufloat(0, 0)
        else:
            h = self.draw_distribution(show=False)
            if h.GetEntries() < 3:
                return None
            m, s = mean_sigma(*get_hist_vecs(h, err=False), err=False)
            fit = h.Fit('gaus', 'sq0', '', m - s, m + s)
            fm, fs = fit.Parameter(1), fit.Parameter(2)
            if .8 * m < fit.Parameter(1) < 1.2 * m and s > 0 and fs < fm and fit.ParError(1) < m:  # only use gauss fit if its not deviating too much from the the mean
                current = ufloat(fm, fit.ParError(1) + self.Precision / 2)
            else:
                current = ufloat(h.GetMean(), h.GetMeanError() + self.Precision / 2)
        self.Draw.server_pickle(self.make_pickle_path('Currents', run=self.Run.Number, ch=self.DUT.Number), current)
        return current

    def get_title(self):
        bias_str = 'at {b} V'.format(b=self.Bias) if self.Bias else ''
        run_str = '{n}'.format(n=self.Run.Number) if not self.IsCollection else 'Plan {rp}'.format(rp=self.Ana.RunPlan)
        return 'Currents of {dia} {b} - Run {r} - {n}'.format(dia=self.DUT.Name, b=bias_str, r=run_str, n=self.Name)

    def get_first_log(self):
        names = sorted(glob(join(self.DataDir, '{}_CH{}'.format(self.Name, self.Channel), '*.log')))
        for i, name in enumerate(names):
            if self.get_log_date(name) > self.Begin:
                return names[i - 1]

    def get_run_number(self):
        return None if self.Ana is None else 'RP{}'.format(self.Ana.RunPlan) if self.IsCollection else self.Ana.Run.Number

    def get_t0(self):
        return self.Data['timestamps'][0]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region PLOTTING
    def draw(self, rel_time=False, ignore_jumps=True, v_range=None, c_range=None, averaging=1, draw_opt='al', with_flux=False, f_range=None, show=True):
        self.reload_data(ignore_jumps)
        t, c, v = (average_list(self.Data[n], averaging) for n in ['timestamps', 'currents', 'voltages'])
        gv = self.Draw.graph(t, v, title=self.get_title(), y_tit='Voltage [nA]', yax_col=self.VCol, color=self.VCol, y_range=choose(v_range, [-1100, 1100]), l_off_x=10, x_ticks=0, show=False, lw=2)
        if with_flux:
            gv = self.Ana.draw_flux(rel_time=rel_time, show=False)
            format_histo(gv, title=self.get_title(), fill_color=4000, fill_style=4000, lw=3, y_range=choose(f_range, [5, 20000]), y_tit='Flux [kHz/cm^{2}]')
        gc = Draw.make_tgrapherrors(t, c)
        format_histo(gc, x_tit='Time [hh:mm]', y_tit='Current [nA]', yax_col=self.CCol, color=self.CCol, y_range=choose(c_range, [round_down_to(min(c)), round_up_to(max(c))]))
        x_range = [gv.GetBinLowEdge(1), gv.GetBinLowEdge(gv.GetNbinsX() + 1)] if with_flux else [t[0], t[-1]]
        for g in [gc, gv]:
            format_histo(g, lab_size=.05, x_off=1.05, tit_size=.06, t_ax_off=t[0] if rel_time else 0, y_off=.8, center_y=True, x_range=x_range, markersize=self.MS, stats=0)
        m = [.09, .09, .2, .1]
        self.Draw(gv, show=show, m=m, w=1.5, h=.75, draw_opt='{}y+'.format('hist' if with_flux else draw_opt), logy=with_flux)
        Draw.tpad('pc', transparent=True, margins=m)
        gc.Draw(draw_opt)
        save_name = '{}_{}{}'.format(self.get_run_number(), self.DUT.Name, 'Flux' if with_flux else '')
        self.Draw.save_plots(save_name, show=show, ftype='png')
        return gc

    def draw_profile(self, bin_width=5, show=True):
        x, y = self.Data['timestamps'], self.Data['currents']
        return self.Draw.profile(x, y, make_bins(x[0], x[-1], bin_width), 'Leakage Current', x_tit='Time [hh:mm]', y_tit='Current [nA]', t_ax_off=0, markersize=.7, w=1.5,
                                 h=.75, lm=.08, y_off=.8, show=show, stats=set_entries())

    def draw_distribution(self, show=True):
        m, s = mean_sigma(self.Data['currents'], err=False)
        xmin, xmax = m - 4 * max(s, .1), m + 4 * max(s, .1)
        return self.Draw.distribution(self.Data['currents'], make_bins(xmin, xmax, self.Precision * 2), 'Current Distribution', show=show, x_tit='Current [nA]')

    def draw_flux_correlation(self, bin_fac=1, show=True):
        p1 = self.draw_profile(show=False)
        p2 = self.Ana.Tel.draw_flux(show=False)
        ybins = log_bins(int(sqrt(p1.GetNbinsX()) * bin_fac) * 3, .1, p1.GetMaximum() * 2)
        xbins = log_bins(int(sqrt(p1.GetNbinsX()) * bin_fac), .1, p2.GetMaximum() * 2)
        h = TH2F('gfcc', 'Correlation of Flux and Current', *(xbins + ybins))
        for i in range(p1.GetNbinsX()):
            if p1.GetBinContent(i) and p2.GetBinContent(i):
                h.Fill(p2.GetBinContent(i), p1.GetBinContent(i))
        format_histo(h, y_tit='Current [nA]', x_tit='Flux [kHz/cm^{2}', stats=0, y_off=1.2)
        self.Draw(h, 'FluxCurrent', draw_opt='colz', rm=.18, show=show, lm=.12, logy=True, logx=True)

    def draw_iv(self, show=True):
        x, y = self.Data['voltages'], self.Data['currents']
        return self.Draw.graph(x, y, 0, self.Precision, title='I-V Curve for {}'.format(self.DUT.Name), x_tit='Voltage [V]', y_tit='Current [nA]', show=show)
    # endregion PLOTTING
    # ----------------------------------------


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
