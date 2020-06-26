# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on May 19th 2016 by M. Reichmann
# --------------------------------------------------------

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gROOT, TF1, TFile, TMath

import pickle
from collections import OrderedDict
from copy import deepcopy
from datetime import datetime, timedelta
from multiprocessing import Pool, cpu_count
from subprocess import call
from sys import stdout
from threading import Thread
from time import time, sleep

from gtts import gTTS
from numpy import sqrt, array, average, mean, arange, log10, concatenate, where, any, count_nonzero, full, ndarray, histogram, searchsorted, cumsum, exp, sin, cos, arctan
from os import makedirs, _exit, remove, devnull
from os import path as pth
from os.path import dirname, realpath
from pytz import timezone, utc
from termcolor import colored
from uncertainties import ufloat
from uncertainties.core import Variable, AffineScalarFunc
from argparse import ArgumentParser
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from ConfigParser import ConfigParser
from scipy.optimize import curve_fit
from scipy import constants
import h5py
from functools import partial
from Queue import Queue
from json import load

OFF = False
ON = True
DEGREE_SIGN = u'\N{DEGREE SIGN}'

M_PI = 139.57018  # MeV/c^2
M_MU = constants.physical_constants['muon mass'][0] / constants.e * constants.c**2 / 1e6
M_E = constants.m_e / constants.e * constants.c**2 / 1e6
M_P = constants.m_p / constants.e * constants.c**2 / 1e6
TAU_PI = 26.033  # ns


# ==============================================
# UTILITY FUNCTIONS
# ==============================================
def get_t_str():
    return datetime.now().strftime('%H:%M:%S')


def log_warning(msg):
    print '{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head=colored('WARNING:', 'yellow'))


def log_critical(msg):
    print '{head} {t} --> {msg}\n'.format(t=get_t_str(), msg=msg, head=colored('CRITICAL:', 'red'))
    _exit(1)


def critical(msg):
    log_critical(msg)


def warning(msg):
    log_warning(msg)


def info(msg, next_line=True, blank_lines=0, prnt=True):
    return log_info(msg, next_line, blank_lines, prnt)


def log_info(msg, next_line=True, blank_lines=0, prnt=True):
    t1 = time()
    if prnt:
        t = datetime.now().strftime('%H:%M:%S')
        print '{bl}\r{head} {t} --> {msg}'.format(head=colored('INFO:', 'cyan', attrs=['dark']), t=t, msg=msg, bl='\n' * blank_lines),
        stdout.flush()
        if next_line:
            print
    return t1


def add_to_info(t, msg='Done', prnt=True):
    if prnt:
        print '{m} ({t:2.2f} s)'.format(m=msg, t=time() - t)


def set_root_warnings(status):
    gROOT.ProcessLine('gErrorIgnoreLevel = {e};'.format(e='0' if status else 'kError'))


def set_root_output(status=True):
    gROOT.SetBatch(not status)
    set_root_warnings(status)


def scale_margins(gr1, gr2):
    ymin1, ymax1 = gr1.GetYaxis().GetXmin(), gr1.GetYaxis().GetXmax()
    ymin2, ymax2 = gr2.GetYaxis().GetXmin(), gr2.GetYaxis().GetXmax()
    ymin = ymin1 if ymin1 < ymin2 else ymin2
    ymax = ymax1 if ymax1 > ymax2 else ymax2
    return ymin, ymax


def untitle(string):
    s = ''
    for word in string.split(' '):
        if word:
            s += word[0].lower() + word[1:] + ' '
    return s.strip(' ')


def draw_frame(pad, x, y, base=False, x_tit='', y_tit='', y_off=1, x_off=1):
    pad.cd()
    fr = pad.DrawFrame(x[0], y[0], x[1], y[1])
    pad.Modified()
    fr.GetXaxis().SetTitleOffset(x_off)
    fr.GetYaxis().SetTitleOffset(y_off)
    format_base_frame(fr, x_tit, y_tit) if base else format_transparent_frame(fr)


def format_transparent_frame(frame):
    fr = frame
    fr.GetXaxis().SetTickLength(0)
    fr.GetYaxis().SetTickLength(0)
    fr.GetXaxis().SetLabelOffset(99)
    fr.GetYaxis().SetLabelOffset(99)
    fr.SetLineColor(0)


def format_base_frame(frame, x_tit, y_tit):
    fr = frame
    # X-axis
    fr.GetXaxis().SetTitle(x_tit)
    # Y-axis
    fr.GetYaxis().SetTitle(y_tit)


def round_down_to(num, val):
    return int(num) / val * val


def round_up_to(num, val):
    return int(num) / val * val + val


def interpolate_two_points(x1, y1, x2, y2, name=''):
    # f = p1*x + p0
    p1 = (y1 - y2) / (x1 - x2)
    p0 = y1 - x1 * p1
    w = abs(x2 - x1)
    fit_range = array(sorted([x1, x2])) + [-w / 3., w / 3.]
    f = TF1('fpol1{}'.format(name), 'pol1', *fit_range)
    f.SetParameters(p0, p1)
    return f


def interpolate_x(x1, x2, y1, y2, y):
    p1 = get_p1(x1, x2, y1, y2)
    p0 = get_p0(x1, y1, p1)
    return (y - p0) / p1 if p1 else 0


def get_p1(x1, x2, y1, y2):
    return (y1 - y2) / (x1 - x2) if x1 != x2 else 0


def get_p0(x1, y1, p1):
    return y1 - x1 * p1


def move_element(odict, thekey, newpos):
    odict[thekey] = odict.pop(thekey)
    for i, (key, value) in enumerate(odict.items()):
        if key != thekey and i >= newpos:
            odict[key] = odict.pop(key)
    return odict


def make_transparent(pad):
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.SetFrameFillStyle(4000)


def hide_axis(axis):
    axis.SetTickLength(0)
    axis.SetLabelOffset(99)
    axis.SetTitleOffset(99)


def set_graph_color(mg, color=None, do_marker=True, marker_size=None):
    try:
        for g in mg.GetListOfGraphs():
            set_graph_color(g, color=color, do_marker=do_marker, marker_size=marker_size)
    except AttributeError:
        mg.SetLineColor(color)
        if do_marker:
            mg.SetMarkerColor(color)
        if marker_size:
            mg.SetMarkerSize(marker_size)


def get_bias_root_string(biases):
    if type(biases) == list:
        retval = ''
        for b in biases:
            if -1 * b in biases:
                b_str = '+-{:4.0f} V'.format(abs(b))
            else:
                b_str = '{:+4.0f} V'.format(b)
            if b_str not in retval:
                if retval != '':
                    retval += ', '
                retval += b_str
    elif type(biases) in [float, int]:
        retval = '{:+5.0f} V'.format(biases)
    else:
        retval = '{}'.format(biases)
    retval = retval.replace('+-', '#pm')
    retval = retval.replace('+/-', '#pm')
    retval = retval.replace('+', '#plus')
    retval = retval.replace('-', '#minus')
    return retval


def resize_markers(mg, default_size=0, marker_sizes=None):
    marker_sizes = {} if marker_sizes is None else marker_sizes
    try:
        for g in mg.GetListOfGraphs():
            resize_markers(g, default_size=default_size, marker_sizes=marker_sizes)
    except AttributeError:
        mg.SetMarkerSize(marker_sizes.get(mg.GetName(), default_size))


def move_legend(leg, x1, y1):
    xdiff = leg.GetX2NDC() - leg.GetX1NDC()
    ydiff = leg.GetY2NDC() - leg.GetY1NDC()
    leg.SetX1NDC(x1)
    leg.SetX2NDC(x1 + xdiff)
    leg.SetY1NDC(y1)
    leg.SetY2NDC(y1 + ydiff)


def scale_legend(leg, txt_size=None, width=None, height=None):
    sleep(.05)
    leg.SetY2NDC(height) if height is not None else do_nothing()
    leg.SetX2NDC(width) if width is not None else do_nothing()
    leg.SetTextSize(txt_size) if txt_size is not None else do_nothing()


def make_irr_string(val):
    if '?' in val:
        return val
    if not float(val):
        return 'nonirradiated'
    val, power = [float(i) for i in val.split('e')]
    return '{v:1.1f}#upoint10^{p} n/cm^{{2}}'.format(v=val, p='{{{}}}'.format(int(power)))


def increased_range(ran, fac_bot=0., fac_top=0.):
    return [(1 + fac_bot) * ran[0] - fac_bot * ran[1], (1 + fac_top) * ran[1] - fac_top * ran[0]]


def calc_weighted_mean(means, sigmas):
    weights = map(lambda x: x ** (-2), sigmas)
    variance = 1 / sum(weights)
    mean_ = sum(map(lambda x, y: x * y, means, weights))
    return mean_ * variance, sqrt(variance)


def mean_sigma(values, weights=None):
    """ Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape. """
    if len(values) == 1:
        value = make_ufloat(values[0])
        return value.n, value.s
    weights = full(len(values), 1) if weights is None else weights
    if type(values[0]) in [Variable, AffineScalarFunc]:
        errors = array([v.s for v in values])
        weights = full(errors.size, 1) if all(errors == errors[0]) else [1 / e if e else 0 for e in errors]
        values = array([v.n for v in values], 'd')
    if all(weight == 0 for weight in weights):
        return [0, 0]
    avrg = average(values, weights=weights)
    variance = average((values - avrg) ** 2, weights=weights)  # Fast and numerically precise
    return avrg, sqrt(variance)


def make_latex_table_row(row, hline=False):
    return '{0}\\\\{1}\n'.format('\t& '.join(row), '\\hline' if hline else '')


def make_latex_table(header, cols, endline=False):
    header = ['\\textbf{0}{1}{2}'.format('{', head, '}') if head else '' for head in header]
    size = max([len(col) for col in cols])
    rows = []
    for col in cols:
        for i in xrange(size):
            if len(rows) < size:
                rows.append([])
            rows[i].append(col[i] if len(col) > i else '')
    out = '\\toprule\n'
    out += make_latex_table_row(header, hline=True)
    for i, row in enumerate(rows, 1):
        out += make_latex_table_row(row) if i != len(rows) else make_latex_table_row(row, endline)
    out += '\\bottomrule' if not endline else ''
    return out


def make_dia_str(dia):
    dia = dia.replace('-', '')
    return '{0}{1}'.format(dia[0].title(), dia[1:])


def make_list(value):
    return array([value]).flatten()


def file_exists(path, warn=False):
    if not pth.isfile(path):
        warning('File "{}" does not exist!'.format(path)) if warn else do_nothing()
        return False
    return True


def dir_exists(path):
    return pth.isdir(path)


def ensure_dir(path):
    if not pth.exists(path):
        info('Creating directory: {d}'.format(d=path))
        makedirs(path)
    return path


def make_col_str(col):
    return '{0:2d}'.format(int(col)) if int(col) > 1 else '{0:3.1f}'.format(col)


def print_banner(msg, symbol='~', new_lines=1, color=None):
    msg = '{} |'.format(msg)
    print colored('{n}{delim}\n{msg}\n{delim}{n}'.format(delim=len(str(msg)) * symbol, msg=msg, n='\n' * new_lines), color)


def print_small_banner(msg, symbol='-', color=None):
    print colored('\n{delim}\n{msg}\n'.format(delim=len(str(msg)) * symbol, msg=msg), color)


def print_elapsed_time(start, what='This', show=True, color=None):
    string = 'Elapsed time for {w}: {t}'.format(t=get_elapsed_time(start), w=what)
    print_banner(string, color=color) if show else do_nothing()
    return string


def get_elapsed_time(start):
    t = datetime.fromtimestamp(time() - start)
    return '{}.{:02.0f}'.format(t.strftime('%M:%S'), t.microsecond / 10000)


def conv_log_time(time_str, strg=False):
    t = datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=utc).astimezone(timezone('Europe/Zurich'))
    return t.strftime('%b %d, %H:%M:%S') if strg else t


def has_bit(num, bit):
    assert (num >= 0 and type(num) is int), 'num has to be non negative int'
    return bool(num & 1 << bit)


def make_tc_str(tc, long_=True, data=False):
    tc_data = str(tc).split('-')
    sub_string = '-{0}'.format(tc_data[-1]) if len(tc_data) > 1 else ''
    if data:
        return datetime.strptime(tc_data[0], '%Y%m').strftime('psi_%Y_%m')
    elif tc_data[0][0].isdigit():
        return '{tc}{s}'.format(tc=datetime.strptime(tc_data[0], '%Y%m').strftime('%B %Y' if long_ else '%b%y'), s=sub_string)
    else:
        return '{tc}{s}'.format(tc=datetime.strptime(tc_data[0], '%b%y').strftime('%Y%m' if long_ else '%B %Y'), s=sub_string)


def tc_to_str(tc, short=True):
    tc_str = str(tc).split('-')[0]
    sub_str = '-{}'.format(tc.split('-')[-1]) if '-' in str(tc) else ''
    return '{tc}{s}'.format(tc=datetime.strptime(tc_str, '%Y%m').strftime('%b%y' if short else '%B %Y'), s=sub_str)


def make_flux_string(rate):
    unit = '{}/cm^{{2}}'.format('MHz' if rate > 1000 else 'kHz')
    rate /= 1000. if rate > 1000 else 1
    return '{: 2.1f} {}'.format(rate, unit)


def make_bias_str(bias):
    return '{s}{bias}V'.format(bias=int(bias), s='+' if bias > 0 else '')


def make_runplan_string(nr):
    nr = str(nr)
    return nr.zfill(2) if len(nr) <= 2 else nr.zfill(4)


def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def isint(x):
    try:
        return float(x) == int(x)
    except (ValueError, TypeError):
        return False


def set_drawing_range(h, legend=True, lfac=None, rfac=None, thresh=10):
    for i in xrange(1, 4):
        h.SetBinContent(i, 0)
    range_ = [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(thresh), h.FindLastBinAbove(thresh)]]
    lfac = lfac if lfac is not None else .2
    rfac = rfac if rfac is not None else .55 if legend else .1
    h.GetXaxis().SetRangeUser(*increased_range(range_, lfac, rfac))


def set_time_axis(histo, form='%H:%M', off=0):
    histo.GetXaxis().SetTimeFormat(form)
    histo.GetXaxis().SetTimeOffset(-off - 3600 if off else 0)
    histo.GetXaxis().SetTimeDisplay(1)


def find_mpv_fwhm(histo, bins=15):
    max_bin = histo.GetMaximumBin()
    fit = TF1('fit', 'gaus', 0, 500)
    histo.Fit('fit', 'qs0', '', histo.GetBinCenter(max_bin - bins), histo.GetBinCenter(max_bin + bins))
    mpv = make_ufloat(fit.GetParameter(1), fit.GetparError(1))
    fwhm = histo.FindLastBinAbove(fit(mpv.n) / 2) - histo.FindFirstBinAbove(fit(mpv.n) / 2)
    return mpv, fwhm, mpv / fwhm


def get_fwhm(h):
    max_v = h.GetMaximum()
    return h.GetBinCenter(h.FindLastBinAbove(max_v / 2)) - h.GetBinCenter(h.FindFirstBinAbove(max_v / 2))


def fit_fwhm(histo, fitfunc='gaus', do_fwhm=True, draw=False):
    h = histo
    if do_fwhm:
        peak_pos = h.GetBinCenter(h.GetMaximumBin())
        bin1 = h.FindFirstBinAbove(h.GetMaximum() / 2)
        bin2 = h.FindLastBinAbove(h.GetMaximum() / 2)
        fwhm = h.GetBinLowEdge(bin2 + 2) - h.GetBinLowEdge(bin1 - 1)
        option = 'qs' if draw else 'qs0'
        fit = h.Fit(fitfunc, option, '', peak_pos - fwhm / 2, peak_pos + fwhm / 2)
    else:
        fit = h.Fit(fitfunc, 'qs')
    return FitRes(fit)


def scale_histo(histo, value=None, to_max=False, x_range=None):
    h = histo
    maximum = h.GetBinContent(h.GetMaximumBin())
    if x_range is not None:
        h.GetXaxis().SetRangeUser(*x_range) if x_range is not None else do_nothing()
        maximum = h.GetBinContent(h.GetMaximumBin())
        h.GetXaxis().UnZoom()
    value = maximum if to_max else value
    if value:
        h.Scale(1. / value)
    return h


def normalise_histo(histo, x_range=None, from_min=False):
    h = histo
    x_axis = h.GetXaxis()
    x_axis.SetRangeUser(*x_range) if x_range is not None else do_nothing()
    min_bin = h.GetMinimumBin() if from_min else 0
    integral = h.Integral(min_bin, h.GetNbinsX() - 1)
    return scale_histo(h, integral)


def make_cut_string(cut, n):
    return '{n}Cuts'.format(n=str(n).zfill(2)) if cut is not None else ''


def get_resolution():
    try:
        from screeninfo import get_monitors
        m = get_monitors()
        return round_down_to(m[0].height, 500)
    except Exception as exc:
        log_warning('Could not get resolution! Using default ...\n\t{e}'.format(e=exc))
        return 1000


def print_table(rows, header=None, prnt=True):
    t = array(rows, dtype=str) if header is None else concatenate((array([header], dtype=str), array(rows, dtype=str)))
    col_width = [len(max(t[:, i], key=len)) for i in xrange(t.shape[1])]
    total_width = sum(col_width) + len(col_width) * 3 + 1
    hline = '{}'.format('~' * total_width)
    if prnt:
        for i, row in enumerate(t):
            if i in [0, 1, t.size]:
                print hline
            print '| {r} |'.format(r=' | '.join(word.ljust(n) for word, n in zip(row, col_width)))
        print '{}\n'.format(hline)
    return rows


def get_base_dir():
    return dirname(dirname(realpath(__file__)))


def do_pickle(path, func, value=None, redo=False, *args, **kwargs):
    if value is not None:
        with open(path, 'w') as f:
            pickle.dump(value, f)
        return value
    try:
        if file_exists(path) and not redo:
            with open(path, 'r') as f:
                return pickle.load(f)
    except ImportError:
        pass
    ret_val = func(*args, **kwargs)
    with open(path, 'w') as f:
        pickle.dump(ret_val, f)
    return ret_val


def do_hdf5(path, func, redo=False, *args, **kwargs):
    if file_exists(path) and redo:
        remove_file(path)
    if file_exists(path) and not redo:
        return h5py.File(path, 'r')['data']
    else:
        data = func(*args, **kwargs)
        f = h5py.File(path, 'w')
        f.create_dataset('data', data=data)
        return f['data']


def fit_poissoni(h, p0=5000, p1=1, name='f_poiss', show=True):
    fit = TF1(name, '[0] * TMath::PoissonI(x, [1])', 0, 30)
    fit.SetParNames('Constant', 'Lambda')
    fit.SetParameters(p0, p1)
    h.Fit(fit, 'q{0}'.format('' if show else 0))
    fit.Draw('same')
    return fit


def int_to_roman(integer):
    """ Convert an integer to Roman numerals. """
    if type(integer) != int:
        raise (TypeError, 'expected integer, got {t}'.format(t=type(integer)))
    if not 0 < integer < 4000:
        raise (ValueError, 'Argument must be between 1 and 3999')
    dic = OrderedDict([(1000, 'M'), (900, 'CM'), (500, 'D'), (400, 'CD'), (100, 'C'), (90, 'XC'), (50, 'L'), (40, 'XL'),
                       (10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')])
    result = ''
    for i, num in dic.iteritems():
        count = int(integer / i)
        result += num * count
        integer -= i * count
    return result


def set_z_range(zmin, zmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetZaxis().SetRangeUser(zmin, zmax)


def set_axes_range(xmin, xmax, ymin, ymax):
    set_x_range(xmin, xmax)
    set_y_range(ymin, ymax)


def set_x_range(xmin, xmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetXaxis().SetRangeUser(xmin, xmax)


def set_y_range(ymin, ymax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetYaxis().SetRangeUser(ymin, ymax)


def remove_letters(string):
    return filter(lambda x: x.isdigit(), string)


def remove_digits(string):
    return filter(lambda x: not x.isdigit(), string)


def get_last_canvas():
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        log_warning('There is no canvas is in the list...')


def close_last_canvas():
    get_last_canvas().Close()


def get_object(name):
    return gROOT.FindObject(name)


def do(fs, pars, exe=-1):
    fs, pars = ([fs], [pars]) if type(fs) is not list else (fs, pars)
    exe = pars if exe == -1 else [exe]
    for f, p, e in zip(fs, pars, exe):
        f(p) if e is not None else do_nothing()


def choose(v, default, decider='None', *args, **kwargs):
    use_default = decider is None if decider != 'None' else v is None
    if callable(default) and use_default:
        default = default(*args, **kwargs)
    return default if use_default else v


def average_list(lst, n):
    return [mean(lst[i:i + n]) for i in arange(0, len(lst), n)] if n > 1 else lst


def log_bins(n_bins, min_val, max_val):
    width = (log10(max_val) - log10(min_val)) / float(n_bins)
    return [n_bins, array([pow(10, log10(min_val) + i * width) for i in xrange(n_bins + 1)])]


def make_ufloat(tup, par=0):
    if type(tup) in [Variable, AffineScalarFunc]:
        return tup
    if isinstance(tup, FitRes):
        return ufloat(tup.Parameter(par), tup.ParError(par))
    if type(tup) in [tuple, list, ndarray]:
        return ufloat(*tup) if type(tup[0]) not in [Variable, AffineScalarFunc] else ufloat(tup[0].n, tup[1].n)
    return ufloat(tup, 0)


def find_graph_margins(graphs):
    graphs = deepcopy(graphs)
    for i, g in enumerate(graphs):
        if g.Class().GetName().startswith('TMulti'):
            graphs[i] = g.GetListOfGraphs()[0]
    return min([min(gr.GetY()[i] for i in xrange(gr.GetN()) if gr.GetY()[i] >= 0.01) for gr in graphs]), max([TMath.MaxElement(gr.GetN(), gr.GetY()) for gr in graphs])


def get_quantiles(values, bins):
    entries, bins = histogram(values, bins)
    bin_width = (bins[1] - bins[0]) / 2.
    return (bins + bin_width)[searchsorted(cumsum(entries), arange(0, 100, .01) * sum(entries))]


def load_root_files(sel, init=True):
    threads = {}
    for run in sel.get_selected_runs():
        thread = MyThread(sel, run, init)
        thread.start()
        threads[run] = thread
    while any([thread.isAlive() for thread in threads.itervalues()]) and init:
        sleep(.1)
    if init:
        pool = Pool(len(threads))
        results = [pool.apply_async(get_time_vec, (thread.Selection, thread.Run)) for thread in threads.itervalues()]
        times = [result.get(60) for result in results]
        for thread, t in zip(threads.itervalues(), times):
            thread.Time = t
        pool.close()
    return threads


class MyThread(Thread):
    def __init__(self, sel, run, init=True):
        Thread.__init__(self)
        self.Load = init
        self.Run = run
        self.Selection = sel
        self.File = None
        self.Tree = None
        self.Tuple = None
        self.Time = None

    def run(self):
        self.load_tree()

    def load_tree(self):
        if not self.Load:
            self.Tuple = False
            return
        info('Loading run {r}'.format(r=self.Run), next_line=False)
        file_path = self.Selection.get_final_file_path(self.Run)
        if file_exists(file_path):
            self.File = TFile(file_path)
            self.Tree = self.File.Get('tree')
            self.Tuple = (self.File, self.Tree)
        self.Tree.SetEstimate(-1)
        return self.Tree


def get_time_vec(sel, run=None):
    if run is None:
        tree = sel
    else:
        t = MyThread(sel, run)
        tree = t.load_tree()
    if tree is None:
        return
    tree.SetEstimate(-1)
    entries = tree.Draw('time', '', 'goff')
    time_vec = get_root_vec(tree, entries)
    fill_empty_time_entries(time_vec)
    time_vec = correct_time(time_vec, run)
    return time_vec


def fill_empty_time_entries(times):
    empty_events = where(times == -1)[0]
    times[empty_events] = times[empty_events.size]
    return times


def time_stamp(dt, off=None):
    t = float(dt.strftime('%s'))
    return t if off is None else t - (off if off > 1 else dt.utcoffset().seconds)


def correct_time(times, run):
    i_off = where(times[:-1] > times[1:])[0]
    i_off = i_off[0] if i_off else None
    if i_off is not None:
        warning('Need to correct timing vector for run {}\n'.format(run))
        diff = times[i_off + 1] - times[i_off]
        i_off += 1  # because range does not include the last index
        return concatenate((times[:i_off], times[i_off:] - diff + 500))  # one TU step should be 500 ms
    return times


def say(txt):
    tts = gTTS(text=txt.decode('utf-8'), lang='en')
    tts.save('good.mp3')
    with open(devnull, 'w') as FNULL:
        call(['mpg321', 'good.mp3'], stdout=FNULL)
    remove('good.mp3')


def remove_file(file_path):
    if file_exists(file_path):
        log_warning('removing {}'.format(file_path))
        remove(file_path)


def get_running_time(t):
    now = datetime.fromtimestamp(time() - t) - timedelta(hours=1)
    return now.strftime('%H:%M:%S')


def del_rootobj(obj):
    if obj is None:
        return
    try:
        if obj.IsA().GetName() != 'TCanvas':
            obj.Delete()
    except AttributeError:
        pass


def get_root_vec(tree, n=0, ind=0, dtype=None, var=None, cut=''):
    if var is not None:
        n = tree.Draw(var, cut, 'goff')
    vec = tree.GetVal(ind)
    vec.SetSize(n)
    return array(vec, dtype=dtype)


def get_root_vecs(tree, n, n_ind, dtype=None):
    return [get_root_vec(tree, n, i, dtype) for i in xrange(n_ind)]


def get_arg(arg, default):
    return default if arg is None else arg


def calc_eff(k=0, n=0, values=None):
    values = array(values) if values is not None else None
    k = float(k if values is None else count_nonzero(values))
    n = float(n if values is None else values.size)
    return make_ufloat((100 * (k + 1) / (n + 2), 100 * sqrt(((k + 1) / (n + 2) * (k + 2) / (n + 3) - ((k + 1) ** 2) / ((n + 2) ** 2)))))


def init_argparser(run=None, tc=None, dut=False, tree=False, has_verbose=False, has_collection=False, return_parser=False):
    p = ArgumentParser()
    p.add_argument('run' if not has_collection else 'runplan', nargs='?', default=run, type=str if has_collection else int, help='run {}'.format('number' if not has_collection else 'plan'))
    p.add_argument('dut', nargs='?', default=dut, type=int, help='diamond number [default: 1] (choose from 1,2,...)') if dut or dut is None else do_nothing()
    p.add_argument('-tc', '--testcampaign', nargs='?', default=tc, help='YYYYMM beam test [default in main.ini]')
    p.add_argument('-v', '--verbose', action='store_false') if has_verbose else do_nothing()
    p.add_argument('-t', '--tree', action='store_false', help='do not load the ROOT TTree') if tree else do_nothing()
    p.add_argument('-c', '--collection', action='store_true', help='begin analysis collection') if has_collection else do_nothing()
    return p if return_parser else p.parse_args()


def load_parser(name):
    p = ConfigParser()
    p.read(name)
    return p


def load_json(name):
    with open(name) as f:
        return load(f)


def measure_time(f, rep=1, *args, **kwargs):
    t = info('Measuring time of method {}:'.format(f.__name__), next_line=False)
    for _ in xrange(int(rep)):
        f(*args, **kwargs)
    add_to_info(t, '')


def u_to_str(v, prec=2):
    return '{{:1.{0}f}} ({{:1.{0}f}})'.format(prec).format(v.n, v.s)


class FitRes:
    def __init__(self, f=None):
        self.Pars = [None]
        self.Errors = [None]
        self.Names = None
        self.vChi2 = None
        self.vNdf = None
        if f is None:
            return
        is_tf1 = 'TF1' in f.ClassName()
        self.NPar = f.GetNpar() if is_tf1 else f.NPar()
        if self.NPar < 1:
            return
        self.Pars = [f.GetParameter(i) for i in xrange(self.NPar)] if is_tf1 else list(f.Parameters())
        self.Errors = [f.GetParError(i) for i in xrange(self.NPar)] if is_tf1 else list(f.Errors())
        self.Names = [f.GetParName(i) if is_tf1 else f.ParName(i) for i in xrange(self.NPar)]
        self.vChi2 = f.GetChisquare() if is_tf1 else f.Chi2()
        self.vNdf = f.GetNDF() if is_tf1 else f.Ndf()

    def Parameter(self, arg):
        return self.Pars[arg]

    def ParError(self, arg):
        return self.Errors[arg]

    def ParName(self, arg):
        return self.Names[arg]

    def Chi2(self):
        return self.vChi2

    def Ndf(self):
        return self.vNdf


class PBar:
    def __init__(self):
        self.PBar = None
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.Step = 0

    def start(self, n):
        self.Step = 0
        self.PBar = ProgressBar(widgets=self.Widgets, maxval=n).start()

    def update(self, i=None):
        i = self.Step if i is None else i
        if i >= self.PBar.maxval:
            return
        self.PBar.update(i + 1)
        self.Step += 1
        if i == self.PBar.maxval - 1:
            self.finish()

    def finish(self):
        self.PBar.finish()


def gauss(x, scale, mean_, sigma, off=0):
    return scale * exp(-.5 * ((x - mean_) / sigma) ** 2) + off


def fit_data(f, y, x=None, p=None):
    x = arange(y.shape[0]) if x is None else x
    return curve_fit(f, x, y, p0=p)


def calc_speed(p, m):
    return 1 / sqrt(1 + m**2 / p**2)


def t_diff(s, p, m1, m2):
    return s * (1 / calc_speed(p, m1) - 1 / calc_speed(p, m2)) / constants.c * 1e9


def e_kin(p, m):
    return sqrt(p**2 + m**2) - m


def gamma_factor(v):
    return 1 / sqrt(1 - v**2)


def momentum(m, v):
    return m * v * gamma_factor(v)


def decay_ratio(p, m, d, tau):
    return exp(-d * m / (tau * 1e-9 * p * constants.c))


def decay_momentum(m, m1, m2=0):
    return sqrt((m**2 + m1**2 + m2**2)**2 - 4 * m**2 * m1**2) / (2 * m)


def decay_energy(m, m1, m2=0):
    return (m**2 + m1**2 - m2**2) / (2 * m)


def decay_angle(theta, p, m, m1, m2=0):
    p1 = decay_momentum(m, m1, m2)
    v = calc_speed(p, m)
    return arctan(p1 * sin(theta) / (gamma_factor(v) * (p1 * cos(theta) + v * decay_energy(m, m1, m2))))


def multi_threading(lst, timeout=60 * 60 * 2):
    """ runs several threads in parallel. [lst] must contain tuples of the methods and the arguments as list."""
    t0 = info('Run multithreading on {} tasks ... '.format(len(lst)), next_line=False)
    lst = [(f, [], {}) for f in lst] if type(lst[0]) not in [list, tuple, ndarray] else lst
    if len(lst[0]) == 2:
        lst = [(f, args, {}) for f, args in lst] if type(lst[0][1]) not in [dict, OrderedDict] else [(f, [], d) for f, d in lst]
    threads = []
    queue = Queue()  # use a queue to get the results
    for f, args, kwargs in lst:
        t = Thread(target=lambda q, a, k: q.put(f(*a, **k)), args=(queue, make_list(args), kwargs))
        t.start()
        threads.append(t)
    for thread in threads:
        thread.join(timeout)
    add_to_info(t0)
    return [queue.get() for _ in range(queue.qsize())]


def get_attribute(instance, string):
    if '.' in string:
        s = string.split('.')
        return getattr(getattr(instance, s[0]), s[1])
    return getattr(instance, string)


def parallelise(f, args_list, timeout=60 * 60):
    t = info('Run parallelisation on {} tasks ... '.format(len(args_list)), next_line=False)
    pool = Pool(cpu_count())
    workers = [pool.apply_async(f, make_list(args)) for args in args_list]
    results = [worker.get(timeout) for worker in workers]
    add_to_info(t)
    return results


def parallelise_instance(instances, method, args, timeout=60 * 60):
    t = info('Run parallelisation on {} tasks ... '.format(len(args)), next_line=False)
    pool = Pool(cpu_count())
    # tasks = [partial(call_it, make_list(instances)[0], method.__name__, *make_list(arg)) for arg in args]
    tasks = [partial(call_it, instance, method.__name__, *make_list(arg)) for instance, arg in zip(instances, args)]
    workers = [pool.apply_async(task) for task in tasks]
    results = [worker.get(timeout) for worker in workers]
    add_to_info(t)
    return results


def call_it(instance, name, *args, **kwargs):
    """indirect caller for instance methods and multiprocessing"""
    return getattr(instance, name)(*args, **kwargs)


def do_nothing():
    pass
