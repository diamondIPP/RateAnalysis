# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on May 19th 2016 by M. Reichmann
# --------------------------------------------------------

from datetime import datetime, timedelta
from termcolor import colored
from numpy import sqrt, array, average, mean, arange, log10
from os import makedirs, _exit, remove, devnull
from os import path as pth
from os.path import basename, join, dirname, expanduser
from time import time, sleep
from collections import OrderedDict
import pickle
from threading import Thread
from multiprocessing import Pool
from uncertainties import ufloat
from uncertainties.core import Variable, AffineScalarFunc
from copy import deepcopy
from gtts import gTTS
from subprocess import call
from pytz import timezone, utc
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gStyle, gROOT, TF1, TColor, TFile, TMath


OFF = False
ON = True


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


def log_message(msg, overlay=False, prnt=True):
    overlay = '\033[1A\r' if overlay else ''
    info = '{ov}{head} {t} --> {msg}{end}'.format(t=get_t_str(), msg=msg, ov=overlay, end=' ' * 10 if overlay else '', head=colored('INFO:', 'cyan', attrs=['dark']))
    if prnt:
        print info
    return info


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


def set_statbox(x=.95, y=.88, w=.16, entries=3, only_fit=False, fit=False, only_entries=False, opt=None, form=None):
    if only_fit or fit:
        gStyle.SetOptStat(1111 if fit else 0011)
        gStyle.SetOptFit(1)
    if only_entries:
        gStyle.SetOptStat(1000000010 if not only_fit else 1000000011)
        entries = 6 if entries == 3 else entries
    gStyle.SetOptStat(opt) if opt is not None else do_nothing()
    gStyle.SetFitFormat(form) if form is not None else do_nothing()
    gStyle.SetStatX(x)
    gStyle.SetStatY(y)
    gStyle.SetStatW(w)
    gStyle.SetStatH(.02 * entries)


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


def scale_multigraph(mg, val=1, to_low_flux=False):
    if val is None:
        return
    g = mg.GetListOfGraphs()[0]
    points = {g.GetX()[i]: g.GetY()[i] for i in xrange(g.GetN())}
    y_vals = [make_ufloat((g.GetY()[i], g.GetEY()[i])) for i in xrange(g.GetN())]
    m, s = mean_sigma(y_vals)
    scale = val / (points[min(points)] if to_low_flux else m)
    for gr in mg.GetListOfGraphs():
        scale_graph(gr, scale)
    for i, l in enumerate(mg.GetListOfGraphs()[0].GetListOfFunctions()):
        y, ey = g.GetY()[i], g.GetErrorY(i)
        l.SetY(y - ey - .003)


def scale_graph(gr, scale=None, val=1, to_low_flux=False):
    if scale is None:
        points = {gr.GetX()[i]: gr.GetY()[i] for i in xrange(gr.GetN())}
        y_vals = [make_ufloat((gr.GetY()[i], gr.GetEY()[i])) for i in xrange(gr.GetN())]
        m, s = mean_sigma(y_vals)
        scale = val / (points[min(points)] if to_low_flux else m)
    for i in xrange(gr.GetN()):
        gr.SetPoint(i, gr.GetX()[i], gr.GetY()[i] * scale)
        try:
            gr.SetPointError(i, gr.GetErrorX(i), gr.GetErrorY(i) * scale)
        except Exception as err:
            log_warning('Error in scale multigraph: {err}'.format(err=err))


def interpolate_two_points(x1, y1, x2, y2):
    # f = p1*x + p0
    p1 = (y1 - y2) / (x1 - x2)
    p0 = y1 - x1 * p1
    f = TF1('fpol1', 'pol1', -1000, 1000)
    f.SetParameters(p0, p1)
    return f


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


def get_graph_data(g):
    n = g.GetN()
    x = [g.GetX()[i] for i in xrange(n)]
    y = [g.GetY()[i] for i in xrange(n)]
    ex = [g.GetEX()[i] for i in xrange(n)]
    ey = [g.GetEY()[i] for i in xrange(n)]
    return n, x, y, ex, ey


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


def move_legend(l, x1, y1):
    xdiff = l.GetX2NDC() - l.GetX1NDC()
    ydiff = l.GetY2NDC() - l.GetY1NDC()
    l.SetX1NDC(x1)
    l.SetX2NDC(x1 + xdiff)
    l.SetY1NDC(y1)
    l.SetY2NDC(y1 + ydiff)


def scale_legend(l, txt_size=None, width=None, height=None):
    sleep(.05)
    l.SetY2NDC(height) if height is not None else do_nothing()
    l.SetX2NDC(width) if width is not None else do_nothing()
    l.SetTextSize(txt_size) if txt_size is not None else do_nothing()


def make_irr_string(val):
    if '?' in val:
        return val
    if not float(val):
        return 'nonirradiated'
    val, power = [float(i) for i in val.split('e')]
    return '{v:1.1f}#upoint10^{p} n/cm^{{2}}'.format(v=val, p='{{{}}}'.format(int(power)))


def increased_range(ran, fac_bot=0., fac_top=0.):
    return [(1 + fac_bot) * ran[0] - fac_bot * ran[1], (1 + fac_top) * ran[1] - fac_top * ran[0]]


def calc_mean(l):
    l = [float(i) for i in l]
    mean_ = sum(l) / len(l)
    mean2 = sum(map(lambda x: x ** 2, l)) / len(l)
    sigma = sqrt(mean2 - mean_ ** 2)
    return mean_, sigma


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
    weights = [1] * len(values) if weights is None else weights
    if type(values[0]) in [Variable, AffineScalarFunc]:
        weights = [1 / v.s for v in values]
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
    l = max([len(col) for col in cols])
    rows = []
    for col in cols:
        for i in xrange(l):
            if len(rows) < l:
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


def file_exists(path, warn=False):
    if not pth.isfile(path):
        warning('File "{}" does not exist!'.format(path)) if warn else do_nothing()
        return False
    return True


def dir_exists(path):
    return pth.isdir(path)


def ensure_dir(path):
    if not pth.exists(path):
        log_message('Creating directory: {d}'.format(d=path))
        makedirs(path)


def make_col_str(col):
    return '{0:2d}'.format(int(col)) if int(col) > 1 else '{0:3.1f}'.format(col)


def print_banner(msg, symbol='=', new_lines=True):
    print '{n}{delim}\n{msg}\n{delim}{n}'.format(delim=len(str(msg)) * symbol, msg=msg, n='\n' if new_lines else '')


def print_small_banner(msg, symbol='-'):
    print '\n{delim}\n{msg}\n'.format(delim=len(str(msg)) * symbol, msg=msg)


def print_elapsed_time(start, what='This', show=True):
    t = '{d}'.format(d=timedelta(seconds=time() - start)).split('.')
    string = 'Elapsed time for {w}: {d1}.{m:2d}'.format(d1=t[0], m=int(round(int(t[1][:3])) / 10.), w=what)
    print_banner(string) if show else do_nothing()
    return string


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


def make_rate_str(rate):
    unit = '{}/cm^{{2}}'.format('MHz' if rate > 1000 else 'kHz')
    rate = round(rate / 1000., 1) if rate > 1000 else int(round(rate, 0))
    return '{rate} {unit}'.format(rate=rate, unit=unit)


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
        a = float(x)
        b = int(a)
        return a == b
    except ValueError:
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


def print_table(rows, header=None):
    col_width = [len(max([header[i]] + [row[i] for row in rows], key=len)) for i in xrange(len(header))]
    total_width = sum(col_width) + len(col_width) * 3 + 1
    hline = '~' * total_width
    if header is not None:
        print hline
        print '| {r} |'.format(r=' | '.join(word.ljust(n) for word, n in zip(header, col_width)))
    print hline
    for row in rows:
        row = row if len(row) > 1 else [row[0]]
        print '| {r} |'.format(r=' | '.join(word.ljust(n) for word, n in zip(row, col_width)).ljust(len(hline) - 4))
    print hline


def get_base_dir():
    return join('/', *__file__.split('/')[1:3])


def server_is_mounted():
    return dir_exists(expanduser(join('~', 'mounts', 'psi', 'Diamonds')))


def do_pickle(path, func, value=None, params=None, redo=False):
    if value is not None:
        with open(path, 'w') as f:
            pickle.dump(value, f)
        return value
    if file_exists(path) and not redo:
        with open(path, 'r') as f:
            return pickle.load(f)
    else:
        ret_val = func() if params is None else func(params)
        with open(path, 'w') as f:
            pickle.dump(ret_val, f)
        return ret_val


def server_pickle(old_path, value):
    if server_is_mounted():
        picklepath = join(get_base_dir(), 'mounts', 'psi', 'Pickles', basename(dirname(old_path)), basename(old_path))
        do_pickle(picklepath, do_nothing, value)


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


def get_last_canvas():
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        log_warning('There is no canvas is in the list...')


def close_last_canvas():
    get_last_canvas().Close()


def get_object(name):
    return gROOT.FindObject(name)


def get_color_gradient(n):
    stops = array([0., .5, 1], 'd')
    green = array([0. / 255., 200. / 255., 80. / 255.], 'd')
    blue = array([0. / 255., 0. / 255., 0. / 255.], 'd')
    red = array([180. / 255., 200. / 255., 0. / 255.], 'd')
    # gStyle.SetNumberContours(20)
    bla = TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 255)
    color_table = [bla + ij for ij in xrange(255)]
    return color_table[0::(len(color_table) + 1) / n]


def do(fs, pars, exe=-1):
    fs, pars = ([fs], [pars]) if type(fs) is not list else (fs, pars)
    exe = pars if exe == -1 else [exe]
    for f, p, e in zip(fs, pars, exe):
        f(p) if e is not None else do_nothing()


def make_bias_str(bias):
    return '{s}{bias}V'.format(bias=int(bias), s='+' if bias > 0 else '')


def markers(i):
    return (range(20, 24) + [29, 33, 34])[i]


def average_list(lst, n):
    return [mean(lst[i:i+n]) for i in arange(0, len(lst), n)] if n > 1 else lst


def log_bins(n_bins, min_val, max_val):
    width = (log10(max_val) - log10(min_val)) / float(n_bins)
    return [n_bins, array([pow(10, log10(min_val) + i * width) for i in xrange(n_bins + 1)])]


def make_ufloat(tup, par=0):
    if type(tup) in [Variable, AffineScalarFunc]:
        return tup
    if isinstance(tup, FitRes):
        return ufloat(tup.Parameter(par), tup.ParError(par))
    return ufloat(tup[0], tup[1]) if type(tup) in [tuple, list] else ufloat(tup, 0)


def find_graph_margins(graphs):
    graphs = deepcopy(graphs)
    for i, g in enumerate(graphs):
        if g.Class().GetName().startswith('TMulti'):
            graphs[i] = g.GetListOfGraphs()[0]
    return min([min(gr.GetY()[i] for i in xrange(gr.GetN()) if gr.GetY()[i] >= 0.01) for gr in graphs]), max([TMath.MaxElement(gr.GetN(), gr.GetY()) for gr in graphs])


def load_root_files(sel, load=True):

    runs = sel.get_selected_runs()
    threads = {}
    for run in runs:
        thread = MyThread(sel, run, load)
        thread.start()
        threads[run] = thread
    while any(thread.isAlive() for thread in threads.itervalues()) and load:
        sleep(.1)
    if load:
        pool = Pool(len(threads))
        results = [pool.apply_async(get_time_vec, (thread.Selection, thread.Run)) for thread in threads.itervalues()]
        times = [result.get(60) for result in results]
        for thread, t in zip(threads.itervalues(), times):
            thread.Time = t
        pool.close()
    return threads


class MyThread (Thread):
    def __init__(self, sel, run, load=True):
        Thread.__init__(self)
        self.Load = load
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
        log_message('Loading run {r}'.format(r=self.Run), overlay=True)
        file_path = self.Selection.get_final_file_path(self.Run)
        if file_exists(file_path):
            self.File = TFile(file_path)
            self.Tree = self.File.Get('tree')
            self.Tuple = (self.File, self.Tree)
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
    time_vec = [tree.GetV1()[i] for i in xrange(entries)]
    fill_empty_time_entries(time_vec)
    if any(time_vec[i + 100] < time_vec[i] for i in xrange(0, len(time_vec) - 100, 100)):
        log_warning('Need to correct timing vector\n')
        time_vec = correct_time(time_vec)
    return time_vec


def time_stamp(t):
    return float(t.strftime('%s'))


def fill_empty_time_entries(times):
    first_valid = 0
    ind = 0
    for i, t in enumerate(times):
        if t != -1:
            first_valid = t
            ind = i
            break
    times[:ind] = [first_valid] * ind


def correct_time(times):
    times = times
    for i in xrange(1, len(times)):
        diff = times[i] - times[i - 1]
        if diff < 0:
            times = times[:i] + list(array(times[i:]) - diff + 500)  # one TU step should be 500 ms
    return list(times)


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


class FitRes:
    def __init__(self, fit_obj=None):
        self.Pars = list(fit_obj.Parameters()) if (fit_obj is not None and len(fit_obj.Parameters()) > 0) else [None]
        self.Errors = list(fit_obj.Errors()) if (fit_obj is not None and len(fit_obj.Parameters()) > 0) else [None]
        self.NPars = len(self.Pars)
        self.Names = [fit_obj.ParName(i) for i in xrange(self.NPars)] if fit_obj is not None else [None]
        self.vChi2 = fit_obj.Chi2() if fit_obj is not None else None
        self.vNdf = fit_obj.Ndf() if fit_obj is not None else None

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


def do_nothing():
    pass