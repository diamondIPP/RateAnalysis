# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on May 19th 2016 by M. Reichmann
# --------------------------------------------------------

from datetime import datetime, timedelta
from termcolor import colored
from ROOT import gStyle, gROOT, TF1, TColor, TFile
from numpy import sqrt, array, average, mean, arange
from os import makedirs
from os import path as pth
from os.path import basename, join, split
from time import time, sleep, mktime
from collections import OrderedDict
import pickle
from threading import Thread
from multiprocessing import Pool


# ==============================================
# UTILITY FUNCTIONS
# ==============================================
def get_t_str():
    return datetime.now().strftime('%H:%M:%S')


def log_warning(msg):
    print '{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head=colored('WARNING:', 'red'))


def log_critical(msg):
    print '{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head=colored('CRITICAL:', 'red'))
    raise Exception


def log_message(msg, overlay=False, prnt=True):
    if prnt:
        print '{ov}{t} --> {msg}{end}'.format(t=get_t_str(), msg=msg, ov='\033[1A\r' if overlay else '', end=' ' * 20 if overlay else '')


def set_root_output(status=True):
    gROOT.SetBatch(not status)
    gROOT.ProcessLine('gErrorIgnoreLevel = {e};'.format(e='0' if status else 'kError'))


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


def set_statbox(x=.95, y=.88, w=.16, entries=3, only_fit=False, only_entries=False, opt=None, form=None):
    if only_fit:
        gStyle.SetOptStat(0011)
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


def scale_multigraph(mg, val=1):
    g = mg.GetListOfGraphs()[0]
    points = {g.GetX()[i]: g.GetY()[i] for i in xrange(g.GetN())}
    y = points[min(points)]
    for gr in mg.GetListOfGraphs():
        for i in xrange(gr.GetN()):
            gr.SetPoint(i, gr.GetX()[i], gr.GetY()[i] / y * val)
            try:
                gr.SetPointError(i, gr.GetErrorX(i), gr.GetErrorY(i) / y * val)
            except Exception as err:
                log_warning('Error in scale multigraph: {err}'.format(err=err))


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
    if not float(val):
        return 'unirradiated'
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


def weighted_avrg_std(values, weights):
    """ Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape. """
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


def file_exists(path):
    return pth.isfile(path)


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
    histo.Fit('fit', 'qs', '', histo.GetBinCenter(max_bin - bins), histo.GetBinCenter(max_bin + bins))
    mpv = fit.GetParameter(1)
    fwhm = histo.FindLastBinAbove(fit(mpv) / 2) - histo.FindFirstBinAbove(fit(mpv) / 2)
    return mpv, fwhm, mpv / fwhm


def get_fwhm(h):
    max_v = h.GetMaximum()
    return h.GetBinCenter(h.FindLastBinAbove(max_v / 2)) - h.GetBinCenter(h.FindFirstBinAbove(max_v / 2))


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
    closing_row = '~' * len('| {r} |'.format(r=' | '.join(rows[-1])))
    if header is not None:
        print closing_row
        print '| {r} |'.format(r=' | '.join(header))
    print closing_row
    for row in rows:
        print '| {r} |'.format(r=(' | '.join(row) if len(row) > 1 else row[0]).ljust(len(closing_row) - 4))
    print closing_row


def get_base_dir():
    return join('/', *__file__.split('/')[1:3])


def server_is_mounted():
    return dir_exists(join(get_base_dir(), 'mounts/psi/Diamonds'))


def do_pickle(path, func, value=None, params=None, redo=False):
    if value is not None:
        f = open(path, 'w')
        pickle.dump(value, f)
        f.close()
        return value
    if file_exists(path) and not redo:
        f = open(path, 'r')
        return pickle.load(f)
    else:
        ret_val = func() if params is None else func(params)
        f = open(path, 'w')
        pickle.dump(ret_val, f)
        f.close()
        return ret_val


def kinder_pickle(old_path, value):
    if server_is_mounted():
        picklepath = join(get_base_dir(), 'mounts/psi/Pickles', basename(split(old_path)[0]), basename(old_path))
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


def scale_axis(xmin, xmax, ymin, ymax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetXaxis().SetRangeUser(xmin, xmax)
    h.GetYaxis().SetRangeUser(ymin, ymax)


def remove_letters(string):
    new_str = ''
    for l in string:
        if l.isdigit():
            new_str += l
    return new_str


def get_last_canvas():
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        log_warning('There is no canvas is in the list...')


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


class FitRes:
    def __init__(self, fit_obj=None):
        self.Pars = list(fit_obj.Parameters()) if (fit_obj is not None and len(fit_obj.Parameters()) > 0) else [None]
        self.Errors = list(fit_obj.Errors()) if (fit_obj is not None and len(fit_obj.Parameters()) > 0) else [None]

    def Parameter(self, arg):
        return self.Pars[arg]

    def ParError(self, arg):
        return self.Errors[arg]


def do_nothing():
    pass
