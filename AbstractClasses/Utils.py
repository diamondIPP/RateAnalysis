# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on May 19th 2016 by M. Reichmann
# --------------------------------------------------------

from datetime import datetime, timedelta
from termcolor import colored
from ROOT import gStyle, gROOT, TF1
from numpy import sqrt
from os import makedirs
from os import path as pth
from time import time


# ==============================================
# UTILITY FUNCTIONS
# ==============================================
def log_warning(msg):
    t = datetime.now().strftime('%H:%M:%S')
    print '{head} {t} --> {msg}'.format(t=t, msg=msg, head=colored('WARNING:', 'red'))


def log_critical(msg):
    t = datetime.now().strftime('%H:%M:%S')
    print '{head} {t} --> {msg}'.format(t=t, msg=msg, head=colored('CRITICAL:', 'red'))
    quit()


def log_message(msg, overlay=False):
    t = datetime.now().strftime('%H:%M:%S')
    print '{ov}{t} --> {msg}{end}'.format(t=t, msg=msg, head=colored('WARNING:', 'red'), ov='\033[1A\r' if overlay else '', end=' ' * 20 if overlay else '')


def set_root_output(status=True):
    if status:
        gROOT.SetBatch(0)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
    else:
        gROOT.SetBatch(1)
        gROOT.ProcessLine("gErrorIgnoreLevel = kError;")


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


def set_statbox(x=.95, y=.88, w=.16, entries=3, only_fit=False, opt=None):
    if only_fit:
        gStyle.SetOptStat(0011)
        gStyle.SetOptFit(1)
    gStyle.SetOptStat(opt) if opt is not None else do_nothing()
    gStyle.SetStatX(x)
    gStyle.SetStatY(y)
    gStyle.SetStatW(w)
    gStyle.SetStatH(.04 * entries)


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


def joinpath(*args):
    return pth.join(*args)


def make_col_str(col):
    return '{0:2d}'.format(int(col)) if int(col) > 1 else '{0:3.1f}'.format(col)


def print_banner(msg, symbol='='):
    print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(msg)) * symbol, msg=msg)


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


def make_tc_str(tc, txt=True, data=False):
    if data:
        return datetime.strptime(tc, '%Y%m').strftime('psi_%Y_%m')
    elif tc[0].isdigit():
        return datetime.strptime(tc, '%Y%m').strftime('%B %Y' if txt else '%b%y')
    else:
        return datetime.strptime(tc, '%b%y').strftime('%Y%m' if txt else '%B %Y')


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


def set_drawing_range(h, legend=True, lfac=None, rfac=None):
    range_ = [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(10), h.FindLastBinAbove(10)]]
    lfac = lfac if lfac is not None else .2
    rfac = rfac if rfac is not None else .55 if legend else .1
    h.GetXaxis().SetRangeUser(*increased_range(range_, lfac, rfac))


def set_time_axis(histo, form='%H:%M', off=3600):
    histo.GetXaxis().SetTimeFormat(form)
    histo.GetXaxis().SetTimeOffset(-off)
    histo.GetXaxis().SetTimeDisplay(1)


def find_mpv_fwhm(histo, bins=15):
    max_bin = histo.GetMaximumBin()
    fit = TF1('fit', 'gaus', 0, 500)
    histo.Fit('fit', 'qs', '', histo.GetBinCenter(max_bin - bins), histo.GetBinCenter(max_bin + bins))
    mpv = fit.GetParameter(1)
    fwhm = histo.FindLastBinAbove(fit(mpv) / 2) - histo.FindFirstBinAbove(fit(mpv) / 2)
    return mpv, fwhm, mpv / fwhm


def make_cut_string(cut, n):
    return '{n}Cuts'.format(n=str(n).zfill(2)) if cut is not None else ''


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
