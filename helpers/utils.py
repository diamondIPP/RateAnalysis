# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on May 19th 2016 by M. Reichmann
# --------------------------------------------------------
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TF1, TSpectrum, TTree

import pickle
from collections import OrderedDict
from copy import deepcopy
from datetime import datetime, timedelta
from multiprocessing import Pool, cpu_count
from subprocess import call
from threading import Thread
from time import time, sleep

from gtts import gTTS
from numpy import sqrt, array, average, mean, arange, log10, concatenate, where, count_nonzero, full, ndarray, exp, sin, cos, arctan, zeros, dot, roll, arctan2, frombuffer, split, cumsum
from numpy import histogram, log2, diff, isfinite, pi, corrcoef, quantile, all
from os import makedirs, remove, devnull, stat, getenv, _exit
from os import path as pth
from os.path import dirname, realpath, join
from pytz import timezone, utc
from termcolor import colored
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.core import Variable, AffineScalarFunc
from argparse import ArgumentParser
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar, SimpleProgress, Widget
from configparser import ConfigParser, NoSectionError, NoOptionError
from scipy.optimize import curve_fit
from scipy import constants
import h5py
from functools import partial, wraps
from queue import Queue
from json import load, loads
from inspect import signature
from types import FunctionType, MethodType
from glob import glob

OFF = False
ON = True
DEGREE_SIGN = u'\N{DEGREE SIGN}'
COUNT = 0
Dir = dirname(dirname(realpath(__file__)))

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


def warning(msg, prnt=True):
    if prnt:
        print(prepare_msg(msg, 'WARNING', 'yellow'))


def critical(msg):
    print(prepare_msg(msg, 'CRITICAL', 'red'))
    _exit(5)


def prepare_msg(msg, head, color=None, attrs=None, blank_lines=0):
    return '{}\r{} {} --> {}'.format('\n' * blank_lines, colored(head, color, attrs=choose(make_list(attrs), None, attrs)), get_t_str(), msg)


def info(msg, endl=True, blank_lines=0, prnt=True):
    if prnt:
        print(prepare_msg(msg, 'INFO', 'cyan', 'dark', blank_lines), flush=True, end='\n' if endl else ' ')
    return time()


def add_to_info(t, msg='Done', prnt=True):
    if prnt:
        print('{m} ({t:2.2f} s)'.format(m=msg, t=time() - t))


def print_check():
    global COUNT
    print('======={}========'.format(COUNT))
    COUNT += 1


def untitle(string):
    s = ''
    for word in string.split(' '):
        if word:
            s += word[0].lower() + word[1:] + ' '
    return s.strip(' ')


def round_down_to(num, val=1):
    return int(num) // val * val


def round_up_to(num, val=1):
    return int(num) // val * val + val


def interpolate_two_points(x1, y1, x2, y2, name=''):
    # f = p1*x + p0
    p1 = (y1 - y2) / (x1 - x2)
    p0 = y1 - x1 * p1
    w = abs(x2 - x1)
    fit_range = array(sorted([x1, x2])) + [-w / 3., w / 3.]
    f = TF1('fpol1{}'.format(name), 'pol1', *fit_range)
    f.SetParameters(p0, p1)
    return f


def get_x(x1, x2, y1, y2, y):
    return (x2 - x1) / (y2 - y1) * (y - y1) + x1


def get_y(x1, x2, y1, y2, x):
    return get_x(y1, y2, x1, x2, x)


def interpolate_x(x1, x2, y1, y2, y):
    p1 = get_p1(x1, x2, y1, y2)
    p0 = get_p0(x1, y1, p1)
    return (y - p0) / p1 if p1 else 0


def interpolate_y(x1, x2, y1, y2, x):
    p1 = get_p1(x1, x2, y1, y2)
    p0 = get_p0(x1, y1, p1)
    return p1 * x + p0


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


def make_irr_string(val):
    if '?' in val:
        return val
    if not float(val):
        return 'nonirradiated'
    val, power = [float(i) for i in val.split('e')]
    return '{v:1.1f}#upoint10^{p} n/cm^{{2}}'.format(v=val, p='{{{}}}'.format(int(power)))


def mean_sigma(values, weights=None, err=True):
    """ Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape. """
    if len(values) == 1:
        value = make_ufloat(values[0])
        return (value, ufloat(value.s, 0)) if err else (value.n, value.s)
    weights = full(len(values), 1) if weights is None else weights
    if is_ufloat(values[0]):
        errors = array([v.s for v in values])
        weights = full(errors.size, 1) if all(errors == errors[0]) else [1 / e if e else 0 for e in errors]
        values = array([v.n for v in values], 'd')
    if all(weights == 0):
        return [0, 0]
    n, avrg = values.size, average(values, weights=weights)
    sigma = sqrt(n / (n - 1) * average((values - avrg) ** 2, weights=weights))  # Fast and numerically precise
    m, s = ufloat(avrg, sigma / (sqrt(len(values)) - 1)), ufloat(sigma, sigma / sqrt(2 * len(values)))
    return (m, s) if err else (m.n, s.n)


def binned_stats(x, values, f, bins):
    return array([f(v) for v in split(values, cumsum(histogram(x, bins)[0].astype('i'))[:-1])])


def make_latex_table_row(row, hline=False):
    row_str, end = ' & '.join(row), '\\hline' if hline else ''
    return f'  {row_str} \\\\{end}'


def make_latex_table(header, rows, hlines=False):
    header, cols = [f'\\textbf{{{head}}}' if head else '' for head in header], array(rows, str).T
    max_width = [len(max(col, key=len)) for col in cols]
    rows = array([[f'{word:<{w}}' for word in col] for col, w in zip(cols, max_width)]).T
    rows = '\n'.join(make_latex_table_row(row, hlines) for row in rows)
    return f'{make_latex_table_row(header, hline=True)}\n{rows}'


def si(v, f='.1f', unit=None):
    unit = '' if unit is None else f'\\{unit}'
    return f'\\SI{{{v:{f}}}}{{{unit}}}'.replace('/', '')


def make_dia_str(dia):
    dia = dia.replace('-', '')
    return '{0}{1}'.format(dia[0].title(), dia[1:])


def make_list(value):
    return array([value], dtype=object).flatten()


def uarr2n(arr):
    return array([i.n for i in arr]) if len(arr) and is_ufloat(arr[0]) else arr


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


def file_is_beeing_written(file_path):
    file_size = stat(file_path)
    sleep(4)
    return file_size != stat(file_path)


def make_col_str(col):
    return '{0:2d}'.format(int(col)) if int(col) > 1 else '{0:3.1f}'.format(col)


def print_banner(msg, symbol='~', new_lines=1, color=None):
    msg = '{} |'.format(msg)
    print(colored('{n}{delim}\n{msg}\n{delim}{n}'.format(delim=len(str(msg)) * symbol, msg=msg, n='\n' * new_lines), color))


def print_small_banner(msg, symbol='-', color=None):
    print(colored('\n{delim}\n{msg}\n'.format(delim=len(str(msg)) * symbol, msg=msg), color))


def print_elapsed_time(start, what='This', show=True, color=None):
    string = 'Elapsed time for {w}: {t}'.format(t=get_elapsed_time(start), w=what)
    print_banner(string, color=color) if show else do_nothing()
    return string


def make_byte_string(v):
    n = int(log2(v) // 10) if v else 0
    return '{:1.1f} {}'.format(v / 2 ** (10 * n), ['B', 'kB', 'MB', 'GB'][n])


def get_elapsed_time(start, hrs=False):
    t = str(timedelta(seconds=round(time() - start, 0 if hrs else 2)))
    return t if hrs else t[2:-4]


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


def make_flux_string(rate, prec=1, term=False):
    unit = f'{"MHz" if rate > 1000 else "kHz"}/cm{"²" if term else "^{2}"}'
    return f'{rate / (1000 if rate > 1000 else 1):2.{prec if rate > 1000 else 0}f} {unit}'


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


def cart2pol(x, y):
    return array([sqrt(x ** 2 + y ** 2), arctan2(y, x)])


def pol2cart(rho, phi):
    return array([rho * cos(phi), rho * sin(phi)])


def make_cut_string(cut, n):
    return '{n}Cuts'.format(n=str(n).zfill(2)) if cut is not None else ''


def get_resolution():
    try:
        from screeninfo import get_monitors
        m = get_monitors()
        return round_down_to(m[0].height, 500)
    except Exception as exc:
        warning('Could not get resolution! Using default ...\n\t{e}'.format(e=exc))
        return 1000


def print_table(rows, header=None, footer=None, prnt=True):
    head, foot = [choose([v], zeros((0, len(rows[0]))), v) for v in [header, footer]]
    t = concatenate([head, rows, foot]).astype('str')
    col_width = [len(max(t[:, i], key=len)) for i in range(t.shape[1])]
    total_width = sum(col_width) + len(col_width) * 3 + 1
    hline = '{}'.format('~' * total_width)
    lines = []
    for i, row in enumerate(t):
        if i in [0] + choose([1], [], header) + choose([t.shape[0] - 1], [], footer):
            lines.append(hline)
        lines.append('| {r} |'.format(r=' | '.join(word.ljust(n) for word, n in zip(row, col_width))))
    lines.append('{}\n'.format(hline))
    if prnt:
        print('\n'.join(lines))
    return '\n'.join(lines)


def get_base_dir():
    return dirname(dirname(realpath(__file__)))


def make_meta_path(main_dir, sub_dir='', name='', ext='pickle', suffix=''):
    ensure_dir(join(main_dir, sub_dir))
    suf = '{}{}'.format('-' if suffix and name else '', '_'.join(make_list(suffix)))
    return join(main_dir, sub_dir, '{}{}.{}'.format(name, suf, ext.strip('.')))


def do_pickle(path, func, value=None, redo=False, *args, **kwargs):
    if value is not None:
        with open(path, 'wb') as f:
            pickle.dump(value, f)
        return value
    try:
        if file_exists(path) and not redo:
            with open(path, 'rb') as f:
                return pickle.load(f)
    except ImportError:
        pass
    ret_val = func(redo=redo, *args, **kwargs) if type(func) in [MethodType, FunctionType] and 'redo' in signature(func).parameters else func(*args, **kwargs)
    with open(path, 'wb') as f:
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


def find_maxima(h, n=3, sigma=2, sort_x=False):
    s = TSpectrum(n)
    n = s.Search(h, sigma)  # return how many maxima were found
    v = array([frombuffer(getattr(s, 'GetPosition{}'.format(i))(), dtype='d', count=n) for i in ['X', 'Y']]).T
    return v[v[:, 0].argsort()] if sort_x else v


def int_to_roman(integer):
    """ Convert an integer to Roman numerals. """
    if type(integer) != int:
        raise TypeError(f'cannot convert {type(integer).__name__} to roman, integer required')
    if not 0 < integer < 4000:
        raise ValueError('Argument must be between 1 and 3999')
    dic = OrderedDict([(1000, 'M'), (900, 'CM'), (500, 'D'), (400, 'CD'), (100, 'C'), (90, 'XC'), (50, 'L'), (40, 'XL'),
                       (10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')])
    result = ''
    for i, num in dic.items():
        count = int(integer // i)
        result += num * count
        integer -= i * count
    return result


def add_spaces(s):
    return ''.join(f' {s[i]}' if i and (s[i].isupper() or s[i].isdigit()) and not s[i - 1].isdigit() and not s[i - 1].isupper() else s[i] for i in range(len(s)))


def remove_letters(s):
    return ''.join(filter(str.isdigit, s))


def remove_digits(string):
    return ''.join(x for x in string if not x.isdigit())


def do(fs, pars, exe=-1):
    fs, pars = ([fs], [pars]) if type(fs) is not list else (fs, pars)
    exe = pars if exe == -1 else [exe]
    for f, p, e in zip(fs, pars, exe):
        f(p) if e is not None else do_nothing()


def choose(v, default, decider='None', *args, **kwargs):
    use_default = decider is None if decider != 'None' else v is None
    if callable(default) and use_default:
        default = default(*args, **kwargs)
    return default if use_default else v(*args, **kwargs) if callable(v) else v


def average_list(lst, n):
    return [mean(lst[i:i + n]) for i in arange(0, len(lst), n)] if n > 1 else lst


def log_bins(n_bins, min_val, max_val):
    width = (log10(max_val) - log10(min_val)) / float(n_bins)
    return [n_bins, array([pow(10, log10(min_val) + i * width) for i in range(n_bins + 1)])]


def fit2u(fit, par):
    return ufloat(fit.Parameter(par), fit.ParError(par))


def eff2u(eff):
    return ufloat(eff[0], mean(eff[1:]))


def make_ufloat(n, s=0):
    return array([ufloat(*v) for v in array([n, s]).T]) if is_iter(n) else n if is_ufloat(n) else ufloat(n, s)


def is_ufloat(value):
    return type(value) in [Variable, AffineScalarFunc]


def is_iter(v):
    try:
        iter(v)
        return True
    except TypeError:
        return False


def get_time_vec(sel, run=None):
    if hasattr(run, 'Number'):
        tree = run.load_rootfile()
        run = run.Number
    elif type(sel) == TTree:
        tree = sel
    else:
        critical('unknown type')
        return
    if tree is None:
        return
    tree.SetEstimate(-1)
    tvec = get_tree_vec(tree, 'time')
    fill_empty_time_entries(tvec)
    time_vec = correct_time(tvec, run)
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
        delta = times[i_off + 1] - times[i_off]
        i_off += 1  # because range does not include the last index
        return concatenate((times[:i_off], times[i_off:] - delta + 500))  # one TU step should be 500 ms
    return times


def say(txt, lang='en'):
    tts = gTTS(text=txt, lang=lang)
    tts.save('good.mp3')
    with open(devnull, 'w') as FNULL:
        call(([] if getenv('SSH_TTY') is None else ['DISPLAY=:0']) + ['mpg321', 'good.mp3'], stdout=FNULL)
    remove('good.mp3')


def remove_file(file_path, prnt=True):
    if file_exists(file_path):
        warning('removing {}'.format(file_path), prnt=prnt)
        remove(file_path)


def remove_files(paths, prnt=True, wildcard=False):
    for name in (glob(paths) if wildcard else paths):
        remove_file(name, prnt)


def get_running_time(t):
    now = datetime.fromtimestamp(time() - t) - timedelta(hours=1)
    return now.strftime('%H:%M:%S')


def get_tree_vec(tree, var, cut='', dtype=None, nentries=None, firstentry=0):
    strings = make_list(var)
    n = tree.Draw(':'.join(strings), cut, 'goff', choose(nentries, tree.kMaxEntries), firstentry)
    dtypes = dtype if type(dtype) in [list, ndarray] else full(strings.size, dtype)
    vals = [get_buf(tree.GetVal(i), n, dtypes[i]) for i in range(strings.size)]
    return vals[0] if len(vals) == 1 else vals


def get_arg(arg, default):
    return default if arg is None else arg


def get_buf(buf, n, dtype=None):
    return frombuffer(buf, dtype=buf.typecode, count=n).astype(dtype)


def calc_eff(k=0, n=0, values=None):
    values = array(values) if values is not None else None
    if n == 0 and (values is None or not values.size):
        return zeros(3)
    k = float(k if values is None else count_nonzero(values))
    n = float(n if values is None else values.size)
    m = (k + 1) / (n + 2)
    mode = k / n
    s = sqrt(((k + 1) / (n + 2) * (k + 2) / (n + 3) - ((k + 1) ** 2) / ((n + 2) ** 2)))
    return array([mode, max(s + (mode - m), 0), max(s - (mode - m), 0)]) * 100


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
    t = info('Measuring time of method {}:'.format(f.__name__), endl=False)
    for _ in range(int(rep)):
        f(*args, **kwargs)
    add_to_info(t, '')


def u_to_str(v, prec=2):
    return '{{:1.{0}f}} ({{:1.{0}f}})'.format(prec).format(v.n, v.s)


def poly_area(x, y):
    return .5 * abs(dot(x, roll(y, 1)) - dot(y, roll(x, 1)))


def discrete_int(x, y):
    """ assume linear interpolation between the points. """
    cut = x.argsort()
    x, y = x[cut], y[cut]
    dx, dy = diff(x), diff(y)
    i = dx * y[:-1] + .5 * dx * dy
    return sum(i[isfinite(i)])


def kramers_kronig(x, y):
    return 1 + 2 / pi * array([discrete_int(x, x * y / (x ** 2 - ix ** 2)) for ix in x])


def freedman_diaconis(x):
    return 2 * (quantile(x, .75) - quantile(x, .25)) / x.size ** (1 / 3)


def p2ecut(n, cut):
    c = zeros(n.size, '?')
    c[arange(n.size).repeat(n)[cut]] = True
    return c


def correlate(l1, l2):
    if len(l1.shape) == 2:
        x, y = l1.flatten(), l2.flatten()
        cut, s = (x > 0) & (y > 0), count_nonzero(x)
        return correlate(x[cut], y[cut]) if count_nonzero(cut) > .6 * s else 0
    return corrcoef(l1, l2)[0][1]


def prep_kw(dic, **default):
    d = deepcopy(dic)
    for kw, value in default.items():
        if kw not in d:
            d[kw] = value
    return d


def make_suffix(ana, values):
    suf_vals = [ana.get_short_name(suf) if type(suf) is str and suf.startswith('TimeIntegralValues') else suf for suf in values]
    return '_'.join(str(int(val) if isint(val) else val.GetName() if hasattr(val, 'GetName') else val) for val in suf_vals if val is not None)


def prep_suffix(f, args, kwargs, suf_args, field=None):
    def_pars = signature(f).parameters
    names, values = list(def_pars.keys()), [par.default for par in def_pars.values()]
    i_arg = (arange(len([n for n in names if n not in ['self', '_redo']])) if suf_args == 'all' else make_list(loads(str(suf_args)))) + 1
    suf_vals = [args[i] if len(args) > i else kwargs[names[i]] if names[i] in kwargs else values[i] for i in i_arg]
    suf_vals += [getattr(args[0], str(field))] if field is not None and hasattr(args[0], field) else []
    return make_suffix(args[0], suf_vals)


def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)


def save_pickle(*pargs, print_dur=False, low_rate=False, high_rate=False, suf_args='[]', field=None, **pkwargs):
    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if '_no_save' in kwargs:
                return func(*args, **kwargs)
            run = args[0].Run.get_high_rate_run(high=not low_rate) if low_rate or high_rate else None
            pickle_path = args[0].make_simple_pickle_path(*pargs, **prep_kw(pkwargs, run=run, suf=prep_suffix(func, args, kwargs, suf_args, field)))
            redo = (kwargs['_redo'] if '_redo' in kwargs else False) or (kwargs['show'] if 'show' in kwargs else False)
            if file_exists(pickle_path) and not redo:
                return load_pickle(pickle_path)
            prnt = print_dur and (kwargs['prnt'] if 'prnt' in kwargs else True)
            t = (args[0].info if hasattr(args[0], 'info') else info)(f'{args[0].__class__.__name__}: {func.__name__.replace("_", " ")} ...', endl=False, prnt=prnt)
            value = func(*args, **kwargs)
            with open(pickle_path, 'wb') as f:
                pickle.dump(value, f)
            (args[0].add_to_info if hasattr(args[0], 'add_to_info') else add_to_info)(t, prnt=prnt)
            return value
        return wrapper
    return inner


def save_hdf5(*pargs, suf_args='[]', **pkwargs):
    def inner(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            file_path = args[0].make_simple_hdf5_path(*pargs, **prep_kw(pkwargs, suf=prep_suffix(f, args, kwargs, suf_args)))
            redo = kwargs['_redo'] if '_redo' in kwargs else False
            if file_exists(file_path) and not redo:
                return h5py.File(file_path, 'r')['data']
            remove_file(file_path)
            data = f(*args, **kwargs)
            hf = h5py.File(file_path, 'w')
            hf.create_dataset('data', data=data)
            return hf['data']
        return wrapper
    return inner


def print_duration(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        prnt = (kwargs['prnt'] if 'prnt' in kwargs else True)
        t = (args[0].info if hasattr(args[0], 'info') else info)(f'{func.__name__.replace("_", " ")} ...', endl=False, prnt=prnt)
        value = func(*args, **kwargs)
        (args[0].add_to_info if hasattr(args[0], 'add_to_info') else add_to_info)(t, prnt=prnt)
        return value
    return wrapper


def quiet(func):
    @wraps(func)
    def wrapper(analysis, *args, **kwargs):
        ana = analysis.Ana if hasattr(analysis, 'Ana') else analysis
        old = ana.Verbose
        ana.set_verbose(False)
        value = func(analysis, *args, **kwargs)
        ana.set_verbose(old)
        return value
    return wrapper


# ----------------------------------------
# region CLASSES


def update_pbar(func):
    @wraps(func)
    def my_func(*args, **kwargs):
        value = func(*args, **kwargs)
        if args[0].PBar is not None and args[0].PBar.PBar is not None and not args[0].PBar.is_finished():
            args[0].PBar.update()
        return value
    return my_func


class PBar(object):
    def __init__(self, start=None, counter=False, t=None):
        self.PBar = None
        self.Widgets = self.init_widgets(counter, t)
        self.Step = 0
        self.N = 0
        self.start(start)

    def __reduce__(self):
        return self.__class__, (None, False, None), (self.Widgets, self.Step, self.N)

    def __setstate__(self, state):
        self.Widgets, self.Step, self.N = state
        if self.N:
            self.PBar = ProgressBar(widgets=self.Widgets, maxval=self.N).start()
            self.update(self.Step) if self.Step > 0 else do_nothing()

    @staticmethod
    def init_widgets(counter, t):
        return ['Progress: ', SimpleProgress('/') if counter else Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed() if t is None else EventSpeed(t)]

    def start(self, n, counter=None, t=None):
        if n is not None:
            self.Step = 0
            self.PBar = ProgressBar(widgets=self.Widgets if t is None and counter is None else self.init_widgets(counter, t), maxval=n).start()
            self.N = n

    def update(self, i=None):
        i = self.Step if i is None else i
        if i >= self.PBar.maxval:
            return
        self.PBar.update(i + 1)
        self.Step += 1
        if i == self.PBar.maxval - 1:
            self.finish()

    def set_last(self):
        if self.PBar:
            self.PBar.currval = self.N
            self.PBar.finished = True

    def finish(self):
        self.PBar.finish()

    def is_finished(self):
        return self.PBar.currval == self.N


class EventSpeed(Widget):
    """Widget for showing the event speed (useful for slow updates)."""

    def __init__(self, t='s'):
        self.unit = t
        self.factor = {'s': 1, 'min': 60, 'h': 60 * 60}[t]

    def update(self, pbar):
        value = 0
        if pbar.seconds_elapsed > 2e-6 and pbar.currval > 2e-6:
            value = pbar.currval / pbar.seconds_elapsed * self.factor
        return f'{value:4.1f} E/{self.unit}'


def load_main_config(config='main', ext='ini'):
    file_name = join(get_base_dir(), 'config', '{}.{}'.format(config.split('.')[0], ext.strip('.')))
    if not file_exists(file_name):
        critical('{} does not exist. Please copy it from the main.default and adapt it to your purpose!'.format(file_name))
    return Config(file_name)


class Config(ConfigParser):

    def __init__(self, file_name, **kwargs):
        super(Config, self).__init__(**kwargs)
        self.FileName = file_name
        self.read(file_name) if type(file_name) is not list else self.read_file(file_name)

    def get_value(self, section, option, dtype: type = str, default=None):
        dtype = type(default) if default is not None else dtype
        try:
            if dtype is bool:
                return self.getboolean(section, option)
            v = self.get(section, option)
            return loads(v) if dtype == list or '[' in v and dtype is not str else dtype(v)
        except (NoOptionError, NoSectionError):
            return default

    def get_values(self, section):
        return [j for i, j in self.items(section)]

    def get_list(self, section, option, default=None):
        return self.get_value(section, option, list, choose(default, []))

    def get_ufloat(self, section, option, default=None):
        return ufloat_fromstr(self.get_value(section, option, default=default))

    def show(self):
        for key, section in self.items():
            print(colored('[{}]'.format(key), 'yellow'))
            for option in section:
                print('{} = {}'.format(option, self.get(key, option)))
            print()

    def write(self, file_name=None, space_around_delimiters=True):
        with open(choose(file_name, self.FileName), 'w') as f:
            super(Config, self).write(f, space_around_delimiters)
# endregion CLASSES
# ----------------------------------------


# ----------------------------------------
# region RELATIVITY
def calc_speed(p, m):
    return 1 / sqrt(1 + m * m / (p * p))


def beta_gamma(p, m):
    v = calc_speed(p, m)
    return lorentz_factor(v) * v


def beta(bg):
    return sqrt(1 / (1 / (bg * bg) + 1))


def gamma(bg):
    return bg / beta(bg)


def t_diff(s, p, m1, m2):
    return s * (1 / calc_speed(p, m1) - 1 / calc_speed(p, m2)) / constants.c * 1e9


def e_kin(p, m):
    return sqrt(p**2 + m**2) - m


def p2e(p, m):
    return e_kin(p, m)


def e2p(e, m):
    return sqrt((e + m) * (e + m) - m * m)


def lorentz_factor(v):
    return 1 / sqrt(1 - v * v)


def momentum(m, v):
    return m * v * lorentz_factor(v)


def decay_ratio(p, m, d, tau):
    return exp(-d * m / (tau * 1e-9 * p * constants.c))


def decay_momentum(m, m1, m2=0):
    return sqrt((m**2 + m1**2 + m2**2)**2 - 4 * m**2 * m1**2) / (2 * m)


def decay_energy(m, m1, m2=0):
    return (m**2 + m1**2 - m2**2) / (2 * m)


def decay_angle(theta, p, m, m1, m2=0):
    p1 = decay_momentum(m, m1, m2)
    v = calc_speed(p, m)
    return arctan(p1 * sin(theta) / (lorentz_factor(v) * (p1 * cos(theta) + v * decay_energy(m, m1, m2))))

# endregion RELATIVITY
# ----------------------------------------


def gauss(x, scale, mean_, sigma, off=0):
    return scale * exp(-.5 * ((x - mean_) / sigma) ** 2) + off


def fit_data(f, y, x=None, p=None):
    x = arange(y.shape[0]) if x is None else x
    return curve_fit(f, x, y, p0=p)


def multi_threading(lst, timeout=60 * 60 * 2):
    """ runs several threads in parallel. [lst] must contain tuples of the methods and the arguments as list."""
    t0 = info('Run multithreading on {} tasks ... '.format(len(lst)), endl=False)
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
    t = info('Run parallelisation on {} tasks ... '.format(len(args_list)), endl=False)
    pool = Pool(cpu_count())
    workers = [pool.apply_async(f, make_list(args)) for args in args_list]
    results = [worker.get(timeout) for worker in workers]
    add_to_info(t)
    return results


def parallelise_instance(instances, method, args, timeout=60 * 60):
    t = info('Run parallelisation on {} tasks ... '.format(len(args)), endl=False)
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


def get_input(msg, default='None'):
    txt = input(f'{msg} (press enter for default: {default}): ')
    return txt if txt else default


def plural(word, pluralise=True):
    return f'{word}s' if pluralise else word


def do_nothing():
    pass
