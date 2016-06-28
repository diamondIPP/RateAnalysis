#
# Created by Micha on 19.05.16.
#

from datetime import datetime
from termcolor import colored
from ROOT import gStyle
from math import sqrt


# ==============================================
# UTILITY FUNCTIONS
# ==============================================
def log_warning(msg):
    t = datetime.now().strftime('%H:%M:%S')
    print '{head} {t} --> {msg}'.format(t=t, msg=msg, head=colored('WARNING:', 'red'))


def log_message(msg):
    t = datetime.now().strftime('%H:%M:%S')
    print '{t} --> {msg}'.format(t=t, msg=msg, head=colored('WARNING:', 'red'))


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


def set_statbox(x=.96, y=.96, w=.16, entries=3, only_fit=False):
    if only_fit:
        gStyle.SetOptStat(0011)
        gStyle.SetOptFit(1)
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
    return num / val * val


def round_up_to(num, val):
    return num / val * val + val


def scale_multigraph(mg, val=1):
    first = mg.GetListOfGraphs()[0].GetY()[0]
    for gr in mg.GetListOfGraphs():
        for i in xrange(gr.GetN()):
            gr.SetPoint(i, gr.GetX()[i], gr.GetY()[i] / first * val)
            try:
                gr.SetPointError(i, gr.GetErrorX(i), gr.GetErrorY(i) / first * val)
            except Exception as err:
                print err
                pass


def make_transparent(pad):
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.SetFrameFillStyle(4000)


def hide_axis(axis):
    axis.SetTickLength(0)
    axis.SetLabelOffset(99)
    axis.SetTitleOffset(99)


def move_legend(l, x1, y1):
    xdiff = l.GetX2NDC() - l.GetX1NDC()
    ydiff = l.GetY2NDC() - l.GetY1NDC()
    l.SetX1NDC(x1)
    l.SetX2NDC(x1 + xdiff)
    l.SetY1NDC(y1)
    l.SetY2NDC(y1 + ydiff)


def calc_mean(l):
    mean = sum(l) / len(l)
    mean2 = sum(map(lambda x: x**2, l)) / len(l)
    sigma = sqrt(mean2 - mean**2)
    return mean, sigma


def calc_weighted_mean(means, sigmas):
    weights = map(lambda x: x**(-2), sigmas)
    variance = 1 / sum(weights)
    mean = sum(map(lambda x, y: x * y, means, weights))
    return mean * variance, sqrt(variance)


def print_banner(msg, symbol='='):
    print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(msg)) * symbol, msg=msg)


def print_small_banner(msg, symbol='-'):
    print '\n{delim}\n{msg}\n'.format(delim=len(str(msg)) * symbol, msg=msg)


def do_nothing():
    pass
