#
# Created by Micha on 19.05.16.
#

from datetime import datetime
from termcolor import colored
from ROOT import gStyle


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


def set_statbox(x=.98, y=.98, w=.16, entries=3):
    gStyle.SetStatX(x)
    gStyle.SetStatY(y)
    gStyle.SetStatW(w)
    gStyle.SetStatH(.05 * entries)


def do_nothing():
    pass
