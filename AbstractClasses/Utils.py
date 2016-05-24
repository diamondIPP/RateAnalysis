#
# Created by Micha on 19.05.16.
#

from datetime import datetime
from termcolor import colored


# ==============================================
# UTILITY FUNCTIONS
# ==============================================

def log_warning(msg):
    t = datetime.now().strftime('%H:%M:%S')
    print '{head} {t} --> {msg}'.format(t=t, msg=msg, head=colored('WARNING:', 'red'))
