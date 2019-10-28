from json import loads

from ROOT import TCut, gROOT, TH1F, TPie
from numpy import zeros, histogram2d, split, histogram

from InfoLegend import InfoLegend
from utils import *
from draw import format_pie
from binning import Bins


class Cut:
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, parent_analysis, skip=False):

        if not skip:
            self.Analysis = parent_analysis
            self.RunNumber = self.Analysis.RunNumber
            self.InfoLegend = InfoLegend(parent_analysis)
            self.Bins = Bins(self.Analysis.Run)

            # config
            self.Config = self.Analysis.Config
            self.CutConfig = {}
            self.NCuts = 0
            self.LowRateRun = None
            self.HighRateRun = None

            # define cut strings
            self.EasyCutStrings = self.init_easy_cutstrings()
            self.CutStrings = self.define_cutstrings()

            self.JumpCut = TCut('JumpCut', '')

            # beam interrupts
            self.Jumps = None
            self.Interruptions = None

            self.load_config()
            # generate cut strings
            self.generate_cut_string()
            self.AllCut = self.generate_all_cut()

    def __call__(self, *args, **kwargs):
        return self.generate_all_cut()

    def set_high_low_rate_run(self, high_run, low_run):
        self.LowRateRun = str(low_run)
        self.HighRateRun = str(high_run)

    def generate_special_cut(self, excluded=None, included=None, name='special_cut', prnt=True):
        cut = TCut(name, '')
        self.NCuts = 0
        for key, value in self.CutStrings.iteritems():
            if excluded and key in excluded:
                continue
            if included and key not in included:
                continue
            if key.startswith('old') or key.startswith('AllCut'):
                continue
            if value.GetTitle() == '':
                continue
            cut += value
            self.NCuts += 1
        self.Analysis.info('generated {name} cut with {num} cuts'.format(name=name, num=self.NCuts)) if prnt else do_nothing()
        return cut

    def generate_all_cut(self):
        cut = TCut('AllCuts', '')
        self.NCuts = 0
        for key, value in self.CutStrings.iteritems():
            if not key.startswith('old') and not key.startswith('AllCut'):
                cut += value
                self.NCuts += 1
        return cut

    @staticmethod
    def init_easy_cutstrings():
        dic = OrderedDict()
        dic['IndividualChCut'] = ''
        dic['EventRange'] = ''
        dic['noPulser'] = ''
        dic['ExcludeFirst'] = ''
        dic['notSaturated'] = ''
        dic['noBeamInter'] = ''
        dic['Tracks'] = ''
        dic['peakPos_high'] = ''
        dic['pedestalsigma'] = ''
        dic['alignment'] = ''
        return dic
    
    @staticmethod
    def define_cutstrings():
        """ Defines the ordered dictionary that contains all the final cuts and the order they are going to be applied."""
        dic = OrderedDict()
        dic['raw'] = TCut('raw', '')
        # waveform
        dic['saturated'] = TCut('saturated', '')
        dic['pulser'] = TCut('pulser', '')
        # general
        dic['event_range'] = TCut('event_range', '')
        dic['beam_interruptions'] = TCut('beam_interruptions', '')
        dic['aligned'] = TCut('aligned', '')
        dic['trigger_phase'] = TCut('trigger_phase', '')
        dic['tracks'] = TCut('tracks', '')
        # waveform
        dic['ped_sigma'] = TCut('ped_sigma', '')
        dic['median'] = TCut('median', '')
        dic['threshold'] = TCut('threshold', '')
        dic['signal_peak_pos'] = TCut('signal_peak_pos', '')
        dic['signal_peak_time'] = TCut('signal_peak_time', '')
        dic['trigger_cell'] = TCut('trigger_cell', '')
        dic['timing'] = TCut('timing', '')
        dic['old_bucket'] = TCut('old_bucket', '')
        dic['bucket'] = TCut('bucket', '')
        # tracks
        dic['hit'] = TCut('hit', '')
        dic['masks'] = TCut('masks', '')
        dic['chi2X'] = TCut('chi2X', '')
        dic['chi2Y'] = TCut('chi2Y', '')
        dic['slope_x'] = TCut('slope_x', '')
        dic['slope_y'] = TCut('slope_y', '')
        dic['rhit'] = TCut('rhit', '')
        dic['fiducial'] = TCut('fiducial', '')
        dic['AllCuts'] = TCut('AllCuts', '')
        return dic

    # ==============================================
    # region GET CONFIG
    def load_config(self):
        self.CutConfig['IndividualChCut'] = ''
        self.CutConfig['JumpExcludeRange'] = loads(self.Config.get('CUT', 'exclude around jump'))
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(self.Config.getfloat('CUT', 'exclude first'))
        self.CutConfig['EventRange'] = self.load_event_range(loads(self.Config.get('CUT', 'event range')))
        self.CutConfig['chi2X'] = self.Config.getint('CUT', 'chi2X')
        self.CutConfig['chi2Y'] = self.Config.getint('CUT', 'chi2Y')
        self.CutConfig['slope'] = self.Config.getint('CUT', 'slope')

    def load_event_range(self, event_range=None):
        """ Gets the event range cut. If the arguments are negative, they are interpreted as time in minutes. Therefore, e.g. load_event_range(-10, 700000) means that only events are considered
        which fulfill: >10 minutes after run start event number < 700000 """
        if event_range is None:
            event_range = [0, 0]
        for i, value in enumerate(event_range):
            if value < 0:
                event_range[i] = self.Analysis.get_event_at_time(seconds=-1 * value * 60)
        if not event_range[1]:
            event_range[1] = self.Analysis.get_event_at_time(-1)
        if not event_range[0]:
            event_range[0] = self.CutConfig['ExcludeFirst']
        return event_range

    def set_event_range(self, event_range):
        self.CutConfig['EventRange'] = self.load_event_range(event_range)

    def load_exclude_first(self, value):
        """ Sets how many events at the very beginning of the run should be excluded. if the argument is negative, it will be interpreted as time in minutes. For a positive argument it is interpreted
         as maximum event number. """
        return value if value >= 0 else self.Analysis.get_event_at_time(-1 * value * 60)

    def set_exclude_first(self, value):
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(value)

    def load_peakpos_high(self, high):
        if high > 0:
            self.EasyCutStrings['peakPos_high'] = 'peakPos<{high}'.format(high=high)
            return high
        else:
            return -1

    def set_peakpos_high(self, value):
        self.CutConfig['peakPos_high'] = self.load_peakpos_high(value)

    def get_fiducial_splits(self):
        return (loads(self.Config.get('SPLIT', 'fiducial')) if self.Config.has_option('SPLIT', 'fiducial') else []) + [int(1e10)]

    # endregion

    def get_event_range(self):
        """
        Returns a the lowest and highest event numbers to consider in the analysis.
        :return: cut eventrange as list, empty if no cut applied
        """
        return self.CutConfig["EventRange"]

    def get_min_event(self):
        """ :return: the smallest event number satisfying the cut conditions. """
        return self.CutConfig["EventRange"][0]

    def get_n_events(self):
        """ :return: number of events in EventRange """
        total_events = self.Analysis.get_event_at_time(-1)
        return total_events if not self.CutConfig["EventRange"] else self.CutConfig["EventRange"][1] - self.CutConfig["EventRange"][0]

    def get_max_event(self):
        """ :return: maximum event number """
        return self.CutConfig["EventRange"][1]

    # ==============================================
    # region GENERATE CUT STRINGS
    def generate_event_range(self):
        cut_string = ''
        if self.CutConfig['EventRange']:
            cut_string = '(event_number<={max}&&event_number>={min})'.format(min=self.CutConfig['EventRange'][0], max=self.CutConfig['EventRange'][1])
        elif self.CutConfig['ExcludeFirst']:
            cut_string = 'event_number>={min}'.format(min=self.CutConfig['ExcludeFirst'])
        return cut_string

    def generate_chi2(self, mode='x', value=None):
        cut_value = self.calc_chi2(mode) if value is None else value
        return 'chi2_{}>=0'.format(mode) + ' && chi2_{mod}<{val}'.format(val=cut_value, mod=mode) if cut_value is not None else ''

    def calc_chi2(self, mode='x'):
        picklepath = self.Analysis.make_pickle_path('Chi2', run=self.RunNumber, suf=mode.title())

        def f():
            t = self.Analysis.info('calculating chi2 cut in {mod} for run {run}...'.format(run=self.Analysis.RunNumber, mod=mode), next_line=False)
            h = TH1F('hc{}'.format(mode), '', 500, 0, 100)
            self.Analysis.Tree.Draw('chi2_{m}>>hc{m}'.format(m=mode), 'n_tracks > 0', 'goff')
            chi2s = zeros(100)
            h.GetQuantiles(100, chi2s, arange(.01, 1.01, .01))
            self.Analysis.add_to_info(t)
            return chi2s

        chi2 = do_pickle(picklepath, f)
        quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
        assert isint(quantile) and 0 < quantile <= 100, 'chi2 quantile has to be and integer between 0 and 100'
        return chi2[quantile] if quantile != 100 else None

    def generate_slope(self, mode='x'):
        cut_variable = '{t}_{m}'.format(t='slope' if self.Analysis.Run.has_branch('slope_x') else 'angle', m=mode)
        angles = self.calc_angle(mode)
        string = '{v}>{min}&&{v}<{max}'.format(v=cut_variable, min=angles[mode][0], max=angles[mode][1])
        return string if self.CutConfig['slope'] > 0 else ''

    def calc_angle(self, mode='x'):
        # take the pickle of the run with a low rate if provided (for ana collection)
        run = self.LowRateRun if self.LowRateRun is not None else self.RunNumber
        picklepath = self.Analysis.make_pickle_path('TrackAngle', mode, run=run)

        def func():
            angle = self.CutConfig['slope']
            t = self.Analysis.info('Generating angle cut in {m} for run {run} ...'.format(run=self.Analysis.RunNumber, m=mode), False)
            set_root_output(False)
            h = self.Analysis.draw_angle_distribution(mode=mode, show=False, print_msg=False)
            fit = fit_fwhm(h)
            mean_ = fit.Parameter(1)
            cut_vals = {mode: [mean_ - angle, mean_ + angle]}
            self.Analysis.add_to_info(t)
            return cut_vals

        return do_pickle(picklepath, func)

    @staticmethod
    def generate_distance(dmin, dmax, thickness=500):
        d_string = '{t}*TMath::Sqrt(TMath::Power(TMath::Sin(TMath::DegToRad()*slope_x), 2) + TMath::Power(TMath::Sin(TMath::DegToRad()*slope_y), 2) + 1)'.format(t=thickness)
        return TCut('distance', '{d}>{min}&&{d}<={max}'.format(d=d_string, min=dmin, max=dmax))

    def generate_cut_string(self):
        """ Creates the cut string. """
        gROOT.SetBatch(1)

        # --TRACKS --
        self.CutStrings['chi2X'] += self.generate_chi2('x')
        self.CutStrings['chi2Y'] += self.generate_chi2('y')
        self.CutStrings['slope_x'] += self.generate_slope('x')
        self.CutStrings['slope_y'] += self.generate_slope('y')
        self.CutStrings['tracks'] += 'n_tracks==1'

        # -- EVENT RANGE CUT --
        self.CutStrings['event_range'] += self.generate_event_range()
        if self.CutConfig['EventRange']:
            self.EasyCutStrings['EventRange'] = 'Evts.{min}k-{max}k'.format(min=int(self.CutConfig['EventRange'][0]) / 1000, max=int(self.CutConfig['EventRange'][1]) / 1000)
            self.EasyCutStrings['ExcludeFirst'] = 'Evts.{min}k+'.format(min=int(self.CutConfig['ExcludeFirst']) / 1000) if self.CutConfig['ExcludeFirst'] > 0 else ''

        # -- BEAM INTERRUPTION CUT --
        self.CutStrings['beam_interruptions'] += self.generate_beam_interruptions()
        self.JumpCut += self.generate_jump_cut()

        gROOT.SetBatch(0)

    def generate_beam_interruptions(self):
        """ This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied. """
        interruptions = self.get_beam_interruptions()
        cut_string = TCut('')
        for interr in interruptions:
            cut_string += TCut('event_number<{low}||event_number>{high}'.format(low=interr[0], high=interr[1]))
        self.EasyCutStrings['noBeamInter'] = 'BeamOn'
        return TCut(cut_string)

    def generate_jump_cut(self):
        cut_string = ''
        start_event = self.CutConfig['EventRange'][0]
        for tup in self.Jumps:
            if tup[1] > start_event:
                low = start_event if tup[0] < start_event else tup[0]
                cut_string += '&&' if cut_string else ''
                cut_string += '!(event_number<={up}&&event_number>={low})'.format(up=tup[1], low=low)
        return TCut(cut_string)

    def generate_flux_cut(self):
        return self.generate_special_cut(included=['beam_interruptions', 'event_range'], name='flux', prnt=False)
    # endregion

    # ==============================================
    # region BEAM INTERRUPTS
    def find_beam_interruptions(self):
        dut_type = self.Analysis.Run.Config.get('BASIC', 'type')
        return self.find_pad_beam_interruptions() if dut_type == 'pad' else self.find_pixel_beam_interruptions()

    def find_pixel_beam_interruptions(self, bin_width=10, threshold=.4):
        """ Finding beam interruptions by incestigation the event rate. """
        t_start = self.Analysis.info('Searching for beam interruptions of run {r} ...'.format(r=self.RunNumber), next_line=False)
        bin_values, time_bins = histogram(self.Analysis.Run.Time / 1000, bins=self.Bins.get_raw_time_binning(bin_width)[1])
        m = mean(bin_values[bin_values.argsort()][-20:-10])  # take the mean of the 20th to the 10th highest bin to get an estimate of the plateau
        deviating_bins = where(abs(1 - bin_values / m) > threshold)[0]
        times = time_bins[deviating_bins] + bin_width / 2 - self.Analysis.Run.Time[0] / 1000  # shift to the center of the bin
        not_connected = where(concatenate([[False], deviating_bins[:-1] != deviating_bins[1:] - 1]))[0]  # find the bins that are not consecutive
        times = split(times, not_connected)
        jumps = [[self.Analysis.get_event_at_time(v) for v in [t[0], t[0] if t.size == 1 else t[-1]]] for t in times]
        interruptions = self.__create_jump_ranges(jumps)
        self.Analysis.add_to_info(t_start)
        return jumps, interruptions

    def find_pad_beam_interruptions(self, bin_width=100, max_thresh=.4):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        t = self.Analysis.info('Searching for beam interruptions of run {r} ...'.format(r=self.RunNumber), next_line=False)
        n = self.Analysis.Tree.Draw('Entry$:pulser', '', 'goff')
        x, y = get_root_vecs(self.Analysis.Tree, n, 2, dtype=int)
        rates, x_bins, y_bins = histogram2d(x, y, bins=[arange(0, n, bin_width, dtype=int), 2])
        rates = rates[:, 1] / bin_width
        thresh = min(max_thresh, mean(rates) + .1)
        events = x_bins[:-1][rates > thresh] + bin_width / 2
        not_connected = where(concatenate([[False], events[:-1] != events[1:] - bin_width]))[0]  # find the events where the previous event is not related to the event (more than a bin width away)
        events = split(events, not_connected)  # events grouped into connecting events
        jumps = [(ev[0], ev[0]) if ev.size == 1 else (ev[0], ev[-1]) for ev in events] if events[0].size else []
        interruptions = self.__create_jump_ranges(jumps)
        self.Analysis.add_to_info(t)
        return jumps, interruptions

    def __create_jump_ranges(self, jumps):
        interruptions = []
        i = 0
        for tup in jumps:
            t_start = self.Analysis.Run.get_time_at_event(tup[0]) - self.Analysis.Run.StartTime - self.CutConfig['JumpExcludeRange'][0]
            t_stop = self.Analysis.Run.get_time_at_event(tup[1]) - self.Analysis.Run.StartTime + self.CutConfig['JumpExcludeRange'][1]
            # if interruptions overlay just set the last stop to the current stop
            if i and t_start <= (interruptions[i - 1][1]) + 10:
                interruptions[i - 1][1] = t_stop
                continue
            interruptions.append([t_start, t_stop])
            i += 1
        return [[self.Analysis.Run.get_event_at_time(t) for t in tup] for tup in interruptions]

    def get_beam_interruptions(self):
        """ The data is stored as a list of jumps, dumped into a pickle file. If no pickle file exists, it will perform a beam interruption analysis in order to identify the beam interruptions.
        :returns: list of interruptions including safety margin from the AnalysisConfig. """
        pickle_path = self.Analysis.make_pickle_path('BeamInterruptions', run=self.RunNumber, suf='_'.join(str(i) for i in self.CutConfig['JumpExcludeRange']))
        interruptions = do_pickle(pickle_path, self.find_beam_interruptions)
        self.Jumps = interruptions[0]
        self.Interruptions = interruptions[1]
        return interruptions[1]
    # endregion

    def get_easy_cutstring(self):
        """
        Returns a short, more user-friendly cut string, which can be used to display the cut configuration as terminal prompt or inside a canvas.
        :return:
        """
        string_ = ""
        for type_ in self.EasyCutStrings.keys():
            if self.EasyCutStrings[type_] != "":
                string_ += self.EasyCutStrings[type_] + ", "
        if string_ != "":
            string_ = string_[:-2]
        return string_

    def reset_cut(self, name):
        if name in self.CutStrings:
            self.CutStrings[name].SetTitle('')
        else:
            print 'There is no cut with the name "{name}"!'.format(name=name)
        self.update_all_cut()

    def update_cut(self, name, value=None):
        if name in self.CutStrings:
            self.CutStrings[name].SetTitle('')
            self.CutStrings[name] += value
            self.update_all_cut()
        else:
            print 'There is no cut with the name "{name}"!'.format(name=name)

    def set_chi2(self, value):
        self.CutConfig['chi2X'] = value
        self.CutConfig['chi2Y'] = value
        self.update_cut('chi2X', self.generate_chi2('x'))
        self.update_cut('chi2Y', self.generate_chi2('y'))

    def update_all_cut(self):
        self.AllCut = self.generate_all_cut()
        self.Analysis.AllCuts = self.AllCut

    def show_cuts(self, easy=True):
        cuts = self.EasyCutStrings if easy else self.CutStrings
        max_len = max(len(key) for key, value in cuts.iteritems() if str(value))
        for key, value in cuts.iteritems():
            if not key == 'AllCuts' and str(value):
                print '{key}:'.format(key=key.rjust(max_len)), value
        return

    @staticmethod
    def get_track_var(num, mode, scale=1):
        return 'dia_track_{m}_local[{n}]*{s}'.format(m=mode, n=num, s=scale)

    def get_track_vars(self, num, scale=1):
        return (self.get_track_var(num, v, scale) for v in ['y', 'x'])

    def generate_consecutive_cuts(self):
        cuts = OrderedDict([('raw', TCut('0', ''))])
        for i, (key, value) in enumerate([(key, value) for key, value in self.CutStrings.iteritems() if str(value) and key != 'AllCuts' and not key.startswith('old')], 1):
            new_cut = cuts.values()[i - 1] + value
            key = 'beam_stops' if 'beam' in key else key
            cuts[key] = TCut('{n}'.format(n=i), str(new_cut))
        return cuts

    def draw_contributions(self, flat=False, short=False, show=True):
        set_root_output(show)
        contr = OrderedDict()
        n_events = self.Analysis.Run.NEntries
        cut_events = 0
        for i, (key, cut) in enumerate(self.generate_consecutive_cuts().iteritems()):
            if key == 'raw':
                continue
            events = n_events - int(self.Analysis.Tree.Draw('1', cut, 'goff'))
            print key.rjust(18), '{0:5d} {1:04.1f}%'.format(events - cut_events, (1. - float(events) / n_events) * 100.)
            contr[key.title().replace('_', ' ')] = (events - cut_events, self.Analysis.get_color())
            cut_events = events
        contr['Good Events'] = (n_events - cut_events, self.Analysis.get_color())
        print contr
        sorted_contr = OrderedDict(sorted(OrderedDict(item for item in contr.iteritems() if item[1][0] >= (.03 * n_events if short else 0)).iteritems(), key=lambda x: x[1]))  # sort by size
        sorted_contr.update({'Other': (n_events - sum(v[0] for v in sorted_contr.values()), self.Analysis.get_color())} if short else {})
        sorted_contr = OrderedDict(sorted_contr.popitem(not i % 2) for i in xrange(len(sorted_contr)))  # sort by largest->smallest->next largest...
        print sorted_contr
        pie = TPie('pie', 'Cut Contributions', len(sorted_contr), array([v[0] for v in sorted_contr.values()], 'f'), array([v[1] for v in sorted_contr.values()], 'i'))
        for i, label in enumerate(sorted_contr.iterkeys()):
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label)
        format_pie(pie, h=.04, r=.2, text_size=.025, angle3d=70, label_format='%txt (%perc)', angle_off=250)
        self.Analysis.save_histo(pie, draw_opt='{0}rsc'.format('3d' if not flat else ''), show=show)
        self.Analysis.reset_colors()
        return sorted_contr

    def draw_fid_cut(self, scale=1):
        cut = get_object('fid{}'.format(self.RunNumber))
        if cut:
            cut = deepcopy(cut)
            cut.SetName('fid{}'.format(scale))
            for i in xrange(cut.GetN()):
                cut.SetPoint(i, scale * cut.GetX()[i], scale * cut.GetY()[i])
            cut.Draw()
            self.Analysis.Objects.append(cut)
