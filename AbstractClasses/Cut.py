import os
import pickle
import json
import ConfigParser
from numpy import mean, array, zeros, arange, delete
from AbstractClasses.Elementary import Elementary
# from newAnalysis import Analysis
# from signal_analysis import SignalAnalysis
from ROOT import TCut, gROOT, TH1F
from collections import OrderedDict


class Cut(Elementary):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, parent_analysis, verbose=True):

        self.analysis = parent_analysis
        self._checklist = {"RemoveBeamInterruptions": False,
                           "GenerateCutString": False}
        # saving stuff
        self.histos = {}

        # config
        self.parser = self.load_parser()
        self.beaminterruptions_folder = self.parser.get('CUT', 'beaminterruptions_folder')
        self.exclude_before_jump = self.parser.getint('CUT', 'excludeBeforeJump')
        self.exclude_after_jump = self.parser.getint('CUT', 'excludeAfterJump')
        self.CutConfig = {}

        # define cut strings
        self.EasyCutStrings = self.init_easy_cutstrings()
        self.CutStrings = self.define_cutstrings()

        self.region_cut = TCut('region_cut', '')

        # beam interrupts
        self.jumps = None
        self.jump_ranges = None

        Elementary.__init__(self, verbose=verbose)

        # generate cut strings
        self.generate_cut_string()
        self.all_cut = self.generate_all_cut()

    def generate_all_cut(self):
        cut = TCut('all_cuts', '')
        for key, value in self.CutStrings.iteritems():
            if not key.startswith('old'):
                cut += value
        return cut

    def get_included_events(self, maxevent=None):
        """
        :param maxevent:
        :return: list of included event numbers not excluded by: excludeFirst, EventRange or BeamInterruptions
        """
        minevent = self.get_min_event()
        maxevent = self.get_max_event() if maxevent is None else maxevent

        excluded = [i for i in arange(0, minevent)]  # first events
        for start, stop in zip(self.jump_ranges['start'], self.jump_ranges['stop']):
            excluded += [i for i in xrange(start, stop + 1)]  # events around jumps
        excluded.sort()
        all_events = arange(0, maxevent)
        included = delete(all_events, excluded)
        return included
    
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
        dic['spread_low'] = ''
        dic['absMedian_high'] = ''
        dic['pedestalsigma'] = ''
        return dic
    
    @staticmethod
    def define_cutstrings():
        dic = OrderedDict()
        dic['raw'] = TCut('raw', '')
        dic['pulser'] = TCut('pulser', '')
        dic['event_range'] = TCut('event_range', '')
        dic['beam_interruptions'] = TCut('beam_interruptions', '')
        dic['ped_sigma'] = TCut('ped_sigma', '')
        dic['spread_low'] = TCut('spread_low', '')
        dic['median'] = TCut('median', '')
        dic['tracks'] = TCut('tracks', '')
        dic['chi2'] = TCut('chi2', '')
        dic['track_angle'] = TCut('track_angle', '')
        dic['saturated'] = TCut('saturated', '')
        dic['old_bucket'] = TCut('old_bucket', '')
        dic['bucket'] = TCut('bucket', '')
        dic['all_cuts'] = TCut('all_cuts', '')
        return dic

    # ==============================================
    # region GET CONFIG
    def load_parser(self):
        parser = ConfigParser.ConfigParser()
        parser.read('Configuration/AnalysisConfig_' + self.analysis.TESTCAMPAIGN + '.cfg')
        return parser

    def load_config(self):
        self.CutConfig['IndividualChCut'] = ''
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(self.parser.getint('CUT', 'excludefirst'))
        self.CutConfig['EventRange'] = self.load_event_range(json.loads(self.parser.get('CUT', 'EventRange')))
        self.CutConfig['spread_low'] = self.load_spread_low(self.parser.getint('CUT', 'spread_low'))
        self.CutConfig['absMedian_high'] = self.load_abs_median_high(self.parser.getint('CUT', 'absMedian_high'))
        self.CutConfig['pedestalsigma'] = self.load_pedestal_sigma(self.parser.getint('CUT', 'pedestalsigma'))
        self.CutConfig['chi2'] = self.parser.getint('CUT', 'chi2')
        self.CutConfig['track_angle'] = self.parser.getint('CUT', 'track_angle')

    def load_event_range(self, event_range=None):
        """
        Gets the event range cut. If the arguments are negative, they are interpreted as time in minutes. Therefore, e.g.
        load_event_range(-10, 700000) means that only events are considered, which fulfill: >10 minutes after run start event number < 700000
        :param event_range:
        :return: event range
        """
        if event_range is None:
            event_range = [0, 0]
        for i, value in enumerate(event_range):
            if value < 0:
                event_range[i] = self.analysis.get_event_at_time(time_sec=-1 * value * 60)
        if not event_range[1]:
            event_range[1] = self.analysis.get_event_at_time(-1)
        if not event_range[0]:
            event_range[0] = self.CutConfig['ExcludeFirst']
        return event_range

    def set_event_range(self, event_range):
        self.CutConfig['EventRange'] = self.load_event_range(event_range)

    def load_exclude_first(self, value):
        """
        Sets how many events at the very beginning of the run should be excluded. if the argument is negative, it will be interpreted as time in minutes. For a positive argument it is interpreted as
        maximum event number.
        :param value: events or time in minutes
        :return:
        """
        if value > 0:
            self.EasyCutStrings['ExcludeFirst'] = str(int(value) / 1000) + 'k+'
            return value
        elif value == 0:
            self.EasyCutStrings['ExcludeFirst'] = ''
            return 0
        else:
            self.EasyCutStrings['ExcludeFirst'] = str(-1 * value) + 'min+'
            seconds = -1 * value * 60
            event = self.analysis.get_event_at_time(seconds)
            return event

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

    def load_spread_low(self, value):
        if value > 0:
            self.EasyCutStrings['spread_low'] = 'spread>{low}'.format(low=value)
            return value
        else:
            return -1

    def set_spread_low(self, low):
        self.CutConfig['spread_low'] = self.load_spread_low(low)

    def load_abs_median_high(self, value):
        if value > 0:
            self.EasyCutStrings['absMedian_high'] = '|median|<{high}'.format(high=value)
            return value
        else:
            return -1

    def set_abs_median_high(self, high):
        self.CutConfig['absMedian_high'] = self.load_abs_median_high(high)

    def load_pedestal_sigma(self, value):
        if value > 0:
            self.EasyCutStrings['pedestalsigma'] = 'PedSigma' + str(value)
            return value
        else:
            self.EasyCutStrings['pedestalsigma'] = ''
            return -1

    def set_pedestal_sigma(self, sigma=-1):
        self.CutConfig['pedestalsigma'] = self.load_pedestal_sigma(sigma)

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
        total_events = self.analysis.get_event_at_time(-1)
        return total_events if not self.CutConfig["EventRange"] else self.CutConfig["EventRange"][1] - self.CutConfig["EventRange"][0]

    def get_max_event(self):
        """ :return: maximum event number """
        return self.CutConfig["EventRange"][1]

    # ==============================================
    # region GENERATE CUT STRINGS
    def generate_event_range(self):
        if self.CutConfig['EventRange']:
            self.CutStrings['event_range'] += '(event_number<={max}&&event_number>={min})'.format(min=self.CutConfig['EventRange'][0], max=self.CutConfig['EventRange'][1])
        elif self.CutConfig['ExcludeFirst']:
            self.CutStrings['event_range'] += 'event_number>={min}'.format(min=self.CutConfig['ExcludeFirst'])

    def generate_chi2(self):
        picklepath = 'Configuration/Individual_Configs/Chi2/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number)

        def func():
            print 'generating chi2 cut for run {run}...'.format(run=self.analysis.run_number)
            gROOT.SetBatch(1)
            h = TH1F('h', '', 200, 0, 100)
            nq = 100
            chi2s = zeros(nq)
            xq = array([(i + 1) / float(nq) for i in range(nq)])
            self.analysis.tree.Draw('chi2_tracks>>h', '', 'goff')
            h.GetQuantiles(nq, chi2s, xq)
            gROOT.SetBatch(0)
            return chi2s

        chi2 = self.do_pickle(picklepath, func)
        assert type(self.CutConfig['chi2']) is int and 0 < self.CutConfig['chi2'] <= 100, 'chi2 quantile has to be and integer between 0 and 100'
        string = 'chi2_tracks<{val}&&chi2_tracks>=0'.format(val=chi2[self.CutConfig['chi2']])
        return TCut(string) if self.CutConfig['chi2'] > 0 else ''

    def generate_slope(self):
        picklepath = 'Configuration/Individual_Configs/Slope/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.lowest_rate_run)
        angle = self.CutConfig['track_angle']

        def func():
            print 'generating slope cut for run {run}...'.format(run=self.analysis.run_number)
            # fit the slope to get the mean
            gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
            gROOT.SetBatch(1)
            h_x = TH1F('hx', '', 70, -4, 4)
            h_y = TH1F('hy', '', 70, -4, 4)
            self.analysis.tree.Draw('slope_x>>hx', '', 'goff')
            self.analysis.tree.Draw('slope_y>>hy', '', 'goff')
            fit_result = h_x.Fit('gaus', 'qs')
            slopes = {'x': [], 'y': []}
            x_mean = fit_result.Parameters()[1]
            slopes['x'] = [x_mean - angle, x_mean + angle]
            fit_result = h_y.Fit('gaus', 'qs')
            y_mean = fit_result.Parameters()[1]
            slopes['y'] = [y_mean - angle, y_mean + angle]
            c = gROOT.FindObject('c1')
            c.Close()
            gROOT.SetBatch(0)
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            return slopes

        slope = self.do_pickle(picklepath, func)
        # create the cut string
        string = 'slope_x>{minx}&&slope_x<{maxx}&&slope_y>{miny}&&slope_y<{maxy}'.format(minx=slope['x'][0], maxx=slope['x'][1], miny=slope['y'][0], maxy=slope['y'][1])
        return TCut(string) if angle > 0 else ''
    
    def generate_cut_string(self):
        """ Creates the cut string. """
        gROOT.SetBatch(1)

        # --TRACKS --
        self.CutStrings['chi2'] = self.generate_chi2()
        self.CutStrings['track_angle'] = self.generate_slope()
        self.CutStrings['tracks'] += 'n_tracks'

        # -- EVENT RANGE CUT --
        self.generate_event_range()
        if self.CutConfig['EventRange']:
            self.EasyCutStrings['EventRange'] = 'Evts.{min}k-{max}k'.format(min=int(self.CutConfig['EventRange'][0]) / 1000, max=int(self.CutConfig['EventRange'][1]) / 1000)
            self.EasyCutStrings['ExcludeFirst'] = 'Evts.{min}k+'.format(min=int(self.CutConfig['ExcludeFirst']) / 1000) if self.CutConfig['ExcludeFirst'] > 0 else ''

        # -- PULSER CUT --
        self.CutStrings['pulser'] += '!pulser'

        # -- BEAM INTERRUPTION CUT --
        self.__generate_beam_interruptions()
        self.EasyCutStrings['noBeamInter'] = 'BeamOn'

        self._checklist['GenerateCutString'] = True
        gROOT.SetBatch(0)

    def __generate_beam_interruptions(self, ):
        """
        This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied.
        """
        self.get_beam_interruptions()

        njumps = len(self.jump_ranges["start"])
        cut_string = ''
        for i in xrange(njumps):
            string = "!(event_number<={upper}&&event_number>={lower})".format(upper=self.jump_ranges["stop"][i], lower=self.jump_ranges["start"][i])
            # new separate strings
            if cut_string != '':
                cut_string += '&&'
            cut_string += string
        self.CutStrings['beam_interruptions'] += cut_string
        self._checklist["RemoveBeamInterruptions"] = True
    # endregion

    # ==============================================
    # region BEAM INTERRUPTS
    def find_beam_interruptions(self):
        """
        Looking for the beam interruptions
        :return: interrupt list
        """
        print 'Searching for beam interruptions...'
        binning = 100
        nbins = int(self.analysis.run.tree.GetEntries()) / binning
        rate = []
        for i in xrange(nbins):
            pulserevents = self.analysis.run.tree.Draw("1", "pulser", "", binning, i * binning)
            rate.append(1. * pulserevents / binning)
        mean_rate = mean(rate)
        interrupts = []
        last_rate = 0
        tup = [0, 0]
        fac = 2
        for i, value in enumerate(rate):
            if value > fac * mean_rate > last_rate:
                tup[0] = i * binning
            elif value < fac * mean_rate < last_rate:
                tup[1] = i * binning
                interrupts.append(tup)
                tup = [0, 0]
            last_rate = value
        return interrupts

    def __save_beaminterrupts(self):
        # check if directories exist
        if not os.path.exists(self.beaminterruptions_folder):
            os.mkdir(self.beaminterruptions_folder)
        if not os.path.exists(self.beaminterruptions_folder + "/data"):
            os.mkdir(self.beaminterruptions_folder + "/data")

        # save jump list to file
        jumpfile = open(self.beaminterruptions_folder + "/data/{testcampaign}Run_{run}.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number), "wb")
        pickle.dump(self.jumps, jumpfile)
        jumpfile.close()

    def __create_jump_ranges(self):
        if self.jump_ranges is None and len(self.jumps) > 0:
            print 'generating jump ranges...'
            start = []
            stop = []
            time_offset = self.analysis.run.get_time_at_event(0)
            t_max = (self.analysis.run.get_time_at_event(-1) - time_offset) / 1000.
            last_stop = 0
            for tup in self.jumps:
                t_start = (self.analysis.run.get_time_at_event(tup[0]) - time_offset) / 1000.
                t_stop = (self.analysis.run.get_time_at_event(tup[1]) - time_offset) / 1000.
                # add offsets from config file
                t_start -= -1 * self.exclude_before_jump if t_start >= -1 * self.exclude_before_jump else 0
                t_stop = t_stop + -1 * self.exclude_after_jump if t_stop + -1 * self.exclude_after_jump <= t_max else t_max
                if t_start < last_stop:
                    stop[-1] = self.analysis.get_event_at_time(t_stop)
                    last_stop = t_stop
                    continue
                start.append(self.analysis.get_event_at_time(t_start))
                stop.append(self.analysis.get_event_at_time(t_stop))
                last_stop = t_stop

            self.jump_ranges = {"start": start,
                                "stop": stop}

        return [self.exclude_before_jump, self.exclude_after_jump, self.jump_ranges]

    def get_beam_interruptions(self):
        """
        If beam interruption data exist in beaminterruptions/data/, it will load it in order to account for beam interruptions. The data is stored as a list of jumps, dumped into a pickle file.
        If no pickle file exists, it will perform a beam interruption analysis in order to identify the beam interruptions. The found interruptions are stored in a list at .jumps and dumped into
        a pickle file.
        :return: list of events where beam interruptions occures
        """
        if self.jump_ranges is None:
            jumps_pickle = self.beaminterruptions_folder + "/data/{testcampaign}Run_{run}.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
            range_pickle = self.beaminterruptions_folder + "/data/{testcampaign}_{run}_Jump_Ranges.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
            self.jumps = self.do_pickle(jumps_pickle, self.find_beam_interruptions)
            ranges = self.do_pickle(range_pickle, self.__create_jump_ranges)
            # redo range pickle if config parameters have changed
            if ranges[0] != self.exclude_before_jump or ranges[1] != self.exclude_after_jump:
                os.remove(range_pickle)
                ranges = self.do_pickle(range_pickle, self.__create_jump_ranges)
            self.jump_ranges = ranges[2]
        return self.jumps
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

    def show_cuts(self, easy=True):
        cuts = self.EasyCutStrings if easy else self.CutStrings
        max_len = max(len(key) for key, value in cuts.iteritems() if str(value))
        for key, value in cuts.iteritems():
            if not key == 'all_cuts' and str(value):
                print '{key}:'.format(key=key.rjust(max_len)), value
        return
