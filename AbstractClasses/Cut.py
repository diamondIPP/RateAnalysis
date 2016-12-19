import json
from numpy import array, zeros
from Elementary import Elementary
from ROOT import TCut, gROOT, TH1F
from collections import OrderedDict
from Utils import *


class Cut(Elementary):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, parent_analysis, skip=False):

        if not skip:
            self.analysis = parent_analysis
            Elementary.__init__(self, verbose=self.analysis.verbose)
            self.RunNumber = self.analysis.run_number

            # saving stuff
            self.RootObjects = []

            # config
            self.DUTType = self.load_dut_type()
            self.BeaminterruptionsDir = self.ana_config_parser.get('CUT', 'beaminterruptions_folder')
            self.CutConfig = {}

            # define cut strings
            self.EasyCutStrings = self.init_easy_cutstrings()
            self.CutStrings = self.define_cutstrings()

            self.region_cut = TCut('region_cut', '')
            self.JumpCut = TCut('JumpCut', '')

            # beam interrupts
            self.Jumps = None
            self.Interruptions = None

            self.load_config()
            # generate cut strings
            self.generate_cut_string()
            self.all_cut = self.generate_all_cut()

    def load_run_config(self):
        return self.load_run_configs(self.analysis.run_number)

    def generate_special_cut(self, excluded=None, included=None, name='special_cut'):
        cut = TCut(name, '')
        n_cuts = 0
        for key, value in self.CutStrings.iteritems():
            if excluded and key in excluded:
                continue
            if included and key not in included:
                continue
            if key.startswith('old') or key.startswith('all_cut'):
                continue
            if value.GetTitle() == '':
                continue
            cut += value
            n_cuts += 1
        self.log_info('generated {name} cut with {num} cuts'.format(name=name, num=n_cuts))
        return cut

    def generate_all_cut(self):
        cut = TCut('all_cuts', '')
        for key, value in self.CutStrings.iteritems():
            if not key.startswith('old') and not key.startswith('all_cut'):
                cut += value
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
        dic['spread_low'] = ''
        dic['absMedian_high'] = ''
        dic['pedestalsigma'] = ''
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
        dic['tracks'] = TCut('tracks', '')
        dic['hit'] = TCut('hit', '')
        dic['masks'] = TCut('masks', '')
        dic['fiducial'] = TCut('fiducial', '')
        dic['chi2X'] = TCut('chi2X', '')
        dic['chi2Y'] = TCut('chi2Y', '')
        dic['track_angle_x'] = TCut('track_angle_x', '')
        dic['track_angle_y'] = TCut('track_angle_y', '')
        dic['rhit'] = TCut('rhit', '')
        dic['all_cuts'] = TCut('all_cuts', '')
        return dic

    # ==============================================
    # region GET CONFIG
    
    def load_dut_type(self):
        dut_type = self.run_config_parser.get("BASIC", "type")
        assert dut_type.lower() in ["pixel", "pad"], "The DUT type {0} should be 'pixel' or 'pad'".format(dut_type)
        return dut_type

    def load_config(self):
        self.CutConfig['IndividualChCut'] = ''
        self.CutConfig['JumpExcludeRange'] = {'before': self.ana_config_parser.getint('CUT', 'excludeBeforeJump'), 'after': self.ana_config_parser.getint('CUT', 'excludeAfterJump')}
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(self.ana_config_parser.getint('CUT', 'excludefirst'))
        self.CutConfig['EventRange'] = self.load_event_range(json.loads(self.ana_config_parser.get('CUT', 'EventRange')))
        self.CutConfig['chi2X'] = self.ana_config_parser.getint('CUT', 'chi2X')
        self.CutConfig['chi2Y'] = self.ana_config_parser.getint('CUT', 'chi2Y')
        self.CutConfig['track_angle'] = self.ana_config_parser.getint('CUT', 'track_angle')

    def load_event_range(self, event_range=None):
        """ Gets the event range cut. If the arguments are negative, they are interpreted as time in minutes. Therefore, e.g. load_event_range(-10, 700000) means that only events are considered
        which fulfill: >10 minutes after run start event number < 700000 """
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
        """ Sets how many events at the very beginning of the run should be excluded. if the argument is negative, it will be interpreted as time in minutes. For a positive argument it is interpreted
         as maximum event number. """
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

    def generate_chi2(self, mode='x'):
        picklepath = 'Configuration/Individual_Configs/Chi2/{tc}_{run}_{mod}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number, mod=mode.title())

        def func():
            print 'generating chi2 cut in {mod} for run {run}...'.format(run=self.analysis.run_number, mod=mode)
            gROOT.SetBatch(1)
            h = TH1F('h', '', 200, 0, 100)
            nq = 100
            chi2s = zeros(nq)
            xq = array([(i + 1) / float(nq) for i in range(nq)])
            self.analysis.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
            h.GetQuantiles(nq, chi2s, xq)
            gROOT.SetBatch(0)
            return chi2s

        chi2 = self.do_pickle(picklepath, func)
        quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
        assert type(quantile) is int and 0 < quantile <= 100, 'chi2 quantile has to be and integer between 0 and 100'
        string = 'chi2_{mod}<{val}&&chi2_{mod}>=0'.format(val=chi2[quantile], mod=mode)
        return string if quantile > 0 else ''

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
            self.analysis.tree.Draw('slope_x>>hx', 'slope_x > -100', 'goff')
            self.analysis.tree.Draw('slope_y>>hy', 'slope_y > -100', 'goff')
            if h_x.GetEntries() > 500 and h_y.GetEntries() > 500:
                fit_result = h_x.Fit('gaus', 'qs')
            else:
                log_warning('Empty slope histogram! Using default values!')
                return {'x': [-4., 4.], 'y': [-4., 4.]}
            slopes = {'x': [], 'y': []}
            x_mean = fit_result.Parameters()[1]
            slopes['x'] = [x_mean - angle, x_mean + angle]
            fit_result = h_y.Fit('gaus', 'qs')
            y_mean = fit_result.Parameters()[1]
            slopes['y'] = [y_mean - angle, y_mean + angle]
            gROOT.SetBatch(0)
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            return slopes

        slope = self.do_pickle(picklepath, func)
        # create the cut string
        string = 'slope_x>{minx}&&slope_x<{maxx}&&slope_y>{miny}&&slope_y<{maxy}'.format(minx=slope['x'][0], maxx=slope['x'][1], miny=slope['y'][0], maxy=slope['y'][1])
        return string if angle > 0 else ''

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
        self.CutStrings['track_angle'] += self.generate_slope()
        self.CutStrings['tracks'] += 'n_tracks'

        # -- EVENT RANGE CUT --
        self.generate_event_range()
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
            cut_string += TCut('event_number<{low}||event_number>{high}'.format(low=interr['i'], high=interr['f']))
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
    # endregion

    # ==============================================
    # region BEAM INTERRUPTS
    def find_beam_interruptions(self):
        dic = self.CutConfig['JumpExcludeRange']
        interruptions = self.find_pad_beam_interruptions() if self.DUTType == 'pad' else self.find_pixel_beam_interruptions()
        dic['interruptions'] = interruptions[0]
        dic['jumps'] = interruptions[1]
        return dic

    def find_pixel_beam_interruptions(self, show=False):
        """ Locates the beam interruptions and cuts some seconds before and some seconds after which are specified in the config file. overlapping segments are fused to one to reduce the size
            of the string. An average of the number of event per time bin is calculated and every bin below 90% or above 120% is excluded """
        # time is in minutes. good results found with bin size of 10 seconds
        bin_size = 10  # seconds
        bins = int(self.analysis.run.totalMinutes * 60 / bin_size)
        self.set_root_output(0)
        h1 = TH1F('h_beam_time_', 'h_beam_time_', bins, 0, self.analysis.run.totalMinutes)
        self.analysis.tree.Draw('(time - {off}) / 60000.>>h_beam_time_'.format(off=self.analysis.run.startTime), '', 'goff')
        self.draw_histo(h1)
        mean = h1.Integral() / float(h1.GetNbinsX())
        interruptions = []
        jumps = []
        n_interr = 0
        ex_range = {key: value / 60. for key, value in self.CutConfig['JumpExcludeRange'].iteritems()}
        for t in xrange(1, h1.GetNbinsX() + 1):
            if h1.GetBinContent(t) < mean * .6 or h1.GetBinContent(t) > mean * 1.3:
                low_edge = h1.GetBinLowEdge(t)
                bin_width = h1.GetBinWidth(t)
                if not n_interr or (not low_edge - ex_range['before'] < interruptions[n_interr - 1]['f'] and n_interr) :
                    interruptions.append({'i': low_edge - ex_range['before'], 'f': low_edge + bin_width + ex_range['after']})
                    jumps.append([low_edge, low_edge + bin_width])
                    n_interr += 1
                else:
                    interruptions[n_interr - 1]['f'] = low_edge + bin_width + ex_range['after']
                    jumps[n_interr - 1][1] = low_edge + bin_width
        interruptions = [{key: self.analysis.get_event_at_time(value * 60) for key, value in dic.iteritems()} for dic in interruptions]
        jumps = [[self.analysis.get_event_at_time(t * 60) for t in jump] for jump in jumps]
        return interruptions, jumps

    def find_pad_beam_interruptions(self):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        print 'Searching for beam interruptions...'
        binning = 200
        nbins = int(self.analysis.run.tree.GetEntries()) / binning
        rate = []
        for i in xrange(nbins):
            pulserevents = self.analysis.run.tree.Draw('1', 'pulser', 'goff', binning, i * binning)
            rate.append(100 * pulserevents / binning)
        jumps = []
        last_rate = 0
        tup = [0, 0]
        cut = 40  # if rate goes higher than n %
        for i, value in enumerate(rate):
            if value >= cut > last_rate:
                tup[0] = i * binning
            elif value < cut <= last_rate:
                tup[1] = i * binning
                jumps.append(tup)
                tup = [0, 0]
            last_rate = value
        interruptions = self.__create_jump_ranges(jumps)
        return interruptions, jumps

    def __create_jump_ranges(self, jumps):
        ex_range = self.CutConfig['JumpExcludeRange']
        interruptions = []
        print 'generating jump ranges...'
        time_offset = self.analysis.run.startTime
        t_max = (self.analysis.run.get_time_at_event(-1) - time_offset) / 1000.
        last_stop = 0
        n_interr = 0
        for tup in jumps:
            t_start = (self.analysis.run.get_time_at_event(tup[0]) - time_offset) / 1000.
            t_stop = (self.analysis.run.get_time_at_event(tup[1]) - time_offset) / 1000.
            # add additional time around jumps to be safe
            t_start -= ex_range['before'] if t_start >= ex_range['before'] else 0
            t_stop = t_stop + ex_range['after'] if t_stop + ex_range['after'] <= t_max else t_max
            if t_start < last_stop:
                interruptions[n_interr - 1]['f'] = self.analysis.get_event_at_time(t_stop)
                last_stop = t_stop
                continue
            interruptions.append({'i': self.analysis.get_event_at_time(t_start), 'f': self.analysis.get_event_at_time(t_stop)})
            last_stop = t_stop
            n_interr += 1
        return interruptions

    def get_beam_interruptions(self):
        """
        If beam interruption data exist in beaminterruptions/data/, it will load it in order to account for beam interruptions. The data is stored as a list of jumps, dumped into a pickle file.
        If no pickle file exists, it will perform a beam interruption analysis in order to identify the beam interruptions. The found interruptions are stored in a list at .jumps and dumped into
        a pickle file.
        """
        # check if directories exist
        ensure_dir(self.BeaminterruptionsDir)
        ensure_dir(joinpath(self.BeaminterruptionsDir, 'data'))
        pickle_path = self.make_pickle_path('BeamInterruptions', run=self.RunNumber)
        interruptions = self.do_pickle(pickle_path, self.find_beam_interruptions)

        # redo range pickle if config parameters have changed
        ex_range = self.CutConfig['JumpExcludeRange']
        if interruptions['before'] != ex_range['before'] or interruptions['after'] != ex_range['after']:
            os.remove(pickle_path)
            interruptions = self.do_pickle(pickle_path, self.find_beam_interruptions())
        self.Jumps = interruptions['jumps']
        self.Interruptions = interruptions['interruptions']
        return interruptions['interruptions']
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

    def update_all_cut(self):
        self.all_cut = self.generate_all_cut()

    def show_cuts(self, easy=True):
        cuts = self.EasyCutStrings if easy else self.CutStrings
        max_len = max(len(key) for key, value in cuts.iteritems() if str(value))
        for key, value in cuts.iteritems():
            if not key == 'all_cuts' and str(value):
                print '{key}:'.format(key=key.rjust(max_len)), value
        return
