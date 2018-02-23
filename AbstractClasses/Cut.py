import json
from numpy import mean
from Elementary import Elementary
from InfoLegend import InfoLegend
from ROOT import TCut, gROOT, TH1F
from Utils import *
from os import remove


class Cut(Elementary):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, parent_analysis, skip=False):

        if not skip:
            self.analysis = parent_analysis
            self.RunNumber = self.analysis.RunNumber
            Elementary.__init__(self, verbose=self.analysis.verbose)
            self.InfoLegend = InfoLegend(parent_analysis)

            # saving stuff
            self.RootObjects = []

            # config
            self.DUTType = self.load_dut_type()
            self.BeaminterruptionsDir = self.ana_config_parser.get('CUT', 'beaminterruptions_folder')
            self.CutConfig = {}
            self.NCuts = 0

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

    def generate_special_cut(self, excluded=None, included=None, name='special_cut', prnt=True):
        cut = TCut(name, '')
        self.NCuts = 0
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
            self.NCuts += 1
        self.log_info('generated {name} cut with {num} cuts'.format(name=name, num=self.NCuts)) if prnt else do_nothing()
        return cut

    def generate_all_cut(self):
        cut = TCut('all_cuts', '')
        self.NCuts = 0
        for key, value in self.CutStrings.iteritems():
            if not key.startswith('old') and not key.startswith('all_cut'):
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
        dic['spread_low'] = ''
        dic['absMedian_high'] = ''
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
        dic['alignment'] = TCut('alignment', '')
        dic['trigger_phase'] = TCut('trigger_phase', '')
        dic['fiducial'] = TCut('fiducial', '')
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
        cut_string = ''
        if self.CutConfig['EventRange']:
            cut_string = '(event_number<={max}&&event_number>={min})'.format(min=self.CutConfig['EventRange'][0], max=self.CutConfig['EventRange'][1])
        elif self.CutConfig['ExcludeFirst']:
            cut_string = 'event_number>={min}'.format(min=self.CutConfig['ExcludeFirst'])
        return cut_string

    def generate_chi2(self, mode='x', value=None):
        picklepath = self.make_pickle_path('Chi2', run=self.analysis.RunNumber, suf=mode.title())

        def func():
            t = self.log_info('generating chi2 cut in {mod} for run {run}...'.format(run=self.analysis.RunNumber, mod=mode), next_line=False)
            gROOT.SetBatch(1)
            h = TH1F('h', '', 200, 0, 100)
            nq = 100
            chi2s = zeros(nq)
            xq = array([(i + 1) / float(nq) for i in range(nq)])
            self.analysis.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
            h.GetQuantiles(nq, chi2s, xq)
            gROOT.SetBatch(0)
            self.add_info(t)
            return chi2s

        chi2 = self.do_pickle(picklepath, func)
        quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
        assert type(quantile) is int and 0 < quantile <= 100, 'chi2 quantile has to be and integer between 0 and 100'
        cut_value = chi2[quantile] if value is None else value
        string = 'chi2_{mod}<{val}&&chi2_{mod}>=0'.format(val=cut_value, mod=mode)
        return string if quantile > 0 else ''

    def generate_track_angle(self, mode='x'):
        cut_variable = '{t}_{m}'.format(t='slope' if self.analysis.run.has_branch('slope_x') else 'angle', m=mode)
        angles = self.calc_angle(mode)
        string = '{v}>{min}&&{v}<{max}'.format(v=cut_variable, min=angles[mode][0], max=angles[mode][1])
        return string if self.CutConfig['track_angle'] > 0 else ''

    def calc_angle(self, mode='x'):
        picklepath = self.make_pickle_path('TrackAngle', mode, run=self.analysis.lowest_rate_run)

        def func():
            angle = self.CutConfig['track_angle']
            t = self.log_info('Generating angle cut in {m} for run {run} ...'.format(run=self.analysis.RunNumber, m=mode), False)
            self.set_root_output(False)
            draw_str = '{t}_{m}'.format(t='slope' if self.analysis.run.has_branch('slope_x') else 'angle', m=mode)
            n = self.analysis.tree.Draw(draw_str, '{s}>-100'.format(s=draw_str), 'goff')
            mean_ = mean([self.analysis.tree.GetV1()[i] for i in xrange(n)])
            cut_vals = {mode: [mean_ - angle, mean_ + angle]}
            self.add_info(t)
            return cut_vals

        return self.do_pickle(picklepath, func)

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
        self.CutStrings['track_angle_x'] += self.generate_track_angle('x')
        self.CutStrings['track_angle_y'] += self.generate_track_angle('y')
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
        interruptions = self.find_pad_beam_interruptions() if self.DUTType == 'pad' else self.find_pixel_beam_interruptions(show=False)
        dic['interruptions'] = interruptions[0]
        dic['jumps'] = interruptions[1]
        return dic

    def find_pixel_beam_interruptions(self, show=True):
        """ Locates the beam interruptions and cuts some seconds before and some seconds after which are specified in the config file. overlapping segments are fused to one to reduce the size
            of the string. An average of the number of event per time bin is calculated and every bin below 90% or above 120% is excluded """
        # time is in minutes. good results found with bin size of 10 seconds
        bin_size = 10  # seconds
        bins = int(self.analysis.run.totalMinutes * 60 / bin_size)
        self.set_root_output(0)
        h = TH1F('h_beam_time_', 'Beam Interruptions', bins, 0, self.analysis.run.totalMinutes)
        self.analysis.tree.Draw('(time - {off}) / 60000.>>h_beam_time_'.format(off=self.analysis.run.StartTime), '', 'goff')
        mean_ = mean(sorted([h.GetBinContent(i) for i in xrange(h.GetNbinsX())])[-10:])  # only take the ten highest values to get an estimate of the plateau
        interruptions = []
        jumps = []
        n_interr = 0
        ex_range = {key: value / 60. for key, value in self.CutConfig['JumpExcludeRange'].iteritems()}
        for t in xrange(1, h.GetNbinsX() + 1):
            if h.GetBinContent(t) < mean_ * .6 or h.GetBinContent(t) > mean_ * 1.3:
                low_edge = h.GetBinLowEdge(t)
                bin_width = h.GetBinWidth(t)
                if not n_interr or (not low_edge - ex_range['before'] < interruptions[n_interr - 1]['f'] and n_interr):
                    interruptions.append({'i': low_edge - ex_range['before'], 'f': low_edge + bin_width + ex_range['after']})
                    jumps.append([low_edge, low_edge + bin_width])
                    n_interr += 1
                else:
                    interruptions[n_interr - 1]['f'] = low_edge + bin_width + ex_range['after']
                    jumps[n_interr - 1][1] = low_edge + bin_width
        interruptions = [{key: self.analysis.get_event_at_time(value * 60) for key, value in dic.iteritems()} for dic in interruptions]
        jumps = [[self.analysis.get_event_at_time(t * 60) for t in jump] for jump in jumps]
        self.format_histo(h, x_tit='Time [min]', y_tit='Number of Events', y_off=1.7, stats=0, fill_color=self.FillColor)
        if show:
            self.save_histo(h, 'BeamInterruptions', show, lm=.125)
        return interruptions, jumps

    def find_pad_beam_interruptions(self):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        t = self.log_info('Searching for beam interruptions of run {r} ...'.format(r=self.RunNumber), next_line=False)
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
        self.add_info(t)
        return interruptions, jumps

    def __create_jump_ranges(self, jumps):
        ex_range = self.CutConfig['JumpExcludeRange']
        interruptions = []
        time_offset = self.analysis.run.StartTime
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
        ensure_dir(join(self.BeaminterruptionsDir, 'data'))
        pickle_path = self.make_pickle_path('BeamInterruptions', run=self.RunNumber)
        interruptions = self.do_pickle(pickle_path, self.find_beam_interruptions)

        # redo range pickle if config parameters have changed
        ex_range = self.CutConfig['JumpExcludeRange']
        if interruptions['before'] != ex_range['before'] or interruptions['after'] != ex_range['after']:
            remove(pickle_path)
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

    def set_cut(self, name, value=None):
        if name in self.CutStrings:
            self.CutStrings[name].SetTitle('')
            self.CutStrings[name] += value
        else:
            print 'There is no cut with the name "{name}"!'.format(name=name)

    def set_chi2(self, value):
        self.CutConfig['chi2X'] = value
        self.CutConfig['chi2Y'] = value
        self.set_cut('chi2X', self.generate_chi2('x'))
        self.set_cut('chi2Y', self.generate_chi2('y'))

    def update_all_cut(self):
        self.all_cut = self.generate_all_cut()
        self.analysis.AllCuts = self.all_cut

    def show_cuts(self, easy=True):
        cuts = self.EasyCutStrings if easy else self.CutStrings
        max_len = max(len(key) for key, value in cuts.iteritems() if str(value))
        for key, value in cuts.iteritems():
            if not key == 'all_cuts' and str(value):
                print '{key}:'.format(key=key.rjust(max_len)), value
        return

    def get_track_var(self, num, mode):
        if self.analysis.run.has_branch('dia_track_x'):
            return 'dia_track_{m}[{n}]'.format(m=mode, n=num)
        else:
            return 'diam{n}_track_{m}'.format(m=mode, n=num + 1)

    def show_cut_contributions(self):
        contributions = {}
        cut_events = 0
        cuts = TCut('consecutive', '')
        total_events = self.analysis.run.n_entries
        output = OrderedDict()
        for cut in self.CutStrings.values():
            name = cut.GetName()
            if not name.startswith('old') and name != 'all_cuts' and name not in contributions and str(cut):
                cuts += cut
                events = self.analysis.tree.GetEntries('!({0})'.format(str(cuts)))
                output[name] = (1. - float(events) / total_events) * 100.
                events -= cut_events
                print name.rjust(18), '{0:5d} {1:04.1f}%'.format(events, output[name])
                contributions[cut.GetName()] = events
                cut_events += events
