import json
from numpy import zeros
from Elementary import Elementary
from InfoLegend import InfoLegend
from ROOT import TCut, gROOT, TH1F, TPie
from Utils import *
from Plots import Plots


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
            self.Plots = Plots(self.analysis.run)

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
        dic['aligned'] = TCut('aligned', '')
        dic['trigger_phase'] = TCut('trigger_phase', '')
        dic['tracks'] = TCut('tracks', '')
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
        dic['hit'] = TCut('hit', '')
        dic['masks'] = TCut('masks', '')
        dic['chi2X'] = TCut('chi2X', '')
        dic['chi2Y'] = TCut('chi2Y', '')
        dic['slope_x'] = TCut('slope_x', '')
        dic['slope_y'] = TCut('slope_y', '')
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
        self.CutConfig['JumpExcludeRange'] = json.loads(self.ana_config_parser.get('CUT', 'exclude_around_jump'))
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(self.ana_config_parser.getfloat('CUT', 'excludefirst'))
        self.CutConfig['EventRange'] = self.load_event_range(json.loads(self.ana_config_parser.get('CUT', 'EventRange')))
        self.CutConfig['chi2X'] = self.ana_config_parser.getint('CUT', 'chi2X')
        self.CutConfig['chi2Y'] = self.ana_config_parser.getint('CUT', 'chi2Y')
        self.CutConfig['slope'] = self.ana_config_parser.getint('CUT', 'slope')

    def load_event_range(self, event_range=None):
        """ Gets the event range cut. If the arguments are negative, they are interpreted as time in minutes. Therefore, e.g. load_event_range(-10, 700000) means that only events are considered
        which fulfill: >10 minutes after run start event number < 700000 """
        if event_range is None:
            event_range = [0, 0]
        for i, value in enumerate(event_range):
            if value < 0:
                event_range[i] = self.analysis.get_event_at_time(seconds=-1 * value * 60)
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
        return value if value >= 0 else self.analysis.get_event_at_time(-1 * value * 60)

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

        chi2 = do_pickle(picklepath, func)
        quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
        assert type(quantile) is int and 0 < quantile <= 100, 'chi2 quantile has to be and integer between 0 and 100'
        cut_value = chi2[quantile] if value is None else value
        string = 'chi2_{mod}<{val}&&chi2_{mod}>=0'.format(val=cut_value, mod=mode)
        return string if quantile > 0 else ''

    def generate_slope(self, mode='x'):
        cut_variable = '{t}_{m}'.format(t='slope' if self.analysis.run.has_branch('slope_x') else 'angle', m=mode)
        angles = self.calc_angle(mode)
        string = '{v}>{min}&&{v}<{max}'.format(v=cut_variable, min=angles[mode][0], max=angles[mode][1])
        return string if self.CutConfig['slope'] > 0 else ''

    def calc_angle(self, mode='x'):
        picklepath = self.make_pickle_path('TrackAngle', mode, run=self.analysis.lowest_rate_run)

        def func():
            angle = self.CutConfig['slope']
            t = self.log_info('Generating angle cut in {m} for run {run} ...'.format(run=self.analysis.RunNumber, m=mode), False)
            set_root_output(False)
            draw_str = '{t}_{m}'.format(t='slope' if self.analysis.run.has_branch('slope_x') else 'angle', m=mode)
            n = self.analysis.tree.Draw(draw_str, '{s}>-100'.format(s=draw_str), 'goff')
            mean_ = mean([self.analysis.tree.GetV1()[i] for i in xrange(n)])
            cut_vals = {mode: [mean_ - angle, mean_ + angle]}
            self.add_info(t)
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
        return self.find_pad_beam_interruptions() if self.DUTType == 'pad' else self.find_pixel_beam_interruptions(show=False)

    def find_pixel_beam_interruptions(self, bin_width=10, rel_time=True, show=True):
        """ Locates the beam interruptions and cuts some seconds before and some seconds after which are specified in the config file. overlapping segments are fused to one to reduce the size
            of the string. An average of the number of event per time bin is calculated and every bin below 90% or above 120% is excluded """
        # time is in minutes. good results found with bin size of 10 seconds
        set_root_output(False)
        h = TH1F('h_beam_time_', 'Beam Interruptions', *self.Plots.get_time_binning(bin_width))
        self.analysis.tree.Draw('time / 1000. >> h_beam_time_'.format(off=self.analysis.run.StartTime), '', 'goff')
        mean_ = mean(sorted([h.GetBinContent(i) for i in xrange(h.GetNbinsX())])[-20:-10])  # only take the ten highest values to get an estimate of the plateau
        jumps = []
        tup = [0, 0]
        cut = .4
        last_rate = 0
        for ibin in xrange(1, h.GetNbinsX() + 1):
            rate = abs(1 - h.GetBinContent(ibin) / mean_)  # normalised deviation from the mean
            if rate > cut:
                tup[0] = h.GetBinLowEdge(ibin)
            elif rate < cut < last_rate:
                tup[1] = h.GetBinLowEdge(ibin) + bin_width
                jumps.append(tup)
                tup = [0, 0]
            last_rate = rate
        if tup[0] != tup[1]:  # if rate did not went down before the run stopped
            jumps.append([tup[0], self.analysis.run.EndTime])
        jumps = [[self.analysis.get_event_at_time(t - self.analysis.run.StartTime) for t in jump] for jump in jumps]
        interruptions = self.__create_jump_ranges(jumps)
        self.format_histo(h, x_tit='Time [min]', y_tit='Number of Events', y_off=1.7, stats=0, fill_color=self.FillColor, t_ax_off=self.analysis.run.StartTime if rel_time else 0)
        if show:
            self.save_histo(h, 'BeamInterruptions', show, lm=.125)
        return jumps, interruptions

    def find_pad_beam_interruptions(self):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        t = self.log_info('Searching for beam interruptions of run {r} ...'.format(r=self.RunNumber), next_line=False)
        bin_width = 200
        rates = [self.analysis.run.tree.Draw('1', 'pulser', 'goff', bin_width, i * bin_width) / float(bin_width) for i in xrange(self.analysis.run.n_entries / bin_width)]
        jumps = []
        tup = [0, 0]
        cut = .4  # if rate goes higher than n %
        for i, rate in enumerate(rates):
            if rate >= cut > rates[i - 1]:
                tup[0] = int((i + .5) * bin_width)
            elif rate < cut <= rates[i - 1]:
                tup[1] = int((i + .5) * bin_width)
                jumps.append(tup)
                tup = [0, 0]
        if tup[0] != tup[1]:  # if rate did not went down before the run stopped
            jumps.append([tup[0], self.analysis.run.n_entries - 1])
        interruptions = self.__create_jump_ranges(jumps)
        self.add_info(t)
        return jumps, interruptions

    def __create_jump_ranges(self, jumps):
        interruptions = []
        i = 0
        for tup in jumps:
            t_start = self.analysis.run.get_time_at_event(tup[0]) - self.analysis.run.StartTime - self.CutConfig['JumpExcludeRange'][0]
            t_stop = self.analysis.run.get_time_at_event(tup[1]) - self.analysis.run.StartTime + self.CutConfig['JumpExcludeRange'][1]
            # if interruptions overlay just set the last stop to the current stop
            if i and t_start <= (interruptions[i - 1][1]) + 10:
                interruptions[i - 1][1] = t_stop
                continue
            interruptions.append([t_start, t_stop])
            i += 1
        return [[self.analysis.run.get_event_at_time(t) for t in tup] for tup in interruptions]

    def get_beam_interruptions(self):
        """
        If beam interruption data exist in beaminterruptions/data/, it will load it in order to account for beam interruptions. The data is stored as a list of jumps, dumped into a pickle file.
        If no pickle file exists, it will perform a beam interruption analysis in order to identify the beam interruptions. The found interruptions are stored in a list at .jumps and dumped into
        a pickle file.
        """
        # check if directories exist
        ensure_dir(self.BeaminterruptionsDir)
        ensure_dir(join(self.BeaminterruptionsDir, 'data'))
        pickle_path = self.make_pickle_path('BeamInterruptions', run=self.RunNumber, suf='_'.join(str(i) for i in self.CutConfig['JumpExcludeRange']))
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
            return 'dia_track_{m}_local[{n}]'.format(m=mode, n=num)
        else:
            return 'diam{n}_track_{m}'.format(m=mode, n=num + 1)

    def generate_consecutive_cuts(self):
        pass

    def draw_cut_contributions(self, flat=False, short=False, show=True):
        set_root_output(show)
        contr = OrderedDict()
        n_events = self.analysis.run.n_entries
        cut_events = 0
        for i, (key, cut) in enumerate(self.generate_consecutive_cuts().iteritems()):
            events = n_events - int(self.analysis.tree.Draw('1', cut, 'goff'))
            print key.rjust(18), '{0:5d} {1:04.1f}%'.format(events - cut_events, (1. - float(events) / n_events) * 100.)
            contr[key.title().replace('_', ' ')] = (events - cut_events, self.get_color())
            cut_events = events
        contr['Good Events'] = (n_events - cut_events, self.get_color())
        print contr
        sorted_contr = OrderedDict(sorted(OrderedDict(item for item in contr.iteritems() if item[1][0] >= (.03 * n_events if short else 0)).iteritems(), key=lambda x: x[1]))  # sort by size
        sorted_contr.update({'Other': (n_events - sum(v[0] for v in sorted_contr.values()), self.get_color())} if short else {})
        sorted_contr = OrderedDict(sorted_contr.popitem(0 if i % 2 else -1) for i in xrange(len(sorted_contr)))  # sort by largest->smallest->next largest...
        print sorted_contr
        pie = TPie('pie', 'Cut Contributions', len(sorted_contr), array([v[0] for v in sorted_contr.values()], 'f'), array([v[1] for v in sorted_contr.values()], 'i'))
        for i, label in enumerate(sorted_contr.iterkeys()):
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label)
        self.format_pie(pie, h=.04, r=.2, text_size=.025, angle3d=70, label_format='%txt (%perc)', angle_off=250)
        self.save_histo(pie, draw_opt='{0}rsc'.format('3d' if not flat else ''), show=show)
        self.reset_colors()
        return sorted_contr
