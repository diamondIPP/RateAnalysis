import copy
import os
import numpy as np
import pickle
import json
import ConfigParser
from numpy import mean, array, zeros

from AbstractClasses.Elementary import Elementary
from ROOT import TCut, gROOT, TH1F
from collections import OrderedDict


# todo: make channel cut subclass
class Cut(Elementary):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, parent_analysis, channel, verbose=True):

        self.analysis = parent_analysis
        self.channel = channel
        self._checklist = {"RemoveBeamInterruptions": False,
                           "GenerateCutString": False}
        # saving stuff
        self.save_canvas = None
        self.histo = None

        # config
        self.parser = self.load_parser()
        self.beaminterruptions_folder = self.parser.get('CUT', 'beaminterruptions_folder')
        self.exclude_before_jump = self.parser.getint('CUT', 'excludeBeforeJump')
        self.exclude_after_jump = self.parser.getint('CUT', 'excludeAfterJump')

        # readable cut types
        self.userCutTypes = {'IndividualChCut': '',
                             'EventRange': 'Evts.50k-1433k',
                             'ExcludeFirst': '50k+',
                             'noPulser': '!pulser',
                             'notSaturated': '!saturated',
                             'noBeamInter': '!BeamInter.',
                             'Tracks': 'Track',
                             'peakPos_high': 'peakPos<250',
                             'spread_low': 'spread>20',
                             'absMedian_high': '|median|<10',
                             'pedestalsigma': 'PedSigma5'}
        # default cut Types
        self.cut_types = {'IndividualChCut': '',
                          'EventRange': [],  # [1234, 123456]
                          'ExcludeFirst': 0,  # 50000 events
                          'spread_low': -1,
                          'absMedian_high': -1,
                          'pedestalsigma': -1,
                          'chi2': -1,
                          'track_angle': -1}

        # define cut string dict
        self.cut_strings = OrderedDict()
        self.cut_strings['raw'] = TCut('raw', '')
        self.cut_strings['event_range'] = TCut('event_range', '')
        self.cut_strings['beam_interruptions'] = TCut('beam_interruptions', '')
        self.cut_strings['ped_sigma'] = TCut('ped_sigma', '')
        self.cut_strings['spread_low'] = TCut('spread_low', '')
        self.cut_strings['median'] = TCut('median', '')
        self.cut_strings['pulser'] = TCut('pulser', '')
        self.cut_strings['tracks'] = TCut('tracks', '')
        self.cut_strings['chi2'] = TCut('chi2', '')
        self.cut_strings['track_angle'] = TCut('track_angle', '')
        self.cut_strings['saturated'] = TCut('saturated', '')
        self.cut_strings['bucket'] = TCut('bucket', '')
        self.cut_strings['all_cuts'] = TCut('all_cuts', '')

        self.__cutstring_settings = None
        self.individualCuts = None

        # beam interrupts
        self.jumps = None
        self.jump_ranges = None

        # pedestal sigma (gets loaded with __load_pedestal_data() )
        self.pedestal_mean = None
        self.pedestal_sigma = None

        # miscellaneous
        self.excludefirst = 0
        self.cut = ""

        Elementary.__init__(self, verbose=verbose)

        # generate cut strings
        self.generate_cut_string()
        self.all_cut = self.generate_all_cut()
        self.cut_strings['all_cuts'] = self.all_cut

    def LoadConfig(self):
        # individual additional cuts:
        self.cut = self.parser.get('CUT', 'cut' + str(self.channel)) if not self.parser.get('CUT', 'cut' + str(self.channel)) in ['-1', '', 'True', 'False'] else ''
        self.cut_types['IndividualChCut'] = copy.deepcopy(self.cut)
        self.userCutTypes['IndividualChCut'] = copy.deepcopy(self.cut)

        # general cuts
        self.cut_types['ExcludeFirst'] = self.load_exclude_first(self.parser.getint('CUT', 'excludefirst'))
        self.cut_types['EventRange'] = self.load_event_range(json.loads(self.parser.get('CUT', 'EventRange')))
        self.cut_types['spread_low'] = self.load_spread_low(self.parser.getint('CUT', 'spread_low'))
        self.cut_types['absMedian_high'] = self.load_abs_median_high(self.parser.getint('CUT', 'absMedian_high'))
        self.cut_types['pedestalsigma'] = self.load_pedestal_sigma(self.parser.getint('CUT', 'pedestalsigma'))
        self.cut_types['chi2'] = self.parser.getint('CUT', 'chi2')
        self.cut_types['track_angle'] = self.parser.getint('CUT', 'track_angle')

        # individual cuts
        self.LoadIndividualCuts()

    def generate_all_cut(self):
        cut = TCut('all', '')
        for value in self.cut_strings.values():
            cut += value
        return cut

    def load_parser(self):
        parser = ConfigParser.ConfigParser()
        parser.read('Configuration/AnalysisConfig_' + self.TESTCAMPAIGN + '.cfg')
        return parser

    def LoadIndividualCuts(self):
        path = "Configuration/Individual_Configs/"
        filename = "{testcp}_Run{run}.json".format(testcp=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
        filepath = path + filename
        if os.path.exists(filepath):
            print "Loading run-specific config file:"
            print "\t" + filepath

            f = open(filepath, "r")
            self.individualCuts = json.load(f)
            f.close()
            print "INDIVIDUAL Cuts:"
            print self.individualCuts
            ch = self.channel
            if self.individualCuts[str(ch)]["EventRange"] is not None:
                self.set_event_range([int(self.individualCuts[str(ch)]["EventRange"][0]), int(self.individualCuts[str(ch)]["EventRange"][1])])
            elif self.individualCuts[str(ch)]["ExcludeFirst"] is not None:
                self.set_exclude_first(value=int(self.individualCuts[str(ch)]["ExcludeFirst"]))

            if self.individualCuts[str(ch)]["peakPos_high"] is not None:
                self.set_peakpos_high(value=int(self.individualCuts[str(ch)]["peakPos_high"]))

            if self.individualCuts[str(ch)]["spread_low"] is not None:
                self.set_spread_low(low=int(self.individualCuts[str(ch)]["spread_low"]))

            if self.individualCuts[str(ch)]["absMedian_high"] is not None:
                self.set_abs_median_high(high=int(self.individualCuts[str(ch)]["absMedian_high"]))

    # ==============================================
    # region GET & SET

    def load_event_range(self, event_range=None):
        """
        Gets the event range cut. If the arguments are negative, they are interpreted as time in minutes. Therefore, e.g.
        SetEventRange(-10, 700000) means that only events are considered, which fulfill: >10 minutes after run start event number < 700000
        :param event_range:
        :return:
        """
        if not event_range:
            event_range = [0, 0]
        for i, value in enumerate(event_range):
            if value < 0:
                event_range[i] = self.analysis.get_event_at_time(time_sec=-1 * value * 60)
        if event_range[0] and event_range[1]:
            pass
        elif event_range[0]:
            event_range[1] = self.analysis.get_event_at_time(-1)
        elif event_range[1]:
            event_range[0] = self.excludefirst
        else:
            event_range = []
        return event_range

    def set_event_range(self, event_range):
        self.cut_types['EventRange'] = self.load_event_range(event_range)

    def GetIncludedEvents(self, maxevent=None):
        """
        Returns a list of event numbers, which are not excluded by the
        following cuts (i.e. event cuts):
            - excludeFirst cutr
            - event-range cut
            - beam-interruptions cut
        :return: list of included event numbers
        """
        minevent = self.GetMinEvent()
        if maxevent is None:
            maxevent = self.GetMaxEvent()

        excluded = [i for i in np.arange(0, minevent)]  # first events
        if self.cut_types["noBeamInter"]:
            self.GetBeamInterruptions()
            for i in xrange(len(self.jump_ranges["start"])):
                excluded += [i for i in np.arange(self.jump_ranges["start"][i], self.jump_ranges["stop"][i] + 1)]  # events around jumps
        excluded.sort()
        all_events = np.arange(0, maxevent)
        included = np.delete(all_events, excluded)
        return included

    def load_exclude_first(self, value):
        """
        Sets how many events at the very beginning of the run should be excluded. if the argument is negative, it will be interpreted as time in minutes. For a positive argument it is interpreted as
        maximum event number.
        :param value: events or time in minutes
        :return:
        """
        if value > 0:
            self.userCutTypes['ExcludeFirst'] = str(int(value) / 1000) + 'k+'
            return value
        elif value == 0:
            self.userCutTypes['ExcludeFirst'] = ''
            return 0
        else:
            self.userCutTypes['ExcludeFirst'] = str(-1 * value) + 'min+'
            seconds = -1 * value * 60
            event = self.analysis.get_event_at_time(seconds)
            return event

    def set_exclude_first(self, value):
        self.cut_types['ExcludeFirst'] = self.load_exclude_first(value)

    def load_peakpos_high(self, high):
        if high > 0:
            self.userCutTypes['peakPos_high'] = 'peakPos<{high}'.format(high=high)
            return high
        else:
            return -1

    def set_peakpos_high(self, value):
        self.cut_types['peakPos_high'] = self.load_peakpos_high(value)

    def load_spread_low(self, value):
        if value > 0:
            self.userCutTypes['spread_low'] = 'spread>{low}'.format(low=value)
            return value
        else:
            return -1

    def set_spread_low(self, low):
        self.cut_types['spread_low'] = self.load_spread_low(low)

    def load_abs_median_high(self, value):
        if value > 0:
            self.userCutTypes['absMedian_high'] = '|median|<{high}'.format(high=value)
            return value
        else:
            return -1

    def set_abs_median_high(self, high):
        self.cut_types['absMedian_high'] = self.load_abs_median_high(high)

    def load_pedestal_sigma(self, value):
        if value > 0:
            self.userCutTypes['pedestalsigma'] = 'PedSigma' + str(value)
            return value
        else:
            self.userCutTypes['pedestalsigma'] = ''
            return -1

    def set_pedestal_sigma(self, sigma=-1):
        self.cut_types['pedestalsigma'] = self.load_pedestal_sigma(sigma)

    # endregion

    def GetEventRange(self):
        """
        Returns a the lowest and highest event numbers to consider in the analysis. This event-range cut is defined either in the Analysis config file or in the indvidual cut file.
        The returned object is a list, which is empty if no event-range cut is applied.
        :return: cut eventrange
        """
        return self.cut_types["EventRange"]

    def GetMinEvent(self):
        """
        Returns the smallest event number satisfying the cut conditions.
        :return:
        """
        if self.cut_types["EventRange"]:
            return self.cut_types["EventRange"][0]
        elif self.cut_types["ExcludeFirst"] > 0:
            return self.cut_types["ExcludeFirst"]
        else:
            return 0

    def GetNEvents(self):
        """
        Returns the number of events satisfying the cut conditions.
        :return:
        """
        totEvents = self.analysis.get_event_at_time(-1)
        if self.cut_types["EventRange"]:
            return self.cut_types["EventRange"][1] - self.cut_types["EventRange"][0]
        elif self.cut_types["ExcludeFirst"] > 0:
            return totEvents - self.cut_types["ExcludeFirst"]
        else:
            return totEvents

    def GetMaxEvent(self):
        """
        Returns the highest event number satisfying the cut conditions.
        :return:
        """
        totEvents = self.analysis.get_event_at_time(-1)
        if self.cut_types["EventRange"]:
            return self.cut_types["EventRange"][1]
        else:
            return totEvents

    # ==============================================
    # region GENERATE CUT STRINGS
    def generate_event_range(self):
        if self.cut_types['EventRange']:
            self.cut_strings['event_range'] += '(event_number<={max}&&event_number>={min})'.format(min=self.cut_types['EventRange'][0], max=self.cut_types['EventRange'][1])
        elif self.cut_types['ExcludeFirst']:
            self.cut_strings['event_range'] += 'event_number>={min}'.format(min=self.cut_types['ExcludeFirst'])

    def generate_chi2(self):
        gROOT.SetBatch(1)
        h = TH1F('h', '', 200, 0, 100)
        nq = 100
        yq = zeros(nq)
        xq = array([(i + 1) / float(nq) for i in range(nq)])
        self.analysis.tree.Draw('chi2_tracks>>h', '', 'goff')
        h.GetQuantiles(nq, yq, xq)
        string = 'chi2_tracks<{val}&&chi2_tracks>=0'.format(val=yq[self.cut_types['chi2']])
        gROOT.SetBatch(0)
        return TCut(string) if self.cut_types['chi2'] > 0 else ''

    def generate_slope(self):
        # fit the slope to get the mean
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        h_x = TH1F('hx', '', 70, -4, 4)
        h_y = TH1F('hy', '', 70, -4, 4)
        self.analysis.tree.Draw('slope_x>>hx', '', 'goff')
        self.analysis.tree.Draw('slope_y>>hy', '', 'goff')
        angle = self.cut_types['track_angle']
        fit_result = h_x.Fit('gaus', 'qs')
        x_mean = fit_result.Parameters()[1]
        x = [x_mean - angle, x_mean + angle]
        fit_result = h_y.Fit('gaus', 'qs')
        y_mean = fit_result.Parameters()[1]
        y = [y_mean - angle, y_mean + angle]
        c = gROOT.FindObject('c1')
        c.Close()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        # create the cut string
        string = 'slope_x>{minx}&&slope_x<{maxx}&&slope_y>{miny}&&slope_y<{maxy}'.format(minx=x[0], maxx=x[1], miny=y[0], maxy=y[1])
        return TCut(string) if angle > 0 else ''

    def generate_bucket(self):
        num = self.analysis.get_signal_numbers('e', 2)
        name = self.analysis.get_signal_names(num)[self.channel]
        string = '{sig2}-{sig1}==0'.format(sig2=name, sig1=self.analysis.signal_names[self.channel])
        return TCut(string)

    def generate_cut_string(self, gen_PulserCut=True, gen_EventRange=True, gen_ExcludeFirst=True, setChannel=True):
        """
        Creates the cut string, which will be stored in self.cut. With the arguments set to False, different cut types can be deactivated in the cut string.
        :param gen_PulserCut:
        :param gen_EventRange:
        :param gen_ExcludeFirst:
        :return:
        """
        gROOT.SetBatch(1)
        if self._checklist["GenerateCutString"]:
            self.LoadConfig()  # re-generate
        cutstring = self.cut

        # --CHI2 --
        self.cut_strings['chi2'] = self.generate_chi2()
        self.cut_strings['track_angle'] = self.generate_slope()
        self.cut_strings['bucket'] = self.generate_bucket()

        # -- EVENT RANGE CUT --
        self.generate_event_range()
        if self.cut_types["EventRange"] != [] and gen_EventRange:
            if cutstring != "":
                cutstring += "&&"
            cutstring += "(event_number<={maxevent}&&event_number>={minevent})".format(minevent=self.cut_types["EventRange"][0], maxevent=self.cut_types["EventRange"][1])
            self.userCutTypes["EventRange"] = "Evts.{min}k-{max}k".format(min=int(self.cut_types["EventRange"][0]) / 1000, max=int(self.cut_types["EventRange"][1]) / 1000)
            self.userCutTypes["ExcludeFirst"] = ""
        elif self.cut_types["ExcludeFirst"] > 0 and gen_ExcludeFirst:
            if cutstring != "":
                cutstring += "&&"
            cutstring += "event_number>={minevent}".format(minevent=self.cut_types["ExcludeFirst"])
            self.userCutTypes["ExcludeFirst"] = "Evts.{min}k+".format(min=int(self.cut_types["ExcludeFirst"]) / 1000)
            self.userCutTypes["EventRange"] = ""
        else:
            self.userCutTypes["EventRange"] = ""
            self.userCutTypes["ExcludeFirst"] = ""

        # -- PULSER CUT --
        self.cut_strings['pulser'] += '!pulser'

        # -- SATURATED CUT --
        self.cut_strings['saturated'] = '!is_saturated[{ch}]'.format(ch=self.channel)

        # -- TRACK CUT --
        self.cut_strings['tracks'] += 'n_tracks'

        # -- SPREAD LOW CUT --
        if self.cut_types['spread_low'] > 0:
            self.cut_strings['spread_low'] += 'sig_spread[{ch}]>{low}'.format(ch=self.channel, low=int(self.cut_types['spread_low']))
            if cutstring != '':
                cutstring += '&&'
            cutstring += 'sig_spread[{channel}]>{low}'.format(channel=self.channel, low=int(self.cut_types['spread_low']))
            self.userCutTypes['spread_low'] = 'spread>{low}'.format(low=int(self.cut_types['spread_low']))
        else:
            self.userCutTypes['spread_low'] = ''

        # -- MEDIAN CUT --
        if self.cut_types['absMedian_high'] > 0:
            self.cut_strings['median'] += 'abs(median[{ch}])<{high}'.format(ch=self.channel, high=int(self.cut_types['absMedian_high']))
            if cutstring != '':
                cutstring += '&&'
            cutstring += 'abs(median[{channel}])<{high}'.format(channel=self.channel, high=int(self.cut_types['absMedian_high']))
            self.userCutTypes['absMedian_high'] = '|median|<{high}'.format(high=int(self.cut_types['absMedian_high']))
        else:
            self.userCutTypes['absMedian_high'] = ''

        # -- PEDESTAL SIGMA CUT --
        if self.cut_types['pedestalsigma'] > 0:
            if cutstring != '':
                cutstring += '&&'
            self.__load_pedestal_data()
            ped_range = self.__calc_pedestal_range()
            string = '{ped}>{min}&&{ped}<{max}'.format(ped=self.analysis.pedestal_names[self.channel], min=ped_range[0], max=ped_range[1])
            cutstring += string
            self.userCutTypes["pedestalsigma"] = "PedSigma" + str(self.cut_types['pedestalsigma'])
            self.cut_strings['ped_sigma'] += string
        else:
            self.userCutTypes["pedestalsigma"] = ""

        # -- set the channel on the cuts --
        if setChannel:
            self.cut = cutstring
            self.cut = self.cut.format(channel=self.channel)

        # -- BEAM INTERRUPTION CUT --
        if self._checklist["GenerateCutString"]:
            self.__generate_beam_interruptions(justDoIt=True)
            self.userCutTypes["noBeamInter"] = "beamOn"
        else:
            self.__generate_beam_interruptions()
            self.userCutTypes["noBeamInter"] = "BeamOn"

        self._checklist["GenerateCutString"] = True
        self.__cutstring_settings = {
            "gen_PulserCut": gen_PulserCut,
            "gen_EventRange": gen_EventRange,
            "gen_ExcludeFirst": gen_ExcludeFirst
        }
        gROOT.SetBatch(0)

    def __calc_pedestal_range(self):
        sigma = self.pedestal_sigma
        ped_mean = self.pedestal_mean
        sigma_range = self.cut_types['pedestalsigma']
        return [ped_mean - sigma_range * sigma, ped_mean + sigma_range * sigma]

    def __generate_beam_interruptions(self, justDoIt=False):
        """
        This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied.
        :return: cut
        """
        if not self._checklist["RemoveBeamInterruptions"] or justDoIt:
            self.GetBeamInterruptions()

            njumps = len(self.jump_ranges["start"])
            cut_string = ''
            for i in xrange(njumps):
                string = "!(event_number<={upper}&&event_number>={lower})".format(upper=self.jump_ranges["stop"][i], lower=self.jump_ranges["start"][i])
                if self.cut != "":
                    self.cut += "&&"
                self.cut += string
                # new seperate strings
                if cut_string != '':
                    cut_string += '&&'
                cut_string += string
            self.cut_strings['beam_interruptions'] += cut_string
            self._checklist["RemoveBeamInterruptions"] = True

        return self.cut

    # endregion

    def _checkCutStringSettings(self, gen_PulserCut, gen_EventRange, gen_ExcludeFirst):
        if self.__cutstring_settings["gen_PulserCut"] == gen_PulserCut and self.__cutstring_settings["gen_EventRange"] == gen_EventRange \
                and self.__cutstring_settings["gen_ExcludeFirst"] == gen_ExcludeFirst:
            return True
        else:
            return False

    def __load_pedestal_data(self):
        picklepath = 'Configuration/Individual_Configs/PedestalPeak/{tc}_{run}_{ch}_PedestalPeak.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number, ch=self.channel)
        if self.pedestal_mean is None:
            # Loading pedestal peak fit data from pickle
            if os.path.exists(picklepath):
                picklefile = open(picklepath, 'rb')
                fitparameters = pickle.load(picklefile)
                picklefile.close()
            # fitting
            else:
                gROOT.SetBatch(1)
                pedestalhisto = TH1F('tmphisto', 'bla', 400, -100, 100)
                self.analysis.tree.Draw(self.analysis.pedestal_names[self.channel] + '>>tmphisto')
                ped_peakpos = pedestalhisto.GetBinCenter(pedestalhisto.GetMaximumBin())
                bin1 = pedestalhisto.FindFirstBinAbove(pedestalhisto.GetMaximum() / 2)
                bin2 = pedestalhisto.FindLastBinAbove(pedestalhisto.GetMaximum() / 2)
                fwhm = pedestalhisto.GetBinCenter(bin2) - pedestalhisto.GetBinCenter(bin1)
                print 'FWHM:', fwhm
                pedestalhisto.Fit('gaus', '', '', ped_peakpos - fwhm / 2, ped_peakpos + fwhm / 2)
                fitfunc = pedestalhisto.GetFunction('gaus')
                fitparameters = [fitfunc.GetParameter(0), fitfunc.GetParameter(1), fitfunc.GetParameter(2)]
                gROOT.SetBatch(0)
                # save to file
                picklefile = open(picklepath, 'wb')
                pickle.dump(fitparameters, picklefile)
                picklefile.close()
            self.pedestal_mean = fitparameters[1]
            self.pedestal_sigma = fitparameters[2]

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
        self.jumps = interrupts
        self.__save_beaminterrupts()
        self.__create_jump_ranges()
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
        return self.jump_ranges

    def GetBeamInterruptions(self):
        """
        If beam interruption data exist in beaminterruptions/data/, it will load it in order to account for beam interruptions. The data is stored as a list of jumps, dumped into a pickle file.
        If no pickle file exists, it will perform a beam interruption analysis in order to identify the beam interruptions. The found interruptions are stored in a list at .jumps and dumped into
        a pickle file.
        :return: list of events where beam interruptions occures
        """
        if self.jump_ranges is None:
            picklepath = self.beaminterruptions_folder + "/data/{testcampaign}Run_{run}.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
            if os.path.exists(picklepath):
                # print "Loading beam interruption data from pickle file: \n\t"+picklepath
                jumpfile = open(picklepath, "rb")
                self.jumps = pickle.load(jumpfile)
                self.__create_jump_ranges()
                jumpfile.close()
            else:
                print "No pickle file found at: ", picklepath, "\n .. analyzing beam interruptions.. "
                self.find_beam_interruptions()
        return self.jumps

    def GetCutFunctionDef(self):
        """
        Returns the cut function definition, which is of type string.
        It is used for applying the cut on event-by-event readout
        loops over a TTree.
        The cut function definition defines the constraints on the
        branches in the ROOT TTree.

        Example: (considering only the branches 'pulser' and
        'is_saturated')
            A possible cut function definition could be:
            'lambda pulser, is_saturated: not pulser and not is_saturated'
        Before the iteration over the TTree branches are executed, a
        boolean function is defined as:
            check = exec(cut_object.GetCutFunctionDef())
        The event selection is then based on this lambda function
        'check':
            for i in xrange(nentries):
                tree.GetEntry(i)
                signal = tree.signal
                pulser = tree.pulser
                saturated = tree.is_saturated
                if check(pulser, saturated):
                    someting.Fill(signal)

        Note: the cut on the event-basis has to be applied by the
        method .GetIncludedEvents(), which returns only the events
        satisfying the event cuts.
        :return:
        """
        defstring_ = "lambda pulser, is_saturated, n_tracks, fft_mean, INVfft_max, sig_time, sig_spread, median: "
        def_ = ""

        if self.cut_types["IndividualChCut"] != "":
            raw_input("WARNING: cut0 and cut3 cannot be applied on Tracking Data! Press Enter to Continue..")

        if self.cut_types["noPulser"] == 1:
            def_ += "not pulser"
        elif self.cut_types["noPulser"] == 0:
            def_ += "pulser"

        if self.cut_types["notSaturated"]:
            if def_ != "":
                def_ += " and "
            def_ += "not is_saturated"

        if self.cut_types["Tracks"]:
            if def_ != "":
                def_ += " and "
            def_ += "n_tracks"

        if self.cut_types["FFT"]:
            if def_ != "":
                def_ += " and "
            assert False, "FFT cut not yet implemented in GetCutFunctionDef() method of Cut class. "
            # to do: FFT entry in _cutTypes should be dict and/or contain a TCutG instance

        if self.cut_types["peakPos_high"] > 0:
            if def_ != "":
                def_ += " and "
            def_ += "sig_time<{high}".format(high=int(self.cut_types["peakPos_high"]))

        if self.cut_types["spread_low"] > 0:
            if def_ != "":
                def_ += " and "
            def_ += "sig_spread>{low}".format(low=int(self.cut_types["spread_low"]))

        if self.cut_types["absMedian_high"] > 0:
            if def_ != "":
                def_ += " and "
            def_ += "abs(median)>{low}".format(low=int(self.cut_types["absMedian_high"]))

        return defstring_ + def_

    def AddCutString(self, cutstring):
        pass

    def GetCut(self, setChannel=True, gen_PulserCut=True, gen_EventRange=True, gen_ExcludeFirst=True, ):
        """
        Returns the cut string.
        If needed, it will re-generate the cut string.
        :param gen_PulserCut:
        :param gen_EventRange:
        :param gen_ExcludeFirst:
        :return:
        """
        # channel = self.channel
        if not self._checkCutStringSettings(gen_PulserCut, gen_EventRange, gen_ExcludeFirst):
            self.generate_cut_string(gen_PulserCut, gen_EventRange, gen_ExcludeFirst, setChannel=setChannel)
        return self.cut

    def GetUserCutString(self):
        """
        Returns a short, more user-friendly cut string, which can be
        used to display the cut configuration as terminal prompt or
        inside a canvas.
        :return:
        """
        string_ = ""
        for type_ in self.userCutTypes.keys():
            if self.userCutTypes[type_] != "":
                string_ += self.userCutTypes[type_] + ", "
        if string_ != "":
            string_ = string_[:-2]
        return string_

    def show_cuts(self):
        for key, value in self.cut_strings.iteritems():
            print key, value
        return

    def SetFFTCut(self):
        pass
