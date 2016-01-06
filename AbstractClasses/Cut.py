import os
import pickle
import json
import ConfigParser
from numpy import mean, array, zeros, arange, delete
from AbstractClasses.Elementary import Elementary
from ROOT import TCut, gROOT, TH1F
from collections import OrderedDict
from Extrema import Extrema2D


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

        # readable cut types
        self.EasyCutStrings = {'IndividualChCut': '',
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

        # config
        self.parser = self.load_parser()
        self.beaminterruptions_folder = self.parser.get('CUT', 'beaminterruptions_folder')
        self.exclude_before_jump = self.parser.getint('CUT', 'excludeBeforeJump')
        self.exclude_after_jump = self.parser.getint('CUT', 'excludeAfterJump')
        self.CutConfig = {}

        self.IndividualCuts = None
        self.load_individual_cuts()

        # define cut string dict
        self.cut_strings = OrderedDict()
        self.cut_strings['raw'] = TCut('raw', '')
        self.cut_strings['pulser'] = TCut('pulser', '')
        self.cut_strings['event_range'] = TCut('event_range', '')
        self.cut_strings['beam_interruptions'] = TCut('beam_interruptions', '')
        self.cut_strings['ped_sigma'] = TCut('ped_sigma', '')
        self.cut_strings['spread_low'] = TCut('spread_low', '')
        self.cut_strings['median'] = TCut('median', '')
        self.cut_strings['tracks'] = TCut('tracks', '')
        self.cut_strings['chi2'] = TCut('chi2', '')
        self.cut_strings['track_angle'] = TCut('track_angle', '')
        self.cut_strings['saturated'] = TCut('saturated', '')
        self.cut_strings['bucket'] = TCut('bucket', '')
        self.cut_strings['all_cuts'] = TCut('all_cuts', '')

        self.region_cut = TCut('region_cut', '')

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

    def load_config(self):
        self.CutConfig['IndividualChCut'] = ''
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(self.parser.getint('CUT', 'excludefirst'))
        self.CutConfig['EventRange'] = self.load_event_range(json.loads(self.parser.get('CUT', 'EventRange')))
        self.CutConfig['spread_low'] = self.load_spread_low(self.parser.getint('CUT', 'spread_low'))
        self.CutConfig['absMedian_high'] = self.load_abs_median_high(self.parser.getint('CUT', 'absMedian_high'))
        self.CutConfig['pedestalsigma'] = self.load_pedestal_sigma(self.parser.getint('CUT', 'pedestalsigma'))
        self.CutConfig['chi2'] = self.parser.getint('CUT', 'chi2')
        self.CutConfig['track_angle'] = self.parser.getint('CUT', 'track_angle')

    def generate_all_cut(self):
        cut = TCut('all_cuts', '')
        for value in self.cut_strings.values():
            cut += value
        return cut

    def load_parser(self):
        parser = ConfigParser.ConfigParser()
        parser.read('Configuration/AnalysisConfig_' + self.analysis.TESTCAMPAIGN + '.cfg')
        return parser

    def load_individual_cuts(self):
        path = 'Configuration/Individual_Configs/'
        filename = '{tc}_Run{run}.json'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
        filepath = path + filename
        if os.path.exists(filepath):
            print 'Loading run-specific config file:', filepath
            f = open(filepath, 'r')
            self.IndividualCuts = json.load(f)
            f.close()
            print 'INDIVIDUAL Cuts:\n', self.IndividualCuts
            ch = str(self.channel)
            if self.IndividualCuts[ch]['ExcludeFirst'] is not None:
                self.set_exclude_first(value=int(self.IndividualCuts[str(ch)]['ExcludeFirst']))
            if self.IndividualCuts[ch]['EventRange'] is not None:
                self.set_event_range(self.IndividualCuts[str(ch)]['EventRange'])
            if self.IndividualCuts[ch]['peakPos_high'] is not None:
                self.set_peakpos_high(value=int(self.IndividualCuts[str(ch)]['peakPos_high']))
            if self.IndividualCuts[ch]['spread_low'] is not None:
                self.set_spread_low(low=int(self.IndividualCuts[str(ch)]['spread_low']))
            if self.IndividualCuts[ch]['absMedian_high'] is not None:
                self.set_abs_median_high(high=int(self.IndividualCuts[str(ch)]['absMedian_high']))

    # ==============================================
    # region GET & SET

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
    def generate_region(self, signal_histo, mean_histo):
        extrema = Extrema2D(signal_histo, mean_histo)
        extrema.region_scan()
        extrema.show_voting_histos()
        all_string = ''
        nr = 2 if self.channel else 1
        for col in xrange(extrema.cols):
            all_val = [bool(extrema.VotingHistos['max'].GetBinContent(col, row)) for row in xrange(extrema.rows)]
            # print col, all_val
            if True not in all_val:
                continue
            all_string += '||' if all_string else ''
            xmin = extrema.VotingHistos['max'].GetXaxis().GetBinLowEdge(col)
            xmax = extrema.VotingHistos['max'].GetXaxis().GetBinUpEdge(col)
            all_string += '(diam{nr}_track_x>{xmin}&&diam{nr}_track_x<{xmax})&&'.format(nr=nr, xmin=xmin, xmax=xmax)
            y_string = ''
            cont = True
            for row in xrange(extrema.rows + 1):
                val = extrema.VotingHistos['max'].GetBinContent(col, row) if not row == extrema.rows else 0
                last_val = extrema.VotingHistos['max'].GetBinContent(col, row - 1) if row else 0
                if val and not last_val:
                    y = extrema.VotingHistos['max'].GetYaxis().GetBinLowEdge(row)
                    if y < abs(1e-10):
                        cont = False
                        continue
                    cont = True
                    y_string += '||' if y_string else '('
                    y_string += 'diam{nr}_track_y>{y}&&'.format(nr=nr, y=y)
                elif not val and last_val and cont:
                    y_string += 'diam{nr}_track_y<{y}'.format(nr=nr, y=extrema.VotingHistos['max'].GetYaxis().GetBinUpEdge(row))
            y_string += ')'
            all_string += y_string
        self.region_cut += all_string
        return extrema
        # print all_string

    def generate_event_range(self):
        if self.CutConfig['EventRange']:
            self.cut_strings['event_range'] += '(event_number<={max}&&event_number>={min})'.format(min=self.CutConfig['EventRange'][0], max=self.CutConfig['EventRange'][1])
        elif self.CutConfig['ExcludeFirst']:
            self.cut_strings['event_range'] += 'event_number>={min}'.format(min=self.CutConfig['ExcludeFirst'])

    def generate_chi2(self):
        picklepath = 'Configuration/Individual_Configs/Chi2/{tc}_{run}_{ch}_Chi2.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number, ch=self.channel)

        def func():
            print 'generating chi2 cut for ch{ch}...'.format(ch=self.channel)
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
        picklepath = 'Configuration/Individual_Configs/Slope/{tc}_{run}_{ch}_Slope.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.lowest_rate_run, ch=self.channel)
        angle = self.CutConfig['track_angle']

        def func():
            print 'generating slope cut for ch{ch}...'.format(ch=self.channel)
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

    def generate_bucket(self):
        num = self.analysis.get_signal_numbers('e', 2)
        name = self.analysis.get_signal_names(num)[self.channel]
        string = '{sig2}-{sig1}==0'.format(sig2=name, sig1=self.analysis.signal_names[self.channel])
        return TCut(string)

    def generate_cut_string(self, set_channel=True):
        """
        Creates the cut string, which will be stored in self.cut. With the arguments set to False, different cut types can be deactivated in the cut string.
        :param set_channel:
        :return:
        """
        gROOT.SetBatch(1)
        if self._checklist["GenerateCutString"]:
            self.load_config()  # re-generate
        cutstring = self.cut

        # --CHI2 --
        self.cut_strings['chi2'] = self.generate_chi2()
        self.cut_strings['track_angle'] = self.generate_slope()
        self.cut_strings['bucket'] = self.generate_bucket()

        # -- EVENT RANGE CUT --
        self.generate_event_range()
        if self.CutConfig["EventRange"]:
            if cutstring != "":
                cutstring += "&&"
            cutstring += "(event_number<={maxevent}&&event_number>={minevent})".format(minevent=self.CutConfig["EventRange"][0], maxevent=self.CutConfig["EventRange"][1])
            self.EasyCutStrings["EventRange"] = "Evts.{min}k-{max}k".format(min=int(self.CutConfig["EventRange"][0]) / 1000, max=int(self.CutConfig["EventRange"][1]) / 1000)
            self.EasyCutStrings["ExcludeFirst"] = ""
        elif self.CutConfig["ExcludeFirst"] > 0:
            if cutstring != "":
                cutstring += "&&"
            cutstring += "event_number>={minevent}".format(minevent=self.CutConfig["ExcludeFirst"])
            self.EasyCutStrings["ExcludeFirst"] = "Evts.{min}k+".format(min=int(self.CutConfig["ExcludeFirst"]) / 1000)
            self.EasyCutStrings["EventRange"] = ""
        else:
            self.EasyCutStrings["EventRange"] = ""
            self.EasyCutStrings["ExcludeFirst"] = ""

        # -- PULSER CUT --
        self.cut_strings['pulser'] += '!pulser'

        # -- SATURATED CUT --
        self.cut_strings['saturated'] = '!is_saturated[{ch}]'.format(ch=self.channel)

        # -- TRACK CUT --
        self.cut_strings['tracks'] += 'n_tracks'

        # -- SPREAD LOW CUT --
        if self.CutConfig['spread_low'] > 0:
            self.cut_strings['spread_low'] += 'sig_spread[{ch}]>{low}'.format(ch=self.channel, low=int(self.CutConfig['spread_low']))
            if cutstring != '':
                cutstring += '&&'
            cutstring += 'sig_spread[{channel}]>{low}'.format(channel=self.channel, low=int(self.CutConfig['spread_low']))
            self.EasyCutStrings['spread_low'] = 'spread>{low}'.format(low=int(self.CutConfig['spread_low']))
        else:
            self.EasyCutStrings['spread_low'] = ''

        # -- MEDIAN CUT --
        if self.CutConfig['absMedian_high'] > 0:
            print 'CONFIG:', self.CutConfig['absMedian_high'], type(self.CutConfig['absMedian_high'])
            self.cut_strings['median'] += 'abs(median[{ch}])<{high}'.format(ch=self.channel, high=int(self.CutConfig['absMedian_high']))
            if cutstring != '':
                cutstring += '&&'
            cutstring += 'abs(median[{channel}])<{high}'.format(channel=self.channel, high=int(self.CutConfig['absMedian_high']))
            self.EasyCutStrings['absMedian_high'] = '|median|<{high}'.format(high=int(self.CutConfig['absMedian_high']))
        else:
            self.EasyCutStrings['absMedian_high'] = ''

        # -- PEDESTAL SIGMA CUT --
        if self.CutConfig['pedestalsigma'] > 0:
            if cutstring != '':
                cutstring += '&&'
            self.__load_pedestal_data()
            ped_range = self.__calc_pedestal_range()
            string = '{ped}>{min}&&{ped}<{max}'.format(ped=self.analysis.pedestal_names[self.channel], min=ped_range[0], max=ped_range[1])
            cutstring += string
            self.EasyCutStrings["pedestalsigma"] = "PedSigma" + str(self.CutConfig['pedestalsigma'])
            self.cut_strings['ped_sigma'] += string
        else:
            self.EasyCutStrings["pedestalsigma"] = ""

        # -- set the channel on the cuts --
        if set_channel:
            self.cut = cutstring
            self.cut = self.cut.format(channel=self.channel)

        # -- BEAM INTERRUPTION CUT --
        self.__generate_beam_interruptions()
        self.EasyCutStrings["noBeamInter"] = "BeamOn"

        self._checklist["GenerateCutString"] = True
        gROOT.SetBatch(0)

    def __calc_pedestal_range(self):
        sigma = self.pedestal_sigma
        ped_mean = self.pedestal_mean
        sigma_range = self.CutConfig['pedestalsigma']
        return [ped_mean - sigma_range * sigma, ped_mean + sigma_range * sigma]

    def __generate_beam_interruptions(self, ):
        """
        This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied.
        :return: cut
        """
        self.get_beam_interruptions()

        njumps = len(self.jump_ranges["start"])
        cut_string = ''
        for i in xrange(njumps):
            string = "!(event_number<={upper}&&event_number>={lower})".format(upper=self.jump_ranges["stop"][i], lower=self.jump_ranges["start"][i])
            if self.cut != "":
                self.cut += "&&"
            self.cut += string
            # new separate strings
            if cut_string != '':
                cut_string += '&&'
            cut_string += string
        self.cut_strings['beam_interruptions'] += cut_string
        self._checklist["RemoveBeamInterruptions"] = True

        return self.cut
    # endregion

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

    def get_function_def(self):
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
        # defstring_ = "lambda pulser, is_saturated, n_tracks, fft_mean, INVfft_max, sig_time, sig_spread, median: "
        # def_ = ""
        #
        # if self.cut_types["IndividualChCut"] != "":
        #     raw_input("WARNING: cut0 and cut3 cannot be applied on Tracking Data! Press Enter to Continue..")
        #
        # if self.cut_types["noPulser"] == 1:
        #     def_ += "not pulser"
        # elif self.cut_types["noPulser"] == 0:
        #     def_ += "pulser"
        #
        # if self.cut_types["notSaturated"]:
        #     if def_ != "":
        #         def_ += " and "
        #     def_ += "not is_saturated"
        #
        # if self.cut_types["Tracks"]:
        #     if def_ != "":
        #         def_ += " and "
        #     def_ += "n_tracks"
        #
        # if self.cut_types["FFT"]:
        #     if def_ != "":
        #         def_ += " and "
        #     assert False, "FFT cut not yet implemented in GetCutFunctionDef() method of Cut class. "
        #     # to do: FFT entry in _cutTypes should be dict and/or contain a TCutG instance
        #
        # if self.cut_types["peakPos_high"] > 0:
        #     if def_ != "":
        #         def_ += " and "
        #     def_ += "sig_time<{high}".format(high=int(self.cut_types["peakPos_high"]))
        #
        # if self.cut_types["spread_low"] > 0:
        #     if def_ != "":
        #         def_ += " and "
        #     def_ += "sig_spread>{low}".format(low=int(self.cut_types["spread_low"]))
        #
        # if self.cut_types["absMedian_high"] > 0:
        #     if def_ != "":
        #         def_ += " and "
        #     def_ += "abs(median)>{low}".format(low=int(self.cut_types["absMedian_high"]))

        return self.all_cut

    def get_easy_cutstring(self):
        """
        Returns a short, more user-friendly cut string, which can be
        used to display the cut configuration as terminal prompt or
        inside a canvas.
        :return:
        """
        string_ = ""
        for type_ in self.EasyCutStrings.keys():
            if self.EasyCutStrings[type_] != "":
                string_ += self.EasyCutStrings[type_] + ", "
        if string_ != "":
            string_ = string_[:-2]
        return string_

    def show_cuts(self):
        for key, value in self.cut_strings.iteritems():
            print key, value
        return
