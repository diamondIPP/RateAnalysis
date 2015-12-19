from RunClass import Run
import types as t
from collections import namedtuple
import json
import copy
# from ConfigParser import ConfigParser
from datetime import datetime as dt


class RunSelection(Run):
    def __init__(self, verbose=False):
        Run.__init__(self, run_number=None, verbose=verbose)
        
        self.runplan_path = self.get_program_dir() + self.run_config_parser.get('BASIC', 'runplaninfofile')
        self.runplan = self.load_runplan()
        self.run_numbers = self.load_run_numbers()
        self.logs = {}
        self._LoadRuns()
        self._InitializeSelections()
        self.run_plan = None

    def __str__(self):
        nr = len(self.run_numbers)
        selected_runs = self.GetSelectedRuns()
        return "RunSelection Object\n" + str(len(selected_runs)) + " Out of " + str(nr) + " runs selected. Selections made:" + self.get_log_string()

    def make_log_entry(self, event):
        time_str = dt.now().strftime('%H:%M:%S')
        self.logs[len(self.logs)] = [event, time_str]

    def get_log_string(self):
        string = '\n'
        for key, log in self.logs.iteritems():
            string += '{key}.)\t{log}'.format(key=key, log=log[0])
        return string

    def print_logs(self):
        for key, log in self.logs.iteritems():
            print '{key}.)\t{time}\t{log}'.format(key=key, time=log[1], log=log[0])

    def load_run_numbers(self):
        run_numbers = []
        f = open(self.runinfofile, 'r')
        data = json.load(f)
        f.close()
        # parser = ConfigParser()
        # parser.read('Configuration/KeyDict_{campaign}.cfg'.format(campaign=self.TESTCAMPAIGN))
        # type_name = parser.get('KEYNAMES', 'type')
        # signal_types = ['signal', 'rate_scan']
        for key in data:
            run_numbers.append(int(key))
        return sorted(run_numbers)

    def _LoadRuns(self):
        """
        loads all the run infos in a dict with the run numbers as keys
        :return:
        """
        self.run_infos = {}
        for runnumber in self.run_numbers:
            self.set_run(runnumber, load_root_file=False)
            self.run_infos[runnumber] = self.current_run

    def _InitializeSelections(self):
        '''
        creates a dict of bools to store selections. dict is filled with False (no run selected)
        the selections log is created/cleared
        :return:
        '''
        self.selections = {}
        self.channel_selections = {}
        self.logs = {}
        for runnumber in self.run_numbers:
            self.selections[runnumber] = False
            self.channel_selections[runnumber] = {
                0: False,
                3: False
            }

    def load_runplan(self):
        f = open(self.runplan_path, 'r')
        runplans = json.load(f)
        f.close()
        try:
            runplan = runplans[self.TESTCAMPAIGN]
        except KeyError:
            print 'No runplan for {tc} available yet, creating an empty one!'.format(tc=self.TESTCAMPAIGN)
            runplan = {}
            self.save_runplan(runplan)
        return runplan

    def save_runplan(self, runplan=None):
        f = open(self.runplan_path, 'r+')
        runplans = json.load(f)
        runplans[self.TESTCAMPAIGN] = self.run_plan if runplan is None else runplan
        f.seek(0)
        json.dump(runplans, f, indent=2, sort_keys=True)
        f.close()

    def SelectAll(self, selectDiamond1=True, selectDiamond2=True):
        '''
        Selects all runs.
        :param selectDiamond1:
        :param selectDiamond2:
        :return:
        '''
        for run_number in self.run_numbers:
            self.selections[run_number] = True
            if selectDiamond1: self.channel_selections[run_number][0] = True
            if selectDiamond2: self.channel_selections[run_number][0] = True
        self.make_log_entry('All runs selected')
        self.VerbosePrint('All runs selected')

    def UnSelectAll(self):
        '''
        Resets all selections made.
        :return:
        '''
        self._InitializeSelections()
        self.VerbosePrint('All runs unselected')

    def UnselectAll(self):
        '''
        Resets all selections made.
        :return:
        '''
        self.UnSelectAll()

    def set_channelss(self, diamond1=True, diamond2=True):
        '''
        Sets the channels (diamonds) of the selected runs to active or
        inactive.
        :param diamond1:
        :param diamond2:
        :return:
        '''
        for run_number in self.GetSelectedRuns():
            self.channel_selections[run_number][0] = diamond1
            self.channel_selections[run_number][3] = diamond2
        change1 = "diamond 1 = " + str(diamond1)
        change2 = "diamond 2 = " + str(diamond2)
        self.make_log_entry("Channels of selected runs set: " + change1 + " " + change2)

    def SelectDataType(self, data_type):  # data type
        '''
        Selects the runs according to the type of run, such as rate_scan,
        test, voltage_scan etc..
        Selects just the runs, NOT the channels (diamonds)!
        To select the channels use SelectDiamondRuns()
        :param data_type:
        :return:
        '''
        types = self._GetValues("type")
        assert (data_type in types), "wrong data type. \n\tSelect type from: " + str(types)
        count = 0
        for run_number in self.run_numbers:
            if self.run_infos[run_number]['type'] == data_type:
                self._SelectRun(run_number)
                count += 1
        self.make_log_entry('Runs of Type ' + data_type + ' selected. +' + str(count) + ' selections')
        self.VerbosePrint('Runs of Type ' + data_type + ' selected. +' + str(count) + ' selections')

    def SelectDiamondRuns(self, diamondname, only_selected_runs=False):
        '''
        Selects all runs, which have the diamond with name 'diamondname'
        in it. It Furthermore selects also the channels corresponding
        to this diamondname.
        :param diamondname:
        :param only_selected_runs:
        :return:
        '''
        diamondnames = self.ShowDiamondNames(getNames=True)
        assert (diamondname in diamondnames), "wrong diamond name. \n\tSelect diamond name from: " + str(diamondnames)
        if not only_selected_runs:
            choice = self.run_numbers
        else:
            choice = self.GetSelectedRuns()
        count = 0
        for run_number in choice:
            if self.run_infos[run_number]['diamond 1'] == diamondname:
                self.selections[run_number] = True
                self.channel_selections[run_number][0] = True
                count += 1
            if self.run_infos[run_number]['diamond 2'] == diamondname:
                self.selections[run_number] = True
                self.channel_selections[run_number][3] = True
                count += 1
        self.make_log_entry('Runs and Channels containing ' + diamondname + ' selected. +' + str(count) + ' runs selected')
        self.VerbosePrint('Runs and Channels containing ' + diamondname + ' selected. +' + str(count) + ' runs selected')

    # def SelectIrradiationRuns(self, irradiated=True, irrtype=None):
    #     count = 0
    #     if irradiated and irrtype == None:
    #         for run_number in self.run_numbers:
    #             self.SetRun(run_number, validate=False)
    #             if self.diamond.Specifications['Irradiation'] == 'proton' or self.diamond.Specifications['Irradiation'] == 'neutron':
    #                 self.selections[run_number] = True
    #                 count += 1
    #         self._Log('Irradiated Runs selected')
    #     elif not irradiated:
    #         for run_number in self.run_numbers:
    #             self.SetRun(run_number, validate=False)
    #             if self.diamond.Specifications['Irradiation'] == 'no':
    #                 self.selections[run_number] = True
    #                 count += 1
    #         self._Log('Non-Irradiated Runs selected')
    #     else:
    #         assert(irrtype in ['no', 'proton', 'neutron']), "wrong irradiation type. Choose irrtype in [`proton`, `neutron`, `no`]"
    #         if irrtype == 'no':
    #             self.SelectIrradiationRuns(irradiated=False)
    #         else:
    #             for run_number in self.run_numbers:
    #                 self.SetRun(run_number, validate=False)
    #                 if self.diamond.Specifications['Irradiation'] == irrtype:
    #                     self.selections[run_number] = True
    #                     count += 1
    #             self._Log('Irradiated Runs selected with '+irrtype+' irradiation. +'+str(count)+' selections')
    #             self.VerbosePrint('Irradiated Runs selected with '+irrtype+' irradiation. +'+str(count)+' selections')

    def UnSelectUnlessDataType(self, data_type):
        '''
        Keeps only runs which are of type 'data_type'.
        :param data_type:
        :return:
        '''
        types = self._GetValues("type")
        assert (data_type in types), "wrong data type. \n\tSelect type from: " + str(types)
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if self.run_infos[run_number]['type'] == data_type:
                    pass
                else:
                    self._UnselectRun(run_number)
                    count += 1
        self.make_log_entry('All Selected Runs unselected if not of Type ' + data_type + '. -' + str(count) + ' selections')
        self.VerbosePrint('All Selected Runs unselected if not of Type ' + data_type + '. -' + str(count) + ' selections')

    # def UnSelectUnlessIrradiation(self, irradiated=True, irrtype=None):
    #     count = 0
    #     if irradiated and irrtype == None:
    #         for run_number in self.run_numbers:
    #             if self.selections[run_number]:
    #                 self.SetRun(run_number, validate=False)
    #                 if self.diamond.Specifications['Irradiation'] == 'proton' or self.diamond.Specifications['Irradiation'] == 'neutron':
    #                     pass
    #                 else:
    #                     self.selections[run_number] = False
    #                     count += 1
    #         self._Log('All Selected Runs unselected if non-irradiated. Only radiated Runs left. -'+str(count)+' selections')
    #         self.VerbosePrint('All Selected Runs unselected if non-irradiated. Only radiated Runs left. -'+str(count)+' selections')
    #     elif not irradiated:
    #         for run_number in self.run_numbers:
    #             if self.selections[run_number]:
    #                 self.SetRun(run_number, validate=False)
    #                 if self.diamond.Specifications['Irradiation'] == 'no':
    #                     pass
    #                 else:
    #                     self.selections[run_number] = False
    #                     count += 1
    #         self._Log('All Selected Runs unselected if irradiated. Only non-radiated Runs left. -'+str(count)+' selections')
    #         self.VerbosePrint('All Selected Runs unselected if irradiated. Only non-radiated Runs left. -'+str(count)+' selections')
    #     else:
    #         assert(irrtype in ['no', 'proton', 'neutron']), "wrong irradiation type. Choose irrtype in [`proton`, `neutron`, `no`]"
    #         if irrtype == 'no':
    #             self.UnSelectUnlessIrradiation(irradiated=False)
    #         else:
    #             for run_number in self.run_numbers:
    #                 if self.selections[run_number]:
    #                     self.SetRun(run_number, validate=False)
    #                     if self.diamond.Specifications['Irradiation'] == irrtype:
    #                         pass
    #                     else:
    #                         self.selections[run_number] = False
    #                         count += 1
    #             self._Log('All Selected Runs unselected if not radiated by '+irrtype+'. Only '+irrtype+'-radiated Runs left. -'+str(count)+' selections')
    #             self.VerbosePrint('All Selected Runs unselected if not radiated by '+irrtype+'. Only '+irrtype+'-radiated Runs left. -'+str(count)+' selections')

    def UnSelectUnlessDiamond(self, diamondname):
        '''
        Keeps only runs which hold the diamond with name 'diamondname'.
        :param diamondname:
        :return:
        '''
        diamondnames = self.ShowDiamondNames(True)
        assert (diamondname in diamondnames), "wrong diamond name. \n\tSelect diamond name from: " + str(diamondnames)
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if diamondname in [self.run_infos[run_number]['diamond 1'], self.run_infos[run_number]['diamond 2']]:
                    pass
                else:
                    self._UnselectRun(run_number)
                    count += 1
        self.make_log_entry('All Selected Runs unselected if not using ' + diamondname + ' diamond. Only runs countaining ' + diamondname + ' left. -' + str(count) + ' selections')
        self.VerbosePrint('All Selected Runs unselected if not using ' + diamondname + ' diamond. Only runs countaining ' + diamondname + ' left. -' + str(count) + ' selections')

    def UnSelectUnlessBias(self, bias):
        '''
        Keeps only runs selected which have a diamond with a given bias
        voltage. Diamonds with a different bias voltage will be un-
        selected.
        :param bias:
        :return:
        '''
        assert (type(bias) == t.IntType), "Bias has to be int-Type"
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                unselectrun = True
                if self.run_infos[run_number]['hv dia1'] == bias:
                    unselectrun = False
                else:
                    self.channel_selections[run_number][0] = False
                if self.run_infos[run_number]['hv dia2'] == bias:
                    unselectrun = False
                else:
                    self.channel_selections[run_number][3] = False
                if unselectrun:
                    self.selections[run_number] = False
                    count += 1
        self.make_log_entry('All Selected Runs unselected if not ' + str(bias) + 'V bias applied. Only ' + str(bias) + 'V Bias Runs left. -' + str(count) + ' selections')
        self.VerbosePrint('All Selected Runs unselected if not ' + str(bias) + 'V bias applied. Only ' + str(bias) + 'V Bias Runs left. -' + str(count) + ' selections')

    def _SelectRun(self, run_number):
        assert (run_number in self.run_numbers), "run number " + str(run_number) + " not found in list of run numbers. Check run_log json file!"
        self.selections[run_number] = True

    def _UnselectRun(self, run_number):
        self.selections[run_number] = False
        self.channel_selections[run_number][0] = False
        self.channel_selections[run_number][3] = False

    def ExcludeRuns(self, run_number):
        '''
        This method will un-select the run with number run_number.
        run_number has to be a single integer or a list of integers in
        order to exclude several runs.
        :param run_number:
        :return:
        '''
        assert (type(run_number) == t.IntType or type(run_number) == t.ListType), "Wrong input type. run_number has to be either integer or list of integer"
        listOfRuns = self.GetSelectedRuns()
        if type(run_number) == t.IntType:
            if run_number in listOfRuns:
                self._UnselectRun(run_number)
                self.make_log_entry("Run " + str(run_number) + " unselected. -1 selection")
        else:
            ListToExclude = run_number
            for run_number in ListToExclude:
                self.ExcludeRuns(int(run_number))

    def SelectRunsInRange(self, minrun, maxrun):
        '''
        Selects all runs with run numbers in range [minrun, .. , maxrun].
        :param minrun:
        :param maxrun:
        :return:
        '''
        for run_number in self.run_numbers:
            if run_number <= maxrun and run_number >= minrun:
                self._SelectRun(run_number)

    def UnSelectUnlessInRange(self, minrun, maxrun):
        '''
        Keeps only runs in range [minrun, .. , maxrun].
        :param minrun:
        :param maxrun:
        :return:
        '''
        for run_number in self.GetSelectedRuns():
            if not (run_number <= maxrun and run_number >= minrun):
                self._UnselectRun(run_number)

    # def ValidateSelectedRuns(self):
    #     #self.ValidateRuns(self.GetSelectedRuns())
    #     pass

    def ShowSelectedRuns(self, show_allcomments=False, commentlength=15):
        '''
        Lists all selected runs.
        :param show_allcomments:
        :param commentlength:
        :return:
        '''
        print len(self.GetSelectedRuns()), " Runs Selected:"

        def multilinetext(text, width):
            length = len(text)
            word_wraps = length / int(width)
            separator = "\n "
            for i in xrange(word_wraps):
                text = text[:((i + 1) * width + len(separator) * i)] + separator + text[(i + 1) * width + len(separator) * i:]
            return text

        if show_allcomments:
            printstring = "{nr} {type} {dia1sel}{dia1} {hv1:>5} {dia2sel}{dia2} {hv2:>5} {rate}"
        else:
            printstring = "{nr} {type} {dia1sel}{dia1} {hv1:>5} {dia2sel}{dia2} {hv2:>5} {rate} {comment}"

        Record = namedtuple("Record", ["runnumber", "type", "dia1", "bias1", "dia2", "bias2", "rate", "comment"])
        listitems = []
        for runnumber in self.GetSelectedRuns():
            listitems += [Record(
                runnumber,
                self.run_infos[runnumber]['type'],
                self.run_infos[runnumber]['diamond 1'],
                self.run_infos[runnumber]['hv dia1'],
                self.run_infos[runnumber]['diamond 2'],
                self.run_infos[runnumber]['hv dia2'],
                int(self.run_infos[runnumber]['measured flux']),
                [self.run_infos[runnumber]['user comments']]
            )]

        print printstring.format(
            nr="Nr.".ljust(3),
            type="Type".ljust(9),
            dia1="Dia 1".ljust(7),
            dia1sel=" ",
            hv1="HV 1".ljust(5),
            dia2="Dia 2".ljust(7),
            dia2sel=" ",
            hv2="HV2 ".ljust(5),
            rate="Rate".ljust(6),
            comment="Comment"
        )
        for item in listitems:
            if len(item.comment[0][:30]) > 0:
                marker = "* "
            else:
                marker = ""
            dia1sel = " "
            dia2sel = " "
            if self.channel_selections[item.runnumber][0]: dia1sel = "*"
            if self.channel_selections[item.runnumber][3]: dia2sel = "*"
            print printstring.format(
                nr=str(item.runnumber).ljust(3),
                type=item.type.ljust(9),
                dia1=item.dia1.ljust(7),
                dia1sel=dia1sel,
                hv1=int(item.bias1),
                dia2=item.dia2.ljust(7),
                dia2sel=dia2sel,
                hv2=int(item.bias2),
                rate=(str(item.rate) + " kHz").rjust(8),
                comment=(marker + item.comment[0][:commentlength]).ljust(commentlength + 2)
            )
            if show_allcomments and len(item.comment[0][:30]) > 0:
                print "   -- COMMENT: -----------------------------------"
                print multilinetext(" " + item.comment[0][:], 50)
                print "--------------------------------------------------"

    def show_run_info(self, runs=None, detailed=False):
        """
        Prints all run infos from the selected runs to the console.
        :param runs:
        :return:
        """

        if detailed:
            if runs is None:
                for runnumber in self.GetSelectedRuns():
                    self._printRunInfo(runnumber)
            else:
                if type(runs) is list:
                    for runnumber in runs:
                        self._printRunInfo(runnumber)
                elif type(runs) is int:
                    self._printRunInfo(runs)
                else:
                    print "Wrong input type"
        else:
            if runs is None:
                self.__print_runinfo_header()
                for runnumber in self.GetSelectedRuns():
                    self.__print_runinfo(runnumber)
            else:
                # todo:
                print 'not yet implemented'

    def _printRunInfo(self, runnumber):
        print "--- RUN ", runnumber, " ---"
        for key in self.run_infos[runnumber].keys():
            print "{key:>20}: {value}".format(key=key, value=self.run_infos[runnumber][key])

    def __print_runinfo(self, run):
        dia1 = self.run_infos[run]['diamond 1']
        dia2 = self.run_infos[run]['diamond 2']
        flux = self.run_infos[run]['measured flux']
        type = self.run_infos[run]['type']
        print '{run}\t{type}\t{dia1}\t{dia2}\t{flux}'.format(run=run, dia1=dia1, dia2=dia2, flux=flux, type=type)

    @staticmethod
    def __print_runinfo_header():
        print 'run\ttype\tdia1\tdia2\tflux'

    def ShowDiamondNames(self, getNames=False):
        '''
        Prints all diamond names from log file into the console.
        If getNames is True, it will return the diamond names as a list
        instead of printing it.
        :param getNames:
        :return:
        '''
        diamondnames = []
        diamond1names = self._GetValues("diamond 1")
        diamond2names = self._GetValues("diamond 2")
        biglist = diamond1names + diamond2names
        for name in biglist:
            if not name in diamondnames:
                diamondnames += [name]
        if not getNames:
            print "Diamondnames:"
            for i in xrange(len(diamondnames)):
                print "\t" + diamondnames[i]
        else:
            return diamondnames

    def ShowTypes(self):
        '''
        Prints all run types from log file into the console.
        :return:
        '''
        types = self._GetValues("type")
        print "Types:"
        for i in types:
            print "\t" + i

    def _GetValues(self, key):
        valuelist = []
        for runnumber in self.run_infos.keys():
            value = self.run_infos[runnumber][key]
            if not value in valuelist:
                valuelist += [value]
        return valuelist

    def GetSelectedRuns(self):
        '''
        Returns the list of selected run numbers.
        :return:
        '''
        selected = []
        for runnumber in self.run_numbers:
            if self.selections[runnumber]:
                selected.append(runnumber)
        if self.verbose and len(selected) == 0:
            print "No Runs Selected"
        return selected

    def GetSelectedDiamonds(self):
        '''
        Returns a list, containing for each selected run an integer
        according to the diamond selection configuration. (i.e. which
        diamonds are selected for analysis).
            1 -> Diamond 1
            2 -> Diamond 2
            3 -> Diamond 1 & 2, or no diamond selection (default: both)
        :return:
        '''
        selected = []
        for runnumber in self.run_numbers:
            if self.selections[runnumber]:
                dia1 = self.channel_selections[runnumber][0]
                dia2 = self.channel_selections[runnumber][3]
                diamonds = int(dia1) * (1 << 0) + int(dia2) * (1 << 1)
                if diamonds == 0: diamonds = 3
                selected.append(diamonds)
        if self.verbose and len(selected) == 0:
            print "No Runs Selected"
        return selected

    def ShowRunPlan(self, detailed=True, type_="rate_scan", show_allcomments=False, commentlength=0):
        '''
        Print a list of all run plans from the current test campaign
        to the console.
        The run plans are defined and saved via
            AddSelectedRunsToRunPlan()
        :param detailed:
        :param type_:
        :param show_allcomments:
        :param commentlength:
        :return:
        '''
        print "RUN PLAN FOR TESTCAMPAIGN:", self.TESTCAMPAIGN
        if not detailed:
            for types_ in self.runplan.keys():
                print types_, " :"
                numbers = map(int, self.runplan[types_].keys())
                numbers.sort()
                for planNr in numbers:
                    print "\t Nr. {nr} : {runs}".format(nr=planNr, runs=self.runplan[types_][str(planNr)])
        else:
            tmp_selections = copy.deepcopy(self.selections)
            tmp_channel_selections = copy.deepcopy(self.channel_selections)
            tmp_selectionLog = copy.deepcopy(self.logs)

            numbers = map(int, self.runplan[type_].keys())
            numbers.sort()
            for i in numbers:
                self.UnselectAll()
                self.SelectRunsFromRunPlan(i, type_=type_)
                print "-----------------------------------------"
                print "RUN PLAN ", i, " ({typ}) : ".format(typ=type_)
                print "-----------------------------------------"
                self.ShowSelectedRuns(show_allcomments=show_allcomments, commentlength=commentlength)
                print "\n"
            self.selections = copy.deepcopy(tmp_selections)
            self.channel_selections = copy.deepcopy(tmp_channel_selections)
            self.logs = copy.deepcopy(tmp_selectionLog)

    def SelectRunsFromRunPlan(self, number, type_="rate_scan"):
        '''
        Selects all runs corresponding to the run plan with key 'number'.
        :param number:
        :param type_:
        :return:
        '''
        runs = self.runplan[type_][str(number)]
        self.run_plan = number
        self.SelectRuns(list_of_runs=runs)

    def AddSelectedRunsToRunPlan(self, key, run_type="rate_scan"):
        '''
        Saves all selected runs as a run plan with name 'key'. This
        run plan can later be imported via
            .SelectRunsFromRunPlan('key')
        :param key:
        :param run_type:
        :return:
        '''
        assert (run_type in ["rate_scan", "voltage_scan", "test"])
        if not type(key) is t.StringType:
            key = str(key)
        if self.runplans[self.TESTCAMPAIGN].has_key(run_type):
            self.runplans[self.TESTCAMPAIGN][run_type][key] = self.GetSelectedRuns()
        else:
            self.runplans[self.TESTCAMPAIGN][run_type] = {}
            self.runplans[self.TESTCAMPAIGN][run_type][key] = self.GetSelectedRuns()
        self.save_runplan()

    def SelectRuns(self, list_of_runs, select_dia1=False, select_dia2=False):
        '''
        Selects all runs with run number in list_of_runs. The diamonds
        can be selected individually.
        :param list_of_runs:
        :param select_dia1:
        :param select_dia2:
        :return:
        '''
        assert (type(list_of_runs) is t.ListType), "list_of_runs not a list"
        for runnumber in list_of_runs:
            self._SelectRun(runnumber)

if __name__ == '__main__':
    z = RunSelection()