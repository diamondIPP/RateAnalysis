from RunClass import Run
from Elementary import Elementary
import types as t
from collections import namedtuple
import json
import copy
# from ConfigParser import ConfigParser
from datetime import datetime as dt


class RunSelection(Elementary):
    def __init__(self, verbose=False):
        Elementary.__init__(self, verbose)
        self.run = Run(run_number=None, verbose=verbose)

        self.runplan_path = self.get_program_dir() + self.run.run_config_parser.get('BASIC', 'runplaninfofile')
        self.run_plan = self.load_runplan()
        self.run_numbers = self.load_run_numbers()
        self.run_infos = self.load_run_infos()
        self.logs = {}
        self.selection = {}
        self.channels = {}

        self.init_selection()

    def __str__(self):
        nr = len(self.run_numbers)
        selected_runs = self.get_selected_runs()
        return 'RunSelection Object\n' + str(len(selected_runs)) + ' Out of ' + str(nr) + ' runs selected. Selections made:' + self.get_log_string()

    # ============================================
    # region LOGGING
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
    # endregion

    # ============================================
    # region INIT
    def load_run_numbers(self):
        run_numbers = []
        f = open(self.run.runinfofile, 'r')
        data = json.load(f)
        f.close()
        # parser = ConfigParser()
        # parser.read('Configuration/KeyDict_{campaign}.cfg'.format(campaign=self.TESTCAMPAIGN))
        # type_name = parser.get('KEYNAMES', 'type')
        # signal_types = ['signal', 'rate_scan']
        for key in data:
            run_numbers.append(int(key))
        return sorted(run_numbers)

    def load_run_infos(self):
        """ loads all the run infos in a dict with the run numbers as keys """
        run_infos = {}
        for runnumber in self.run_numbers:
            self.run.set_run(runnumber, load_root_file=False)
            run_infos[runnumber] = self.run.RunInfo
        return run_infos

    def init_selection(self):
        self.reset_selection()

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
    # endregion

    def reset_selection(self):
        """ Creates a dict of bools to store the selection, which is filled with False (no run selected). Resets the logs. """
        self.logs = {}
        for run in self.run_numbers:
            self.selection[run] = False
            self.channels[run] = {}
            for ch in self.run.channels:
                self.channels[run][ch] = False

    def save_runplan(self, runplan=None):
        f = open(self.runplan_path, 'r+')
        runplans = json.load(f)
        runplans[self.TESTCAMPAIGN] = self.run_plan if runplan is None else runplan
        f.seek(0)
        json.dump(runplans, f, indent=2, sort_keys=True)
        f.close()

    def select_all_runs(self, dia1=True, dia2=True):
        for run in self.run_numbers:
            self.selection[run] = True
            self.channels[run][self.run.channels[0]] = dia1
            self.channels[run][self.run.channels[1]] = dia2
        self.make_log_entry('All runs selected')
        self.verbose_print('All runs selected')

    def unselect_all_runs(self):
        self.reset_selection()
        self.verbose_print('All runs unselected')

    def set_channels(self, dia1=True, dia2=True):
        """
        Sets the channels (diamonds) of the selected runs to active or inactive.
        :param dia1:
        :param dia2:
        """
        dias = [dia1, dia2]
        for run_number in self.get_selected_runs():
            for i, ch in enumerate(self.run.channels):
                self.channels[run_number][ch] = dias[i]
        self.make_log_entry('Channels of selected runs set: diamond1 to {dia1}, diamond2 to {dia2}'.format(dia1=dia1, dia2=dia2))

    def reset_channels(self, run):
        for ch in self.run.channels:
            self.channels[run][ch] = False

    def select_runs_of_type(self, run_type, unselect=False, only_selected=False):
        """
        Selects the runs according to the type of run, such as rate_scan, test, voltage_scan etc..
        :param run_type:
        :param unselect:
        :param only_selected:
        """
        types = self.get_runinfo_values('type')
        assert run_type in types, 'wrong data type.\n\t-->Select type from: {types}'.format(types=types)
        runs = self.get_selected_runs() if only_selected else self.run_numbers
        selected_runs = 0
        for run in runs:
            if self.run_infos[run]['type'] == run_type:
                self.select_run(run, False) if not unselect else self.unselect_run(run, False)
                selected_runs += 1
        prefix = 'un' if unselect else ''
        self.make_log_entry('Runs of type {type} {pref}selected ({nr} {pref}selections).'.format(type=run_type, pref=prefix, nr=selected_runs))
        self.verbose_print('Runs of type {type} {pref}selected ({nr} {pref}selections).'.format(type=run_type, pref=prefix, nr=selected_runs))

    def unselect_runs_of_type(self, run_type):
        self.select_runs_of_type(run_type, unselect=True)

    def select_diamond_runs(self, diamondname, only_selected_runs=False):
        """
        Selects all runs, which have the diamond with name 'diamondname' in it. It Furthermore selects also the channels corresponding to this diamondname.
        :param diamondname:
        :param only_selected_runs:
        """
        diamondnames = self.get_diamond_names()
        assert diamondname in diamondnames, 'wrong diamond name. \n\t-->Select diamond name from: {dias}'.format(dias=diamondnames)
        runs = self.get_selected_runs() if only_selected_runs else self.run_numbers
        selected_runs = 0
        dia_keys = ['diamond 1', 'diamond 2']
        for i, run in enumerate(runs):
            if self.run_infos[run][dia_keys[i]] == diamondname:
                self.selection[run] = True
                self.channels[run][0] = True
                selected_runs += 1
            if self.run_infos[run]['diamond 2'] == diamondname:
                self.selection[run] = True
                self.channels[run][3] = True
                selected_runs += 1
        self.make_log_entry('Runs and Channels containing ' + diamondname + ' selected. +' + str(selected_runs) + ' runs selected')
        self.verbose_print('Runs and Channels containing ' + diamondname + ' selected. +' + str(selected_runs) + ' runs selected')

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
    #         assert(irrtype in ['no', 'proton', 'neutron']), 'wrong irradiation type. Choose irrtype in [`proton`, `neutron`, `no`]'
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
        types = self.get_runinfo_values('type')
        assert (data_type in types), 'wrong data type. \n\tSelect type from: ' + str(types)
        count = 0
        for run_number in self.run_numbers:
            if self.selection[run_number]:
                if self.run_infos[run_number]['type'] == data_type:
                    pass
                else:
                    self.unselect_run(run_number)
                    count += 1
        self.make_log_entry('All Selected Runs unselected if not of Type ' + data_type + '. -' + str(count) + ' selections')
        self.verbose_print('All Selected Runs unselected if not of Type ' + data_type + '. -' + str(count) + ' selections')

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
    #         assert(irrtype in ['no', 'proton', 'neutron']), 'wrong irradiation type. Choose irrtype in [`proton`, `neutron`, `no`]'
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
        diamondnames = self.show_diamond_names(True)
        assert (diamondname in diamondnames), 'wrong diamond name. \n\tSelect diamond name from: ' + str(diamondnames)
        count = 0
        for run_number in self.run_numbers:
            if self.selection[run_number]:
                if diamondname in [self.run_infos[run_number]['diamond 1'], self.run_infos[run_number]['diamond 2']]:
                    pass
                else:
                    self.unselect_run(run_number)
                    count += 1
        self.make_log_entry('All Selected Runs unselected if not using ' + diamondname + ' diamond. Only runs countaining ' + diamondname + ' left. -' + str(count) + ' selections')
        self.verbose_print('All Selected Runs unselected if not using ' + diamondname + ' diamond. Only runs countaining ' + diamondname + ' left. -' + str(count) + ' selections')

    def UnSelectUnlessBias(self, bias):
        '''
        Keeps only runs selected which have a diamond with a given bias
        voltage. Diamonds with a different bias voltage will be un-
        selected.
        :param bias:
        :return:
        '''
        assert (type(bias) == t.IntType), 'Bias has to be int-Type'
        count = 0
        for run_number in self.run_numbers:
            if self.selection[run_number]:
                unselectrun = True
                if self.run_infos[run_number]['hv dia1'] == bias:
                    unselectrun = False
                else:
                    self.channels[run_number][0] = False
                if self.run_infos[run_number]['hv dia2'] == bias:
                    unselectrun = False
                else:
                    self.channels[run_number][3] = False
                if unselectrun:
                    self.selection[run_number] = False
                    count += 1
        self.make_log_entry('All Selected Runs unselected if not ' + str(bias) + 'V bias applied. Only ' + str(bias) + 'V Bias Runs left. -' + str(count) + ' selections')
        self.verbose_print('All Selected Runs unselected if not ' + str(bias) + 'V bias applied. Only ' + str(bias) + 'V Bias Runs left. -' + str(count) + ' selections')

    def select_run(self, run_number, do_assert=True, unselect=False):
        if do_assert:
            assert run_number in self.run_numbers, 'run {run} not found in list of run numbers. Check run_log json file!'.format(run=run_number)
        self.selection[run_number] = True if not unselect else False
        if unselect:
            self.reset_channels(run_number)

    def unselect_run(self, run_number, do_assert=True):
        self.select_run(run_number, do_assert, unselect=True)

    def ExcludeRuns(self, run_number):
        '''
        This method will un-select the run with number run_number.
        run_number has to be a single integer or a list of integers in
        order to exclude several runs.
        :param run_number:
        :return:
        '''
        assert (type(run_number) == t.IntType or type(run_number) == t.ListType), 'Wrong input type. run_number has to be either integer or list of integer'
        listOfRuns = self.get_selected_runs()
        if type(run_number) == t.IntType:
            if run_number in listOfRuns:
                self.unselect_run(run_number)
                self.make_log_entry('Run ' + str(run_number) + ' unselected. -1 selection')
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
                self.select_run(run_number)

    def UnSelectUnlessInRange(self, minrun, maxrun):
        '''
        Keeps only runs in range [minrun, .. , maxrun].
        :param minrun:
        :param maxrun:
        :return:
        '''
        for run_number in self.get_selected_runs():
            if not (run_number <= maxrun and run_number >= minrun):
                self.unselect_run(run_number)

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
        print len(self.get_selected_runs()), ' Runs Selected:'

        def multilinetext(text, width):
            length = len(text)
            word_wraps = length / int(width)
            separator = '\n '
            for i in xrange(word_wraps):
                text = text[:((i + 1) * width + len(separator) * i)] + separator + text[(i + 1) * width + len(separator) * i:]
            return text

        if show_allcomments:
            printstring = '{nr} {type} {dia1sel}{dia1} {hv1:>5} {dia2sel}{dia2} {hv2:>5} {rate}'
        else:
            printstring = '{nr} {type} {dia1sel}{dia1} {hv1:>5} {dia2sel}{dia2} {hv2:>5} {rate} {comment}'

        Record = namedtuple('Record', ['runnumber', 'type', 'dia1', 'bias1', 'dia2', 'bias2', 'rate', 'comment'])
        listitems = []
        for runnumber in self.get_selected_runs():
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
            nr='Nr.'.ljust(3),
            type='Type'.ljust(9),
            dia1='Dia 1'.ljust(7),
            dia1sel=' ',
            hv1='HV 1'.ljust(5),
            dia2='Dia 2'.ljust(7),
            dia2sel=' ',
            hv2='HV2 '.ljust(5),
            rate='Rate'.ljust(6),
            comment='Comment'
        )
        for item in listitems:
            if len(item.comment[0][:30]) > 0:
                marker = '* '
            else:
                marker = ''
            dia1sel = ' '
            dia2sel = ' '
            if self.channels[item.runnumber][0]: dia1sel = '*'
            if self.channels[item.runnumber][3]: dia2sel = '*'
            print printstring.format(
                nr=str(item.runnumber).ljust(3),
                type=item.type.ljust(9),
                dia1=item.dia1.ljust(7),
                dia1sel=dia1sel,
                hv1=int(item.bias1),
                dia2=item.dia2.ljust(7),
                dia2sel=dia2sel,
                hv2=int(item.bias2),
                rate=(str(item.rate) + ' kHz').rjust(8),
                comment=(marker + item.comment[0][:commentlength]).ljust(commentlength + 2)
            )
            if show_allcomments and len(item.comment[0][:30]) > 0:
                print '   -- COMMENT: -----------------------------------'
                print multilinetext(' ' + item.comment[0][:], 50)
                print '--------------------------------------------------'

    def show_run_info(self, runs=None, detailed=False):
        '''
        Prints all run infos from the selected runs to the console.
        :param runs:
        :return:
        '''

        if detailed:
            if runs is None:
                for runnumber in self.get_selected_runs():
                    self._printRunInfo(runnumber)
            else:
                if type(runs) is list:
                    for runnumber in runs:
                        self._printRunInfo(runnumber)
                elif type(runs) is int:
                    self._printRunInfo(runs)
                else:
                    print 'Wrong input type'
        else:
            if runs is None:
                self.__print_runinfo_header()
                for runnumber in self.get_selected_runs():
                    self.__print_runinfo(runnumber)
            else:
                # todo:
                print 'not yet implemented'

    def _printRunInfo(self, runnumber):
        print '--- RUN ', runnumber, ' ---'
        for key in self.run_infos[runnumber].keys():
            print '{key:>20}: {value}'.format(key=key, value=self.run_infos[runnumber][key])

    def __print_runinfo(self, run):
        dia1 = self.run_infos[run]['diamond 1']
        dia2 = self.run_infos[run]['diamond 2']
        flux = self.run_infos[run]['measured flux']
        type = self.run_infos[run]['type']
        print '{run}\t{type}\t{dia1}\t{dia2}\t{flux}'.format(run=run, dia1=dia1, dia2=dia2, flux=flux, type=type)

    @staticmethod
    def __print_runinfo_header():
        print 'run\ttype\tdia1\tdia2\tflux'

    def get_diamond_names(self):
        """
        :return: all diamond names from the logfile
        """
        names = self.get_runinfo_values('diamond 1')
        for name in self.get_runinfo_values('diamond 2'):
            if name not in names:
                names.append(name)
        return sorted(names)

    def show_diamond_names(self):
        print 'Diamondnames:'
        for name in self.get_diamond_names():
            print '  ' + name

    def show_run_types(self):
        print 'Types:'
        for type_ in self.get_runinfo_values('type'):
            print '  ' + type_

    def get_runinfo_values(self, key):
        '''
        :param key: key of run info
        :return: all different values of the run info dict
        '''
        values = []
        for run, info in self.run_infos.iteritems():
            value = info[key]
            if not value in values:
                values.append(value)
        return sorted(values)

    def get_selected_runs(self):
        """ :return: list of selected run numbers. """
        selected = []
        for runnumber in self.run_numbers:
            if self.selection[runnumber]:
                selected.append(runnumber)
        if not len(selected):
            print 'No runs selected!'
        return sorted(selected)

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
            if self.selection[runnumber]:
                dia1 = self.channels[runnumber][0]
                dia2 = self.channels[runnumber][3]
                diamonds = int(dia1) * (1 << 0) + int(dia2) * (1 << 1)
                if diamonds == 0: diamonds = 3
                selected.append(diamonds)
        if self.verbose and len(selected) == 0:
            print 'No Runs Selected'
        return selected

    def ShowRunPlan(self, detailed=True, type_='rate_scan', show_allcomments=False, commentlength=0):
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
        print 'RUN PLAN FOR TESTCAMPAIGN:', self.TESTCAMPAIGN
        if not detailed:
            for types_ in self.run_plan.keys():
                print types_, ' :'
                numbers = map(int, self.run_plan[types_].keys())
                numbers.sort()
                for planNr in numbers:
                    print '\t Nr. {nr} : {runs}'.format(nr=planNr, runs=self.run_plan[types_][str(planNr)])
        else:
            tmp_selections = copy.deepcopy(self.selection)
            tmp_channel_selections = copy.deepcopy(self.channels)
            tmp_selectionLog = copy.deepcopy(self.logs)

            numbers = map(int, self.run_plan[type_].keys())
            numbers.sort()
            for i in numbers:
                self.UnselectAll()
                self.SelectRunsFromRunPlan(i, type_=type_)
                print '-----------------------------------------'
                print 'RUN PLAN ', i, ' ({typ}) : '.format(typ=type_)
                print '-----------------------------------------'
                self.ShowSelectedRuns(show_allcomments=show_allcomments, commentlength=commentlength)
                print '\n'
            self.selection = copy.deepcopy(tmp_selections)
            self.channels = copy.deepcopy(tmp_channel_selections)
            self.logs = copy.deepcopy(tmp_selectionLog)

    def SelectRunsFromRunPlan(self, number, type_='rate_scan'):
        '''
        Selects all runs corresponding to the run plan with key 'number'.
        :param number:
        :param type_:
        :return:
        '''
        runs = self.run_plan[type_][str(number)]
        self.run_plan = number
        self.SelectRuns(list_of_runs=runs)

    def AddSelectedRunsToRunPlan(self, key, run_type='rate_scan'):
        '''
        Saves all selected runs as a run plan with name 'key'. This
        run plan can later be imported via
            .SelectRunsFromRunPlan('key')
        :param key:
        :param run_type:
        :return:
        '''
        assert (run_type in ['rate_scan', 'voltage_scan', 'test'])
        if not type(key) is t.StringType:
            key = str(key)
        if self.runplans[self.TESTCAMPAIGN].has_key(run_type):
            self.runplans[self.TESTCAMPAIGN][run_type][key] = self.get_selected_runs()
        else:
            self.runplans[self.TESTCAMPAIGN][run_type] = {}
            self.runplans[self.TESTCAMPAIGN][run_type][key] = self.get_selected_runs()
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
        assert (type(list_of_runs) is t.ListType), 'list_of_runs not a list'
        for runnumber in list_of_runs:
            self.select_run(runnumber)

if __name__ == '__main__':
    z = RunSelection()