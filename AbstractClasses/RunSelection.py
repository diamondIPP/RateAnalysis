from RunClass import Run
import types as t
from collections import namedtuple
import json
import copy
from ROOT import TString


class RunSelection(Run):

    def __init__(self, verbose = False):
        self.run_numbers = []
        self._selectionLog = {}
        Run.__init__(self, validate=False, run_number=None, verbose=verbose)
        self.run_numbers = [int(self.allRunKeys[i][-3:]) for i in xrange(len(self.allRunKeys))] # list of all the run numbers found in json file
        self.run_numbers.sort()
        self._GetRunPlanFromFile()
        self._LoadRuns()
        self._InitializeSelections()

    def __str__(self):
        nr = len(self.run_numbers)
        selected_runs = self.GetSelectedRuns()
        return "RunSelection Object\n"+str(len(selected_runs))+" Out of "+str(nr)+" runs selected. Selections made:"+self._Log()

    def _Log(self, event=None):
        if event is not None:
            i = len(self._selectionLog)
            self._selectionLog[i] = event
        string = '\n'
        for i in self._selectionLog.keys():
            string += str(i+1)+'.) '+self._selectionLog[i]+'\n'
        return string

    def _LoadRuns(self):
        '''
        loads all the run infos in a dict with the run numbers as keys
        :return:
        '''
        self.runs = {}
        for runnumber in self.run_numbers:
            self.SetRun(runnumber, validate=False, loadROOTFile = False)
            self.runs[runnumber] = self.current_run

    def _InitializeSelections(self):
        '''
        creates a dict of bools to store selections. dict is filled with False (no run selected)
        the selections log is created/cleared
        :return:
        '''
        self.selections = {}
        self.channel_selections = {}
        self._selectionLog = {}
        for runnumber in self.run_numbers:
            self.selections[runnumber] = False
            self.channel_selections[runnumber] = {
                0 : False,
                3 : False
            }

    def _GetRunPlanFromFile(self):
        runplanfile = open(self.runplaninfofile, "r")
        self.runplans = json.load(runplanfile)
        runplanfile.close()
        campaigns = self.runplans.keys()
        if not self.TESTCAMPAIGN in campaigns:
            self.runplans[self.TESTCAMPAIGN] = {}
            self._SaveRunPlanInfoFile()
        self.runplan = self.runplans[self.TESTCAMPAIGN]

    def _SaveRunPlanInfoFile(self):
        f = open(self.runplaninfofile, "w")
        json.dump(self.runplans, f, indent=2, sort_keys=True)
        f.close()
        self._GetRunPlanFromFile()

    def SelectAll(self, selectDiamond1=True, selectDiamond2=True):
        for run_number in self.run_numbers:
            self.selections[run_number] = True
            if selectDiamond1: self.channel_selections[run_number][0] = True
            if selectDiamond2: self.channel_selections[run_number][0] = True
        self._Log('All runs selected')
        self.VerbosePrint('All runs selected')

    def UnSelectAll(self):
        self._InitializeSelections()
        self.VerbosePrint('All runs unselected')

    def UnselectAll(self):
        self.UnSelectAll()

    def SetChannels(self, diamond1=True, diamond2=True):
        '''
        Sets the channels (diamonds) of the selected runs to active or inactive
        :param diamond1:
        :param diamond2:
        :return:
        '''
        for run_number in self.GetSelectedRuns():
            self.channel_selections[run_number][0] = diamond1
            self.channel_selections[run_number][3] = diamond2
        change1 = "diamond 1 = "+str(diamond1)
        change2 = "diamond 2 = "+str(diamond2)
        self._Log("Channels of selected runs set: "+change1+" "+change2)

    def SelectDataType(self, data_type): # data type
        '''
        Selects just the runs, NOT the channels (diamonds)!
        To select the channels use SelectChannels() or SelectDiamonds()
        :param data_type:
        :return:
        '''
        types = self._GetValues("type")
        assert(data_type in types), "wrong data type. \n\tSelect type from: "+str(types)
        count = 0
        for run_number in self.run_numbers:
            if self.runs[run_number]['type'] == data_type:
                self._SelectRun(run_number)
                count += 1
        self._Log('Runs of Type '+data_type+' selected. +'+str(count)+' selections')
        self.VerbosePrint('Runs of Type '+data_type+' selected. +'+str(count)+' selections')

    def SelectDiamondRuns(self, diamondname, only_selected_runs=False):
        diamondnames = self.ShowDiamondNames(getNames=True)
        assert(diamondname in diamondnames), "wrong diamond name. \n\tSelect diamond name from: "+str(diamondnames)
        if not only_selected_runs:
            choice = self.run_numbers
        else:
            choice = self.GetSelectedRuns()
        count = 0
        for run_number in choice:
            if self.runs[run_number]['diamond 1'] == diamondname:
                self.selections[run_number] = True
                self.channel_selections[run_number][0] = True
                count += 1
            if self.runs[run_number]['diamond 2'] == diamondname:
                self.selections[run_number] = True
                self.channel_selections[run_number][3] = True
                count += 1
        self._Log('Runs and Channels containing '+diamondname+' selected. +'+str(count)+' runs selected')
        self.VerbosePrint('Runs and Channels containing '+diamondname+' selected. +'+str(count)+' runs selected')

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
        types = self._GetValues("type")
        assert(data_type in types), "wrong data type. \n\tSelect type from: "+str(types)
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if self.runs[run_number]['type'] == data_type:
                    pass
                else:
                    self._UnselectRun(run_number)
                    count += 1
        self._Log('All Selected Runs unselected if not of Type '+data_type+'. -'+str(count)+' selections')
        self.VerbosePrint('All Selected Runs unselected if not of Type '+data_type+'. -'+str(count)+' selections')

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
        diamondnames = self.ShowDiamondNames(True)
        assert(diamondname in diamondnames), "wrong diamond name. \n\tSelect diamond name from: "+str(diamondnames)
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if diamondname in [self.runs[run_number]['diamond 1'], self.runs[run_number]['diamond 2']]:
                    pass
                else:
                    self._UnselectRun(run_number)
                    count += 1
        self._Log('All Selected Runs unselected if not using '+diamondname+' diamond. Only runs countaining '+diamondname+' left. -'+str(count)+' selections')
        self.VerbosePrint('All Selected Runs unselected if not using '+diamondname+' diamond. Only runs countaining '+diamondname+' left. -'+str(count)+' selections')

    def UnSelectUnlessBias(self, bias):
        assert(type(bias) == t.IntType), "Bias has to be int-Type"
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                unselectrun = True
                if self.runs[run_number]['hv dia1'] == bias:
                    unselectrun = False
                else:
                    self.channel_selections[run_number][0] = False
                if self.runs[run_number]['hv dia2'] == bias:
                    unselectrun = False
                else:
                    self.channel_selections[run_number][3] = False
                if unselectrun:
                    self.selections[run_number] = False
                    count += 1
        self._Log('All Selected Runs unselected if not '+str(bias)+'V bias applied. Only '+str(bias)+'V Bias Runs left. -'+str(count)+' selections')
        self.VerbosePrint('All Selected Runs unselected if not '+str(bias)+'V bias applied. Only '+str(bias)+'V Bias Runs left. -'+str(count)+' selections')

    def _SelectRun(self, run_number):
        assert(run_number in self.run_numbers), "run number "+str(run_number)+" not found in list of run numbers. Check run_log json file!"
        self.selections[run_number] = True

    def _UnselectRun(self, run_number):
        self.selections[run_number] = False
        self.channel_selections[run_number][0] = False
        self.channel_selections[run_number][3] = False

    def ExcludeRuns(self, run_number):
        assert(type(run_number) == t.IntType or type(run_number) == t.ListType), "Wrong input type. run_number has to be either integer or list of integer"
        listOfRuns = self.GetSelectedRuns()
        if type(run_number) == t.IntType:
            if run_number in listOfRuns:
                self._UnselectRun(run_number)
                self._Log("Run "+str(run_number)+" unselected. -1 selection")
        else:
            ListToExclude = run_number
            for run_number in ListToExclude:
                self.ExcludeRuns(int(run_number))

    def SelectRunsInRange(self, minrun, maxrun):
        for run_number in self.run_numbers:
            if run_number <= maxrun and run_number >= minrun:
                self._SelectRun(run_number)

    def UnSelectUnlessInRange(self, minrun, maxrun):
        for run_number in self.GetSelectedRuns():
            if not (run_number <= maxrun and run_number >= minrun):
                self._UnselectRun(run_number)

    # def ValidateSelectedRuns(self):
    #     #self.ValidateRuns(self.GetSelectedRuns())
    #     pass

    def ShowSelectedRuns(self, show_allcomments=False, commentlength = 15):
        print len(self.GetSelectedRuns()), " Runs Selected:"
        def multilinetext(text, width):
            length = len(text)
            word_wraps = length/int(width)
            separator = "\n "
            for i in xrange(word_wraps):
                text = text[:((i+1)*width+len(separator)*i)]+separator+text[(i+1)*width+len(separator)*i:]
            return text

        if show_allcomments:
            printstring = "{nr} {type} {dia1sel}{dia1} {hv1:>5} {dia2sel}{dia2} {hv2:>5} {rate}"
        else:
            printstring = "{nr} {type} {dia1sel}{dia1} {hv1:>5} {dia2sel}{dia2} {hv2:>5} {rate} {comment}"

        Record = namedtuple("Record",["runnumber", "type", "dia1", "bias1", "dia2", "bias2", "rate", "comment"])
        listitems = []
        for runnumber in self.GetSelectedRuns():
            listitems += [Record(
                runnumber,
                self.runs[runnumber]['type'],
                self.runs[runnumber]['diamond 1'],
                self.runs[runnumber]['hv dia1'],
                self.runs[runnumber]['diamond 2'],
                self.runs[runnumber]['hv dia2'],
                int(self.runs[runnumber]['measured flux']),
                [self.runs[runnumber]['user comments']]
                )]

        print printstring.format(
                nr = "Nr.".ljust(3),
                type = "Type".ljust(9),
                dia1 = "Dia 1".ljust(7),
                dia1sel = " ",
                hv1 = "HV 1".ljust(5),
                dia2 = "Dia 2".ljust(7),
                dia2sel=" ",
                hv2 = "HV2 ".ljust(5),
                rate = "Rate".ljust(6),
                comment = "Comment"
            )
        for item in listitems:
            if len(item.comment[0][:30])>0:
                marker = "* "
            else:
                marker = ""
            dia1sel = " "
            dia2sel = " "
            if self.channel_selections[item.runnumber][0]: dia1sel = "*"
            if self.channel_selections[item.runnumber][3]: dia2sel = "*"
            print printstring.format(
                nr = str(item.runnumber).ljust(3),
                type = item.type.ljust(9),
                dia1 = item.dia1.ljust(7),
                dia1sel = dia1sel,
                hv1 = int(item.bias1),
                dia2 = item.dia2.ljust(7),
                dia2sel=dia2sel,
                hv2 = int(item.bias2),
                rate = (str(item.rate)+" kHz").rjust(8),
                comment = (marker+item.comment[0][:commentlength]).ljust(commentlength+2)
            )
            if  show_allcomments and len(item.comment[0][:30])>0:
                print "   -- COMMENT: -----------------------------------"
                print multilinetext(" "+item.comment[0][:], 50)
                print "--------------------------------------------------"

    def ShowRunInfo(self, runs=[]):
        if runs == []:
            for runnumber in self.GetSelectedRuns():
                self._printRunInfo(runnumber)
        else:
            if type(runs) is t.ListType:
                for runnumber in runs:
                    self._printRunInfo(runnumber)
            elif type(runs) is t.IntType:
                self._printRunInfo(runs)
            else:
                print "Wrong input type"

    def _printRunInfo(self, runnumber):
        print "--- RUN ", runnumber, " ---"
        for key in self.runs[runnumber].keys():
            print "{key:>20}: {value}".format(key=key, value=self.runs[runnumber][key])

    def ShowDiamondNames(self, getNames=False):
        diamondnames = []
        diamond1names = self._GetValues("diamond 1")
        diamond2names = self._GetValues("diamond 2")
        biglist = diamond1names+diamond2names
        for name in biglist:
            if not name in diamondnames:
                diamondnames += [name]
        if not getNames:
            print "Diamondnames:"
            for i in xrange(len(diamondnames)):
                print "\t"+diamondnames[i]
        else:
            return diamondnames

    def ShowTypes(self):
        types = self._GetValues("type")
        print "Types:"
        for i in types:
            print "\t"+i

    def _GetValues(self, key):
        valuelist = []
        for runnumber in self.runs.keys():
            value = self.runs[runnumber][key]
            if not value in valuelist:
                valuelist += [value]
        return valuelist

    def GetSelectedRuns(self):
        selected = []
        for runnumber in self.run_numbers:
            if self.selections[runnumber]:
                selected.append(runnumber)
        if self.verbose and len(selected) == 0:
            print "No Runs Selected"
        return selected

    def GetSelectedDiamonds(self):
        selected = []
        for runnumber in self.run_numbers:
            if self.selections[runnumber]:
                dia1 = self.channel_selections[runnumber][0]
                dia2 = self.channel_selections[runnumber][3]
                diamonds = int(dia1)*(1<<0) + int(dia2)*(1<<1)
                if diamonds == 0: diamonds = 3
                selected.append(diamonds)
        if self.verbose and len(selected) == 0:
            print "No Runs Selected"
        return selected

    def ShowRunPlan(self, detailed=True, type_="rate_scan", show_allcomments=False, commentlength=0):
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
            tmp_selectionLog = copy.deepcopy(self._selectionLog)


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
            self._selectionLog = copy.deepcopy(tmp_selectionLog)

    def SelectRunsFromRunPlan(self, number, type_="rate_scan"):
        runs = self.runplan[type_][str(number)]
        self.SelectRuns(list_of_runs=runs)

    def AddSelectedRunsToRunPlan(self, key, run_type="rate_scan"):
        assert(run_type in ["rate_scan", "voltage_scan", "test"])
        if not type(key) is t.StringType:
            key = str(key)
        if self.runplans[self.TESTCAMPAIGN].has_key(run_type):
            self.runplans[self.TESTCAMPAIGN][run_type][key] = self.GetSelectedRuns()
        else:
            self.runplans[self.TESTCAMPAIGN][run_type] = {}
            self.runplans[self.TESTCAMPAIGN][run_type][key] = self.GetSelectedRuns()
        self._SaveRunPlanInfoFile()

    def SelectRuns(self, list_of_runs, select_dia1=False, select_dia2=False):
        assert(type(list_of_runs) is t.ListType), "list_of_runs not a list"
        for runnumber in list_of_runs:
            self._SelectRun(runnumber)
