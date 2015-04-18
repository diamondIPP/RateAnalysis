from RunClass import Run
import types as t


data_types = {
    0 : "DATA",
    1 : "PEDESTAL",
    2 : "VOLTAGESCAN",
    3 : "LONG_RUN",
    4 : "OTHER",
  }

class RunSelection(Run):

    def __init__(self, verbose = False):
        Run.__init__(self, validate=False, run_number=None, verbose=verbose)
        self.run_numbers = self.runinfo.keys()
        self.LoadRuns()
        self.InitializeSelections()
        self.SelectionLog = {}

    def __str__(self):
        nr = len(self.run_numbers)
        selected_runs = self.GetSelectedRuns()
        return "RunSelection Object\n"+str(len(selected_runs))+" Out of "+str(nr)+" runs selected. Selections made:"+self.Log()

    def Log(self, event=None):
        if event is not None:
            i = len(self.SelectionLog)
            self.SelectionLog[i] = event
        string = '\n'
        for i in self.SelectionLog.keys():
            string += str(i+1)+'.) '+self.SelectionLog[i]+'\n'
        return string

    def LoadRuns(self):
        self.runs = {}
        for runnumber in self.run_numbers:
            self.SetRun(runnumber, validate=False)
            self.runs[runnumber] = self.current_run

    def InitializeSelections(self):
        self.selections = {}
        self.SelectionLog = {}
        for runnumber in self.run_numbers:
            self.selections[runnumber] = False

    def SelectAll(self):
        for run_number in self.run_numbers:
            self.selections[run_number] = True
        self.Log('All runs selected')
        self.VerbosePrint('All runs selected')

    def UnSelectAll(self):
        self.InitializeSelections()
        self.VerbosePrint('All runs unselected')

    def SelectDataType(self, data_type): # data type
        assert(type(data_type) is t.IntType and 0<=data_type<=4), "wrong data type"
        count = 0
        for run_number in self.run_numbers:
            if self.runs[run_number]['data_type'] == data_type:
                self.selections[run_number] = True
                count += 1
        self.Log('Runs of Type '+data_types[data_type]+' selected. +'+str(count)+'selections')
        self.VerbosePrint('Runs of Type '+data_types[data_type]+' selected. +'+str(count)+' selections')

    def SelectDiamondRuns(self, diamondname):
        count = 0
        for run_number in self.run_numbers:
            if self.runs[run_number]['diamond'] == diamondname:
                self.selections[run_number] = True
                count += 1
        self.Log('Runs with diamond '+diamondname+' selected. +'+str(count)+' selections')
        self.VerbosePrint('Runs with diamond '+diamondname+' selected. +'+str(count)+' selections')

    def SelectIrradiationRuns(self, irradiated=True, irrtype=None):
        count = 0
        if irradiated and irrtype is None:
            for run_number in self.run_numbers:
                self.SetRun(run_number, validate=False)
                if self.diamond.Specifications['Irradiation'] == 'proton' or self.diamond.Specifications['Irradiation'] == 'neutron':
                    self.selections[run_number] = True
                    count += 1
            self.Log('Irradiated Runs selected')
        elif not irradiated:
            for run_number in self.run_numbers:
                self.SetRun(run_number, validate=False)
                if self.diamond.Specifications['Irradiation'] == 'no':
                    self.selections[run_number] = True
                    count += 1
            self.Log('Non-Irradiated Runs selected')
        else:
            assert(irrtype in ['no', 'proton', 'neutron']), "wrong irradiation type. Choose irrtype in [`proton`, `neutron`, `no`]"
            if irrtype == 'no':
                self.SelectIrradiationRuns(irradiated=False)
            else:
                for run_number in self.run_numbers:
                    self.SetRun(run_number, validate=False)
                    if self.diamond.Specifications['Irradiation'] == irrtype:
                        self.selections[run_number] = True
                        count += 1
                self.Log('Irradiated Runs selected with '+irrtype+' irradiation. +'+str(count)+' selections')
                self.VerbosePrint('Irradiated Runs selected with '+irrtype+' irradiation. +'+str(count)+' selections')

    def UnSelectUnlessDataType(self, data_type):
        assert(type(data_type) is t.IntType and 0<=data_type<=4), "wrong data type"
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if self.runs[run_number]['data_type'] == data_type:
                    pass
                else:
                    self.selections[run_number] = False
                    count += 1
        self.Log('All Selected Runs unselected if not of Type '+data_types[data_type]+'. -'+str(count)+' selections')
        self.VerbosePrint('All Selected Runs unselected if not of Type '+data_types[data_type]+'. -'+str(count)+' selections')

    def UnSelectUnlessIrradiation(self, irradiated=True, irrtype=None):
        count = 0
        if irradiated and irrtype is None:
            for run_number in self.run_numbers:
                if self.selections[run_number]:
                    self.SetRun(run_number, validate=False)
                    if self.diamond.Specifications['Irradiation'] == 'proton' or self.diamond.Specifications['Irradiation'] == 'neutron':
                        pass
                    else:
                        self.selections[run_number] = False
                        count += 1
            self.Log('All Selected Runs unselected if non-irradiated. Only radiated Runs left. -'+str(count)+' selections')
            self.VerbosePrint('All Selected Runs unselected if non-irradiated. Only radiated Runs left. -'+str(count)+' selections')
        elif not irradiated:
            for run_number in self.run_numbers:
                if self.selections[run_number]:
                    self.SetRun(run_number, validate=False)
                    if self.diamond.Specifications['Irradiation'] == 'no':
                        pass
                    else:
                        self.selections[run_number] = False
                        count += 1
            self.Log('All Selected Runs unselected if irradiated. Only non-radiated Runs left. -'+str(count)+' selections')
            self.VerbosePrint('All Selected Runs unselected if irradiated. Only non-radiated Runs left. -'+str(count)+' selections')
        else:
            assert(irrtype in ['no', 'proton', 'neutron']), "wrong irradiation type. Choose irrtype in [`proton`, `neutron`, `no`]"
            if irrtype == 'no':
                self.UnSelectUnlessIrradiation(irradiated=False)
            else:
                for run_number in self.run_numbers:
                    if self.selections[run_number]:
                        self.SetRun(run_number, validate=False)
                        if self.diamond.Specifications['Irradiation'] == irrtype:
                            pass
                        else:
                            self.selections[run_number] = False
                            count += 1
                self.Log('All Selected Runs unselected if not radiated by '+irrtype+'. Only '+irrtype+'-radiated Runs left. -'+str(count)+' selections')
                self.VerbosePrint('All Selected Runs unselected if not radiated by '+irrtype+'. Only '+irrtype+'-radiated Runs left. -'+str(count)+' selections')

    def UnSelectUnlessDiamond(self, diamondname):
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if self.runs[run_number]['diamond'] == diamondname:
                    pass
                else:
                    self.selections[run_number] = False
                    count += 1
        self.Log('All Selected Runs unselected if not using '+diamondname+' diamond. Only '+diamondname+' diamond Runs left. -'+str(count)+' selections')
        self.VerbosePrint('All Selected Runs unselected if not using '+diamondname+' diamond. Only '+diamondname+' diamond Runs left. -'+str(count)+' selections')

    def UnSelectUnlessBias(self, bias):
        assert(type(bias) is t.IntType), "Bias has to be int-Type"
        count = 0
        for run_number in self.run_numbers:
            if self.selections[run_number]:
                if self.runs[run_number]['bias_voltage'] == bias:
                    pass
                else:
                    self.selections[run_number] = False
                    count += 1
        self.Log('All Selected Runs unselected if not '+str(bias)+'V bias applied. Only '+str(bias)+'V Bias Runs left. -'+str(count)+' selections')
        self.VerbosePrint('All Selected Runs unselected if not '+str(bias)+'V bias applied. Only '+str(bias)+'V Bias Runs left. -'+str(count)+' selections')

    def ValidateSelectedRuns(self):
        self.ValidateRuns(self.GetSelectedRuns())

    def ShowSelectedRuns(self, get_string=False):
        string = 'Detailed View of Selected Runs:\n'
        string += 'Nr.\tDiamond\t\tCrystal\tIrr.\tType\t\tBias\trate_ps\trate_raw\trate_trigger\n'
        for runnumber in self.GetSelectedRuns():
            self.SetRun(runnumber, validate=False)
            diamond = self.runs[runnumber]['diamond']
            irr = self.diamond.Specifications['Irradiation']
            crystal = self.diamond.Specifications['CrystalStructure']
            data_type = self.runs[runnumber]['data_type']
            bias_voltage = self.runs[runnumber]['bias_voltage']
            rate_ps = self.runs[runnumber]['rate_ps']
            rate_raw = self.runs[runnumber]['rate_raw']
            rate_trigger = self.runs[runnumber]['rate_trigger']
            if data_type == 1 or data_type == 2 or data_type == 3:
                dt_tab = '\t'
            else:
                dt_tab = '\t\t'
            if diamond == 'S129-old-box':
                d_tab = '\t'
            else:
                d_tab = '\t\t'
            string += str(runnumber)+'\t'+diamond+d_tab+crystal+'\t'+irr+'\t'+data_types[data_type]+dt_tab+str(bias_voltage)+'\t'+str(rate_ps)+'\t'+str(rate_raw)+'\t\t'+str(rate_trigger)+'\n'
        if not get_string:
            print string
        else:
            return string

    def GetSelectedRuns(self):
        selected = []
        for runnumber in self.run_numbers:
            if self.selections[runnumber]:
                selected.append(runnumber)
        if self.verbose and len(selected) == 0:
            print "No Runs Selected"
        return selected
