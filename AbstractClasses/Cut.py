from AbstractClasses.Elementary import Elementary
import ConfigParser


class Cut(Elementary):

    def __init__(self, parent_analysis, ):

        self.analysis = parent_analysis
        self.userCutTypes = {
            "EventCut": "",
            "ExcludeFirst": "",
            "FFT": ""

        }
        self._cutTypes = {

        }
        self.excludefirst = 0
        self.cut1 = ""
        self.cut2 = ""
        self.cut = {
            0: "", # cut diamond 1
            3: ""  # cut diamond 2
        }
        self.pulseHeightDefinition = {
            0: "", # pulseHeightDefinition diamond 1
            3: ""  # pulseHeightDefinition diamond 2
        }
        Elementary.__init__()

    def LoadConfig(self):
        configfile = "Configuration/AnalysisConfig_"+self.TESTCAMPAIGN+".cfg"
        parser = ConfigParser.ConfigParser()
        parser.read(configfile)

        self.cut[0] = parser.get("BASIC", "cut")
        self.cut[3] = parser.get("BASIC", "cut")
        self._signalname = parser.get("BASIC", "signalname")
        self.pedestal_correction = parser.getboolean("BASIC", "pedestal_correction")
        self._pedestalname = parser.get("BASIC", "pedestalname")
        self.excludefirst = parser.getint("BASIC", "excludefirst")
        self.excludeBeforeJump = parser.getint("BASIC", "excludeBeforeJump")
        self.excludeAfterJump = parser.getint("BASIC", "excludeAfterJump")

        channels = [0,3]
        for i in channels:
            self.cut[i] = self.cut[i].format(channel=1)
            if not self.pedestal_correction:
                self.pulseHeightDefinition[i] = (self._signalname+"[{channel}]").format(channel=i)
                self.userSignalDefinition = self._signalname
            else:
                self.pulseHeightDefinition[i] = (self._signalname+"[{channel}]-"+self._pedestalname+"[{channel}]").format(channel=i)
                self.userSignalDefinition = self._signalname+"-"+self._pedestalname

    def GetCut(self, channel):
        return self.cut[channel]

    def GetPulseHeightDefinition(self, channel):
        return self.pulseHeightDefinition[channel]

    def GetUserCutString(self):
        string_ = ""
        return string_

    def SetCut(self):