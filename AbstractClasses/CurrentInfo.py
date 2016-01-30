# ====================================
# IMPORTS
# ====================================
from glob import glob
from datetime import datetime, timedelta
from configparser import ConfigParser
from Elementary import Elementary


# ====================================
# CLASS FOR THE DATA
# ====================================
class KeithleyInfo(Elementary):
    """reads in information from the keithley log file"""

    def __init__(self, analysis, averaging=False, points=10):
        Elementary.__init__(self)

        self.DoAveraging = averaging
        self.Points = points
        self.analysis = analysis

        # analysis/run info
        self.RunInfo = analysis.run.RunInfo
        self.RunNumber = analysis.run_number
        self.StartEvent = analysis.StartEvent
        self.EndEvent = analysis.EndEvent
        self.StartTime = datetime.strptime(self.RunInfo['start time'], '%Y-%m-%dT%H:%M:%SZ') + timedelta(hours=1)
        self.StopTime = datetime.strptime(self.RunInfo['stop time'], '%Y-%m-%dT%H:%M:%SZ') + timedelta(hours=1)
        self.Channel = analysis.channel
        self.DiamondName = analysis.diamond_name

        # config
        self.RunConfigParser = analysis.run.run_config_parser
        self.ConfigParser = self.load_parser()

        # device info
        self.Number = self.get_device_nr()
        self.Channel = self.get_device_channel()
        self.Brand = self.ConfigParser.get('HV' + self.Number, 'name').split('-')[0].strip('0123456789')
        self.Model = self.ConfigParser.get('HV' + self.Number, 'model')
        self.Name = '{0} {1}'.format(self.Brand, self.Model)

        self.DataPath = self.find_data_path()
        self.LogNames = self.get_logs_from_start()

        # data
        self.Currents = []
        self.Voltages = []
        self.Time = []

        # self.log_names = self.logs_from_start()
        # self.mean_curr = 0
        # self.mean_volt = 0
        # self.update_data()

    # ==========================================================================
    # region INIT
    def load_parser(self):
        parser = ConfigParser()
        parser.read(self.RunConfigParser.get('BASIC', 'hvconfigfile'))
        return parser

    def get_device_nr(self):
        full_str = self.RunInfo['dia{dia}supply'.format(dia=1 if not self.Channel else 2)]
        return str(full_str.split('-')[0])

    def get_device_channel(self):
        full_str = self.RunInfo['dia{dia}supply'.format(dia=1 if not self.Channel else 2)]
        return full_str.split('-')[1] if len(full_str) > 1 else '0'

    def find_data_path(self):
        hv_datapath = self.RunConfigParser.get('BASIC', 'hvdatapath')
        return '{data}{dev}_CH{ch}/'.format(data=hv_datapath, dev=self.ConfigParser.get('HV' + self.Number, 'name'), ch=self.Channel)
    # endregion

    def get_logs_from_start(self):
        log_names = sorted([name for name in glob(self.DataPath + '*')])
        start_log = None
        for i, name in enumerate(log_names):
            log_date = self.get_log_date(name)
            if log_date >= self.StartTime:
                break
            start_log = i
        return log_names[start_log:]

    @staticmethod
    def get_log_date(name):
        log_date = name.split('/')[-1].split('_')
        log_date = ''.join(log_date[3:])
        return datetime.strptime(log_date, '%Y%m%d%H%M%S.log')

    def find_data(self):
        stop = False
        for i, name in enumerate(self.LogNames):
            self.mean_curr = 0
            self.mean_volt = 0
            log_date = self.get_log_date(name)
            data = open(name, 'r')
            # jump to the correct line of the first file
            if not i:
                self.find_start(data, log_date)
            index = 0
            for line in data:
                info = line.split()
                if self.isfloat(info[1]) and len(info) > 2:
                    now = datetime.strptime(log_date.strftime('%Y%m%d') + info[0], '%Y%m%d%H:%M:%S')
                    if self.StartTime < now < self.StopTime and float(info[2]) < 1e30:
                        self.save_data(now, info, index)
                        index += 1
                    if self.StopTime < now:
                        stop = True
                        break
            data.close()
            if stop:
                break

    def save_data(self, now, info, index, shifting=False):
        total_seconds = (now - datetime(now.year, 1, 1)).total_seconds()
        if self.StartTime < now < self.StopTime and float(info[2]) < 1e30:
            index += 1
            if self.DoAveraging:
                if not shifting:
                    self.mean_curr += float(info[2]) * 1e9
                    self.mean_volt += float(info[1])
                    if index % self.Points == 0:
                        self.Currents.append(self.mean_curr / self.Points)
                        self.Time.append(total_seconds)
                        self.Voltages.append(self.mean_volt / self.Points)
                        self.mean_curr = 0
                        self.mean_volt = 0
                # else:
                #     if index <= self.Points:
                #         self.mean_curr += float(info[2]) * 1e9
                #         dicts[1][key].append(self.mean_curr / index)
                #         if index == self.Points:
                #             self.mean_curr /= self.Points
                #     else:
                #         mean_curr = self.mean_curr * weight + (1 - weight) * float(info[2]) * 1e9
                #         dicts[1][key].append(mean_curr)
                #     dicts[0][key].append(convert_time(now))
                #     dicts[2][key].append(float(info[1]))
            else:
                self.Currents.append(float(info[2]) * 1e9)
                self.Time.append(total_seconds)
                self.Voltages.append(float(info[1]))

    def find_start(self, data, log_date):
        lines = len(data.readlines())
        data.seek(0)
        if lines < 10000:
            return
        was_lines = 0
        for i in range(6):
            lines /= 2
            for j in xrange(lines):
                data.readline()
            while True:
                info = data.readline().split()
                if not info:
                    break
                if self.isfloat(info[1]):
                    now = datetime.strptime(log_date.strftime('%Y%m%d') + info[0], '%Y%m%d%H:%M:%S')
                    if now < self.StartTime:
                        was_lines += lines
                        break
                    else:
                        data.seek(0)
                        for k in xrange(was_lines):
                            data.readline()
                        break

    def convert_to_relative_time(self):
        zero = self.Time[0]
        for i in xrange(len(self.Time)):
            self.Time[i] = self.Time[i] - zero
