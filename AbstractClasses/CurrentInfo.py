# ====================================
# IMPORTS
# ====================================
import glob
from datetime import datetime
from collections import OrderedDict
from RunClass import Run
from functions1 import *

weight = 0.93
# ====================================
# CLASS FOR THE DATA
# ====================================
class KeithleyInfo():
    """reads in information from the keithley log file"""

    def __init__(self, log, jsonfile, start, stop, number, averaging, points):
        self.RunInfo =

        self.single_mode = (True if number == "1" else False)
        self.do_averaging = (True if averaging else False)
        self.points = int(points) if is_float(points) else 10
        self.number = number
        if start == -1 and stop == -1:

            RunInfo.__init__(self, jsonfile, start, time()+24)
        else:
            RunInfo.__init__(self, jsonfile, start, stop)
        self.keithleys = OrderedDict([("Keithley1", "Silicon"),
                                      ("Keithley2", str(self.dia1)),
                                      ("Keithley3", str(self.dia2))])
        self.log_dir = [str(log[0]) + "/HV*log", str(log[1]) + "/HV*log"]
        self.keithley_name = [self.get_keithley_name(0), self.get_keithley_name(1)]
        if self.single_mode:
            self.keithleys = OrderedDict([("Keithley1", self.keithley_name[0])])
        elif number == "2":
            self.keithleys = OrderedDict([("Keithley1", self.keithley_name[0]),
                                          ("Keithley2", self.keithley_name[1])])
        self.dias = self.make_dict(self.dia1, self.dia2)
        self.log_names = self.logs_from_start()
        self.mean_curr = 0
        self.mean_volt = 0
        self.update_data()

    def update_data(self):
        self.data = self.find_data()
        self.time_x = self.data[0]
        self.current_y = self.data[1]
        self.voltage_y = self.data[2]

    def get_keithley_name(self, num):
        name = self.log_dir[num].split('/')
        for i in name:
            if i.lower().startswith('keithley') and not i.lower().endswith('client'):
                name = i
                break
        return str(name)

    def get_lognames(self, num):
        log_names = []
        for name in glob.glob(self.log_dir[num]):
            log_names.append(name)
        log_names = sorted(log_names)
        return log_names

    def get_start_log(self, num):
        start_log = 0
        log_names = self.get_lognames(num)
        for i in range(len(log_names)):
            name = log_names[i].strip('.log').split('/')
            name = name[-1].split('_')
            log_date = ""
            for j in range(3, 9):
                log_date += name[j] + " "
            log_date = log_date.strip(' ')
            log_date = datetime.strptime(log_date, "%Y %m %d %H %M %S")
            if log_date >= self.start:
                break
            start_log = i
        return start_log

    def logs_from_start(self):
        log_names = {}
        ind = 0
        for key in self.keithleys:
            log_names[key] = []
            for i in range(self.get_start_log(ind), len(self.get_lognames(ind))):
                log_names[key].append(self.get_lognames(ind)[i])
            ind +=1
        return log_names

    def create_dicts(self):
        dicts = {}
        for key in self.keithleys:
            dicts[key] = []
        return dicts

    @staticmethod
    def get_log_date(name):
        name = name.strip('.log').split('/')
        name = name[-1].split('_')
        log_date = ""
        for i in range(3, 6):
            log_date += name[i] + "-"
        log_date = log_date.strip('-')
        return log_date

    def find_data(self):
        self.log_names = self.logs_from_start()
        dicts = [self.create_dicts(), self.create_dicts(), self.create_dicts()]
        stop = [0,0,0]
        ind = 0
        stop_ind = 0
        for key in self.keithleys:
            for name in self.log_names[key]:
                self.mean_curr = 0
                self.mean_volt = 0
                log_date = self.get_log_date(name)
                data = open(name, 'r')
                # print 'reading', name
                if ind == 0:
                    self.find_start(data, log_date)
                index = 0
                i = 0
                info = None
                for line in data:
                    info = line.split()
                    if i == 0:
                        # print info
                        pass
                    i += 1
                    if is_float(info[1]):
                        now = datetime.strptime(log_date + " " + info[0], "%Y-%m-%d %H:%M:%S")
                        index = self.averaging(dicts, now, info, index, key)
                        if self.stop != -1:
                            if self.stop < now:
                                stop[stop_ind] = True
                                break
                if info:
                    # print i,info, len(data.readlines())
                    # print
                    pass
                data.close()
                if stop[stop_ind]:
                    stop_ind += 1
                    break
        self.check_empty(dicts)
        ind += 1
        return dicts

    def averaging(self, dicts, now, info, index, key, shifting=False):
        # for key in self.keithleys:
        if len(info) > 2:
            # print self.start, now, self.stop
            if self.start < now < self.stop and float(info[2]) < 1e30:
                index += 1
                if self.do_averaging:
                    if not shifting:
                        self.mean_curr += float(info[2]) * 1e9
                        self.mean_volt += float(info[1])
                        if index % self.points == 0:
                            dicts[1][key].append(self.mean_curr / self.points)
                            dicts[0][key].append(convert_time(now))
                            dicts[2][key].append(self.mean_volt / self.points)
                            self.mean_curr = 0
                            self.mean_volt = 0
                    else:
                        if index <= self.points:
                            self.mean_curr += float(info[2]) * 1e9
                            dicts[1][key].append(self.mean_curr / index)
                            if index == self.points:
                                self.mean_curr /= self.points
                        else:
                            mean_curr = self.mean_curr * weight + (1 - weight) * float(info[2]) * 1e9
                            dicts[1][key].append(mean_curr)
                        dicts[0][key].append(convert_time(now))
                        dicts[2][key].append(float(info[1]))
                else:
                    dicts[1][key].append(float(info[2]) * 1e9)
                    dicts[0][key].append(convert_time(now))
                    dicts[2][key].append(float(info[1]))
        return index

    def find_start(self, data, log_date):
        lines = len(data.readlines())
        was_lines = 0
        data.seek(0)
        if lines > 10000:
            for i in range(6):
                lines /= 2
                for j in range(lines):
                    data.readline()
                while True:
                    info = data.readline().split()
                    if not info:
                        break
                    if is_float(info[1]):
                        now = datetime.strptime(log_date + " " + info[0], "%Y-%m-%d %H:%M:%S")
                        if now < self.start:
                            was_lines += lines
                            break
                        else:
                            data.seek(0)
                            for k in range(was_lines):
                                data.readline()
                            break

    def check_empty(self, dicts):
        for i in range(len(dicts)):
            for key in self.keithleys:
                if len(dicts[i][key]) == 0:
                    if i == 0:
                        if self.stop != -1:
                            dicts[i][key] = [convert_time(self.start), convert_time(self.stop)]
                        else:
                             dicts[i][key] = [convert_time(self.start), convert_time(time())]
                    else:
                        dicts[i][key] = [0]
                        dicts[i][key] = [0]

    def relative_time(self):
        for key in self.keithleys:
            zero = self.time_x[key][0]
            for i in range(len(self.time_x[key])):
                self.time_x[key][i] = self.time_x[key][i] - zero
        return self.time_x

    def make_dict(self, arg1, arg2):
        if self.single_mode:
            return OrderedDict([("Keithley1", arg1)])
        elif self.number == "2":
            return OrderedDict([("Keithley1", arg1),
                                ("Keithley2", arg2)])

if __name__ == '__main__':
    test = KeithleyInfo('logs_237', 'test.json', '2015-06-29.10:50', '2015-06-29.12:00', '1', False, 10)
