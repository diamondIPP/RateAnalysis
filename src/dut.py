# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from analysis import load_json, get_base_dir, OrderedDict, critical, Analysis, join
from json import load, loads


class DUT:
    """ Class with all information about a single DUT. """
    def __init__(self, number, run_info):

        self.Config = Analysis.load_main_config()
        self.Dir = get_base_dir()

        # Info
        self.Number = 1 if number is None else number
        self.Name = run_info['dia{}'.format(self.Number)]
        self.Bias = run_info['dia{}hv'.format(self.Number)]
        self.Attenuator = run_info['att_dia{}'.format(self.Number)] if 'att_dia{}'.format(self.Number) in run_info else None

        # Specs
        self.Specs = load_json(join(self.Config.get('MAIN', 'data directory'), 'dia_info.json'))[self.Name]
        self.Irradiation = self.load_spec('irradiation')
        self.Thickness = self.load_spec('thickness', typ=int, default=500)
        self.CCD = self.load_spec('CCD', typ=int)
        self.Size = self.load_spec('size', lst=True)
        self.Metal = self.load_spec('metal', typ=float)
        self.GuardRing = self.load_spec('guard ring', typ=float)

    def __str__(self):
        return 'DUT {}, {}, Bias: {:1.0f}V'.format(self.Number, self.Name, self.Bias)

    def __repr__(self):
        return self.__str__()

    def load_irradiation(self):
        with open(join(self.Dir, self.Config.get('MISC', 'irradiation file'))) as f:
            data = load(f)
            return OrderedDict([(key, dic[self.Name]) for key, dic in sorted(data.iteritems()) if self.Name in dic])

    def get_irradiation(self, tc):
        return self.Irradiation[tc] if tc in self.Irradiation else critical('Please add "{}" to the irradiation file for {}'.format(self.Name, tc))

    def load_spec(self, section, typ=None, lst=False, default=None):
        spec = default if section not in self.Specs or self.Specs[section] == 'None' else self.Specs[section] if typ is None else typ(self.Specs[section])
        return loads(spec) if lst and spec is not None else spec

    def set_number(self, value):
        self.Number = value
