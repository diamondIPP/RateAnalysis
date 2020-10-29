# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.analysis import load_json, get_base_dir, OrderedDict, critical, load_main_config, join, ufloat, load, loads


class DUT:
    """ Class with all information about a single DUT. """
    def __init__(self, number, run_info):

        self.Config = load_main_config()
        self.Dir = get_base_dir()

        # Info
        self.Number = 1 if number is None else number
        self.Name = run_info['dia{}'.format(self.Number)]
        self.Bias = run_info['dia{}hv'.format(self.Number)]
        self.Attenuator = run_info['att_dia{}'.format(self.Number)] if 'att_dia{}'.format(self.Number) in run_info else None

        # Specs
        self.Specs = load_json(join(self.Config.get('MAIN', 'data directory'), 'dia_info.json'))[self.Name.upper()]
        self.Irradiation = self.load_spec('irradiation')
        self.Thickness = self.load_spec('thickness', typ=int, default=500)
        self.CCD = self.load_spec('CCD', typ=int)
        self.Size = self.load_spec('size', lst=True, default=[5, 5])
        self.PadSize = self.load_spec('metal', typ=float, error=.02, default=ufloat(3.5, .2))
        self.ActiveArea = self.PadSize ** 2 if self.PadSize is not None else None
        self.GuardRing = self.load_spec('guard ring', typ=float)

    def __str__(self):
        return 'DUT {}, {}, Bias: {:1.0f}V'.format(self.Number, self.Name, self.Bias)

    def __repr__(self):
        return self.__str__()

    def load_irradiation(self):
        with open(join(self.Dir, self.Config.get('MISC', 'irradiation file'))) as f:
            data = load(f)
            return OrderedDict([(key, dic[self.Name]) for key, dic in sorted(data.items()) if self.Name in dic])

    def get_irradiation(self, tc):
        return self.Irradiation[tc] if tc in self.Irradiation else critical('Please make an irradiation entry in the dia_info.json for "{}" in {}'.format(self.Name, tc))

    def load_spec(self, section, typ=None, lst=False, error=None, default=None):
        if section not in self.Specs or self.Specs[section] == 'None':
            return default
        spec = self.Specs[section] if typ is None else typ(self.Specs[section])
        return loads(spec) if lst and spec is not None else ufloat(spec, error) if error is not None and spec is not None else spec

    def set_number(self, value):
        self.Number = value
