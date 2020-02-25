# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from analysis import load_parser, get_base_dir, OrderedDict, critical, Analysis, join
from json import load, loads


class DUT:
    """ Class with all information about a single DUT. """
    def __init__(self, number, run_info):

        self.Config = Analysis.load_main_config()
        self.Dir = get_base_dir()

        # Info
        self.Number = number
        self.Name = run_info['dia{}'.format(number)]
        self.Bias = run_info['dia{}hv'.format(number)]
        self.Attenuator = run_info['att_dia{}'.format(number)] if 'att_dia{}'.format(number) in run_info else None

        # Specs
        self.Specs = load_parser(join(self.Config.get('MAIN', 'data directory'), 'diaSpecs.ini'))
        self.Irradiation = self.load_irradiation()
        self.Thickness = self.load_spec('Thickness', typ=int)
        self.CCD = self.load_spec('CCD', typ=int)
        self.Size = loads(self.load_spec('Size')) if self.load_spec('Size') is not None else None

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

    def load_spec(self, section, typ=None):
        return None if not self.Specs.has_option(section, self.Name) else self.Specs.get(section, self.Name) if typ is None else typ(self.Specs.get(section, self.Name))

    def set_number(self, value):
        self.Number = value
