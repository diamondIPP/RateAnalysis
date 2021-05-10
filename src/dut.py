# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.utils import load_json, get_base_dir, OrderedDict, critical, load_main_config, join, ufloat, load, loads, sqrt, pi


class DUT:
    """ Class with all information about a single DUT. """
    def __init__(self, number, run_info):

        self.Config = load_main_config()
        self.Dir = get_base_dir()

        # Info
        self.Number = 1 if number is None else number
        self.Name = run_info['dia{}'.format(self.Number)]
        self.Bias = float(run_info['dia{}hv'.format(self.Number)])
        self.Attenuator = run_info['att_dia{}'.format(self.Number)] if 'att_dia{}'.format(self.Number) in run_info else None

        # Specs
        self.Specs = load_json(join(self.Config.get('Directories', 'data'), 'dia_info.json'))[self.Name.upper()]
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

    def get_area(self, bcm=False):
        """ :returns: area of the DUT in cm^2"""
        return self.get_bcm_area() if bcm else self.ActiveArea * .01

    def get_bcm_area(self):
        """ :returns: total area of the BCM' pad sizes """
        i = int(self.Name.split('-')[-1]) - 1
        base_length = 0.0928125  # [cm]
        spacing = 0.0025
        radius = 0.0049568
        rounded_edge = radius ** 2 * (4 - pi)
        return 2 ** i * base_length ** 2 + get_spacings(i, spacing, base_length) - rounded_edge


def get_spacings(i, spacing, length):
    """ :returns: the additional spacings for the BCM' pad sizes """
    if not i:
        return 0
    j = 2 ** (i / 2)
    return spacing * (j * length + (j - 1) * spacing) + 2 * get_spacings(i - 1, spacing, length)


class Plane(object):
    """ Class with all information about a single pixel plane. """

    conf = load_main_config()
    Type = conf.get('PLANE', 'name')
    NCols, NRows = conf.get_list('PLANE', 'pixel')
    NPixels = NCols * NRows
    PX, PY = conf.get_list('PLANE', 'pitch')  # mm
    PixArea = PX * PY
    Area = PixArea * NPixels
    WX, WY = PX * NCols, PY * NRows
    WMax = max(WX, WY)
    R0 = 2 / sqrt(12)
    Frequency = conf.getfloat('PLANE', 'clock frequency')

    def __str__(self):
        return '{} Plane'.format(Plane.Type.title())

    def __repr__(self):
        return '{} Plane with {}x{} pixels of a size {:1.1f} x {:1.1f}um'.format(Plane.Type.upper(), Plane.NCols, Plane.NRows, Plane.PX * 1e3, Plane.PY * 1e3)

    def __call__(self, number=None):
        return self

    @staticmethod
    def get_area(v=None):
        return Plane.Area / 100 if v is None else Plane.PixArea / 100 * (v[2] - v[0] + 1) * (v[3] - v[1] + 1)

    @staticmethod
    def get_xpix(aspect_ratio=False):
        return round((max(Plane.WX, Plane.WY) - Plane.WX) / 2 / Plane.PX) if aspect_ratio else 0

    @staticmethod
    def get_ypix(aspect_ratio=False):
        return round((max(Plane.WX, Plane.WY) - Plane.WY) / 2 / Plane.PY) if aspect_ratio else 0
