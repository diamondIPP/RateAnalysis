# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.utils import load_json, get_base_dir, OrderedDict, critical, load_main_config, join, ufloat, load, loads, sqrt, pi, array, choose, add_err, prep_kw, update_pbar, PBar, do_nothing
from numpy import multiply, deg2rad, tan, rad2deg, arange, cos, sin, arctan, linspace, mean, invert, isnan
from plotting.draw import Draw, bins_from_vec


class DUT:
    """ Class with all information about a single DUT. """

    K = add_err(ufloat(2.65, .13), .18) * 1e-18  # damage constant for neutrons in pCVD diamond [cm²/um]

    def __init__(self, number, run_info):

        self.Config = load_main_config()
        self.Dir = get_base_dir()

        # Info
        self.Number = 1 if number is None else number
        self.Name = run_info['dia{}'.format(self.Number)]
        self.Key = self.Name.split('_')[0]
        self.Bias = float(run_info['dia{}hv'.format(self.Number)])
        self.Attenuator = run_info['att_dia{}'.format(self.Number)] if 'att_dia{}'.format(self.Number) in run_info else None

        # Specs
        self.Specs = load_json(join(self.Dir, 'Runinfos', 'dia_info.json'))[self.Key]
        self.Type = self.load_spec('type', default='pad')
        self.Irradiation = self.load_spec('irradiation')
        self.Version = self.load_spec('version')
        self.Thickness = self.load_spec('thickness', typ=int, default=500)
        self.CCD = self.load_spec('CCD', typ=int)
        self.Size = self.load_spec('size', lst=True, default=[5, 5])
        self.PadSize = self.load_spec('metal', typ=float, error=.02, default=ufloat(3.5, .2))
        self.ActiveArea = self.PadSize ** 2 if self.PadSize is not None else None
        self.GuardRing = self.load_spec('guard ring', typ=float)

    def __str__(self):
        return self.Name

    def __repr__(self):
        return f'{self.__class__.__name__} {self.Number}, {self.Name}, Bias: {self.Bias:1.0f}V'

    def full_name(self, tc):
        return str(self) if self.Version is None else f'{self}-V{next((i for i, v in enumerate(self.Version) if v > tc), len(self.Version))}'

    def load_irradiation(self):
        with open(join(self.Dir, self.Config.get('MISC', 'irradiation file'))) as f:
            data = load(f)
            return OrderedDict([(key, dic[self.Key]) for key, dic in sorted(data.items()) if self.Key in dic])

    def get_irradiation(self, tc):
        return self.Irradiation[tc] if tc in self.Irradiation else critical(f'Please make an irradiation entry in the dia_info.json for "{self.Key}" in {tc}')

    def load_spec(self, section, typ=None, lst=False, error=None, default=None):
        if section not in self.Specs or self.Specs[section] == 'None':
            return default
        spec = self.Specs[section] if typ is None else typ(self.Specs[section])
        return loads(spec) if lst and spec is not None else ufloat(spec, error) if error is not None and spec is not None else spec

    def set_number(self, value):
        self.Number = value

    def get_area(self, bcm=False):
        """ :returns: area of the DUT in cm²"""
        return self.get_bcm_area() if bcm else self.ActiveArea * .01

    def get_bcm_area(self):
        """ :returns: total area of the BCM' pad sizes """
        i = int(self.Name.split('-')[-1]) - 1
        base_length = 0.0928125  # [cm]
        spacing = 0.0025
        radius = 0.0049568
        rounded_edge = radius ** 2 * (4 - pi)
        return 2 ** i * base_length ** 2 + get_spacings(i, spacing, base_length) - rounded_edge

    def get_e_field(self, bias):
        return choose(bias, self.Bias) / self.max_drift_distance

    @property
    def max_drift_distance(self):
        return self.Thickness

    @property
    def max_fluence(self):
        return 1 / (self.K * self.max_drift_distance)


def get_spacings(i, spacing, length):
    """ :returns: the additional spacings for the BCM' pad sizes """
    if not i:
        return 0
    j = 2 ** (i / 2)
    return spacing * (j * length + (j - 1) * spacing) + 2 * get_spacings(i - 1, spacing, length)


class PixelDUT(DUT):

    def __init__(self, number, run_info):
        super().__init__(number, run_info)

        self.PX, self.PY = self.load_spec('cell size', lst=True, default=[Plane.PX, Plane.PY])
        self.A = self.PX * self.PY
        self.GX, self.GY = self.PX / Plane.PX / 1e3, self.PY / Plane.PY / 1e3  # ganging
        self.ActivePixels = self.load_spec('active pixels', lst=True, default=[Plane.NCols, Plane.NRows])
        self.W, self.H = self.PX / self.GX * self.ActivePixels[0] / 1e3, self.PY / self.GY * self.ActivePixels[1] / 1e3
        self.ColDia = self.load_spec('column diameter', typ=float)
        self.Is3D = self.ColDia is not None
        self.ColArea = (self.ColDia / 2) ** 2 * pi if self.Is3D else None
        self.ColRatio = 2 * self.ColArea / self.A if self.Is3D else None  # two columns per cell

        self.Draw = Draw()
        self.PBar = PBar()

    def get_area(self, bcm=False):
        return Plane.PixArea * multiply(*self.Size) * .01

    @property
    def max_drift_distance(self):
        return sqrt(self.PX ** 2 + self.PY ** 2) / 2 if self.Is3D else self.Thickness

    @property
    def r_col2area(self):
        return f'{self.ColRatio * 100:.2f} %'
    
    def crit_angle(self, n=1, x=None):
        return arctan(n * choose(x, self.PX) / self.Thickness)

    def n_cells(self, a, x=None):
        return tan(a) * self.Thickness / choose(x, self.PX)

    def draw_n_cells(self, amax=30, x=None, **dkw):
        f = Draw.make_tf1('n cells', lambda a: self.n_cells(deg2rad(a), x), 0, amax, npx=200)
        return self.Draw(f, **prep_kw(dkw, x_tit='Rotation Angle', y_tit='Number of Intersected Cells', y_range=[0, f.GetMaximum() * 1.1]))

    def draw_n_cellsx(self, amax=30, x=None, **dkw):
        x = choose(x, [50, 150])
        f = [self.draw_n_cells(amax, ix, show=False) for ix in x]
        self.Draw.functions(f, [f'{ix} #mum' for ix in x], **prep_kw(dkw, grid=True, left=True))

    def draw_crit_angle(self, n=1, **dkw):
        f = Draw.make_tf1('crit angle', lambda x: rad2deg(self.crit_angle(n, x)), 0, 160, npx=160)
        return self.Draw(f, **prep_kw(dkw, x_tit='Cell Size', y_tit='Rotation Angle', grid=True, ndivx=503))

    def draw_crit_angles(self, n=4, **dkw):
        f = [self.draw_crit_angle(i, show=False) for i in arange(n) + 1]
        self.Draw.functions(f, [f'n pixels = {i + 1}' for i in arange(n)], **prep_kw(dkw, ndivx=503, grid=True, left=True))

    def path_per_cell(self, a, x=None):
        a = deg2rad(a)
        return self.Thickness / cos(a) if a < self.crit_angle(x=x) else (choose(x, self.PX) - self.ColDia) / sin(a)

    def draw_path_per_cell(self, amax=30, x=None, **dkw):
        f = Draw.make_tf1('cellpath', self.path_per_cell, 0, amax, x=x, npx=200)
        return self.Draw(f, **prep_kw(dkw, x_tit='Rotation Angle', y_tit='Max Path Length Per Cell [#mum]', y_range=[0, f.GetMaximum() * 1.1]))

    def draw_path_per_x(self, amax=30, x=None, **dkw):
        x = choose(x, [50, 150])
        f = [self.draw_path_per_cell(amax, ix, show=False) for ix in x]
        self.Draw.functions(f, [f'{ix} #mum' for ix in x], **prep_kw(dkw, grid=True, y_range=[0, f[1].GetMaximum() * 1.1]))

    def min_path(self, a):
        a, t, d, s = deg2rad(a), self.Thickness, self.ColDia, self.PX
        if not a:
            return 0
        l, l0 = t / cos(a), d / sin(a)
        if a < self.crit_angle(1):
            return max(0, l - l0) / t
        if a < self.crit_angle(1, s + d):
            l1 = (t * tan(a) - s) / sin(a)
            return (l - l0 - l1) / t
        if a < arctan(2 * s / t):
            return l - 2 * l0
        if a < arctan((2 * s + d) / t):
            l2 = (t * tan(a) - 2 * s) / sin(a)
            return l - 2 * l0 - l2
        if a < arctan(3 * s / t):
            return l - 3 * l0
        return 1

    def draw_min_path(self, max_a=5, **dkw):
        f = Draw.make_tf1('path', self.min_path, 0, max_a, npx=200)
        self.Draw(f, **prep_kw(dkw, x_tit='Rotation Angle', y_tit='Relative Minimum Path Length', file_name='MinPath', y_range=[0, 1.1], gridy=True))

    # ----------------------------------------
    # region PATH
    @update_pbar
    def path_length(self, a=0, x=0, y=0, cols=1, _no_update=False):
        """ returns path lenght of a track through 3D detectors. coordinate centre is in centre of 3D column"""
        a, t, d, s = deg2rad(a), self.Thickness, self.ColDia, self.PX
        lt = t / cos(a)                                                 # total path length
        d *= sqrt(cols)                                                 # account for several cols by increasing area of a single one
        f = sqrt(1 - (2 * y / d) ** 2) if abs(y) < d / 2 else 0         # y dependence (chord of the circle / d)
        if f == 0:
            return lt
        d *= f
        x = (x - d / 2) % s                                             # move col from centre to right edge
        ld = min(lt, d / sin(a) if a else lt)                           # max path through one column
        n = int(self.n_cells(a))                                        # min number of intersected columns
        xd = t * tan(a) % s                                             # distance in x travelled through the detector
        xc = s - xd - d                                                 # x pos when next col is hit
        max_path = lt - n * ld
        if x < xc:
            return max_path
        xf = min(xd, d)                                                 # max distance in x travelled through column
        if x < xc + xf:
            return max_path - (x - xc) / xf * ld
        if x < s - xf:
            return max_path - ld
        return max_path - (1 - (x - s + xf) / xf) * ld

    def px(self, a=0, y=0, n=1000):
        return mean([self.path_length(a, x, y) for x in linspace(-self.PX / 2, self.PX / 2, n)])

    def p_mean(self, a=0, n=1000, pbar=True):
        x, y = [linspace(-i / 2, i / 2, n + 1) for i in [self.PX, self.PY]]
        self.PBar.start(x.size ** 2) if pbar else do_nothing()
        d = array([self.path_length(a, ix, iy, cols=2) for ix in x for iy in y])
        return mean(d[invert(isnan(d))])

    def draw_pa(self, max_a=30, x=0, y=0, cols=1, **dkw):
        f = Draw.make_tf1('xpath', self.path_length, 0, max_a, npx=200, cols=cols, y=y, x=x)
        return self.Draw(f, **prep_kw(dkw, x_tit='Rotation Angle', y_tit='Path Length [#mum]', y_range=[0, f.GetMaximum() * 1.1]))

    def draw_px(self, a=0, y=0, cols=1, xmax=None, **dkw):
        xmax = choose(xmax, self.PX)
        f = Draw.make_tf1('xpath', lambda x: self.path_length(a=a, x=x, y=y, cols=cols), -xmax / 2, xmax / 2, npx=500)
        return self.Draw(f, **prep_kw(dkw, x_tit='X [#mum]', y_tit='Path Length [#mum]', y_range=[0, f.GetMaximum() * 1.1]))

    def draw_pxy(self, a=0, n=100, rx=None, ry=None, **dkw):
        x, y = [linspace(-i / 2, i / 2, n + 1) for i in [choose(rx, self.PX), choose(ry, self.PY)]]
        self.PBar.start(x.size ** 2)
        d = array([[ix, iy, self.path_length(a, ix, iy)] for ix in x for iy in y]).T
        bins = sum([bins_from_vec(i, centre=True) for i in [x, y]], [])
        return self.Draw.prof2d(*d, bins, f'Eff at {a:.1f} deg', **prep_kw(dkw, x_tit='X [#mum]', y_tit='Y [#mum]', z_tit='Path Length [#mum]', stats=False, z_range=[0, max(d[2])]))

    def draw_pxa(self, y=0, amax=30, n=20, **dkw):
        x = linspace(0, amax, n)
        self.Draw.graph(x, [self.px(a, y) for a in x], **prep_kw(dkw, x_tit='Rotation Angle [deg]', y_tit='Average Path Length', draw_opt='al', lw=2, color=2))

    def draw_path_length(self, amax=10, n=10, ns=200, **dkw):
        sx, sy = [linspace(-i / 2, i / 2, ns + 1) for i in [self.PX, self.PY]]  # noqa
        self.PBar.start(sx.size ** 2 * (n + 1))
        x = linspace(0, amax, n + 1)
        y = array([self.p_mean(a, ns, pbar=False) for a in x]) / self.Thickness
        return self.Draw.graph(x, y, **prep_kw(dkw, x_tit='Rotation Angle [deg]', y_tit='Normalised Average Path Length', draw_opt='al', lw=2, color=2, lm=.14, y_off=1.5))

    def draw_path_lengths(self, ana, amax=10, n=10, ns=200, **dkw):
        g = [f(amax, n, ns, show=False) for f in [self.draw_path_length, ana.DUT.draw_path_length]]  # noqa
        self.Draw.multigraph(g, 'b', [str(self), str(ana.DUT)], **prep_kw(dkw, lm=.14, draw_opt='l'))
    # endregion PATH
    # ----------------------------------------

    # ----------------------------------------
    # region EFFICIENCY
    @update_pbar
    def eff(self, thresh=.2, a=0, x=0, y=0, cols=1):
        return self.path_length(a, x, y, cols, _no_update=True) > thresh * self.Thickness if abs(y) < self.ColDia / 2 else True

    def ex(self, thresh=.2, a=0, y=0, cols=1, n=10000):
        return mean([self.eff(thresh, a, x, y, cols) for x in linspace(-self.PX / 2, self.PX / 2, n)])

    def mean_eff(self, thresh=.2, a=0, n=1000, cols=2, pbar=True):
        x, y = [linspace(-i / 2, i / 2, n + 1) for i in [self.PX, self.PY]]
        self.PBar.start(x.size ** 2) if pbar else do_nothing()
        e = mean([self.eff(thresh, a, ix, iy) for ix in x for iy in y])
        return 1 - (1 - e) * cols

    def draw_ea(self, thresh=.2, x=0, y=0, amax=30, **dkw):
        f = Draw.make_tf1('xpath', lambda a: 100 * self.eff(thresh, a, x, y), 0, amax, npx=500)
        return self.Draw(f, **prep_kw(dkw, x_tit='Rotation Angle [deg]', y_tit='Efficiency [%]', y_range=[0, 105]))

    def draw_ex(self, thresh=.2, a=0, y=0, xmax=None, **dkw):
        xmax = choose(xmax, self.PX)
        f = Draw.make_tf1('xpath', lambda x: 100 * self.eff(thresh, a, x, y), -xmax / 2, xmax / 2, npx=500)
        return self.Draw(f, **prep_kw(dkw, x_tit='X [#mum]', y_tit='Efficiency [%]', y_range=[0, 105]))

    def draw_exy(self, thresh=.2, a=0, n=100, bw=1, **dkw):
        x, y = [linspace(-i / 2, i / 2, n + 1) for i in [self.PX, self.PY]]
        self.PBar.start(x.size ** 2)
        d = array([[ix, iy, self.eff(thresh, a, ix, iy)] for ix in x for iy in y]).T
        bins = [bins_from_vec(linspace(-i / 2, i / 2, n // bw + 1), centre=True) for i in [self.PX, self.PY]]
        return self.Draw.prof2d(*d, bins[0] + bins[1], f'Eff at {a:.1f} deg', **prep_kw(dkw, x_tit='X [#mum]', y_tit='Y [#mum]', z_tit='Efficiency', stats=False))

    def draw_eax(self, thresh=.2, y=0, cols=1, amax=1, n=20, **dkw):
        x = linspace(0, amax, n)
        self.Draw.graph(x, [self.ex(thresh, a, y, cols) for a in x], **prep_kw(dkw, draw_opt='al', lw=2, color=2))

    def draw_eff(self, thresh=.2, amax=1, n=10, ns=200, pbar=True, **dkw):
        sx, sy = [linspace(-i / 2, i / 2, ns + 1) for i in [self.PX, self.PY]]  # noqa
        self.PBar.start(sx.size ** 2 * (n + 1)) if pbar else do_nothing()
        x = linspace(0, amax, n + 1)
        y = array([self.mean_eff(thresh, a, ns, pbar=False) for a in x]) * 100
        return self.Draw.graph(x, y, **prep_kw(dkw, x_tit='Rotation Angle [deg]', y_tit='Efficiency [%]', draw_opt='al', lw=2, color=2))

    def draw_effs(self, tmax=.6, amax=1, n=10, np=10, ns=200, **dkw):
        self.PBar.start((ns + 1) ** 2 * n * (np + 1))
        x = linspace(tmax / n, tmax, n)
        g = [self.draw_eff(t, amax, np, ns, pbar=False, show=False) for t in x]
        self.Draw.multigraph(g, 'EffThresh', [f'{t:.2f}' for t in x], **prep_kw(dkw, draw_opt='l', wleg=.1))
    # endregion EFFICIENCY
    # ----------------------------------------


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
    AxTits = {'x_tit': 'Column', 'y_tit': 'Row'}

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
    def get_mask_dim(v=None, mm=True):
        return array([Plane.WX, Plane.WY] if v is None else [Plane.PX * (v[2] - v[0] + 1), Plane.PY * (v[3] - v[1] + 1)]) / (1 if mm else 10)

    @staticmethod
    def get_xpix(aspect_ratio=False):
        return round((max(Plane.WX, Plane.WY) - Plane.WX) / 2 / Plane.PX) if aspect_ratio else 0

    @staticmethod
    def get_ypix(aspect_ratio=False):
        return round((max(Plane.WX, Plane.WY) - Plane.WY) / 2 / Plane.PY) if aspect_ratio else 0
