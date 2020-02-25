from utils import *
from cut import Cut, CutString, loads, invert, TCut, linspace


class CutPix(Cut):
    """ Generates all cut strings which correspond to a single DUT with pixel readout. """
    def __init__(self, analysis):
        Cut.__init__(self, analysis)

        self.DUTNumber = analysis.DUTNumber + 3
        self.DUTName = analysis.DUTName

        self.Bins = self.Analysis.Bins

        self.update_config()
        self.generate_dut()
        self.ConsecutiveCuts = self.generate_consecutive()

    @classmethod
    def from_parent(cls, parent):
        return cls(parent.Analysis)

    # ----------------------------------------
    # region CONFIG
    def update_config(self):
        """ Loads the pixel configuration parameters from the config file. """
        self.CutConfig['rhit'] = self.get_config('rhit', typ=int)
        self.CutConfig['trigger_phase'] = self.load_dut_config('trigger phase')
        self.CutConfig['row_mask'] = self.load_mask('row')
        self.CutConfig['col_mask'] = self.load_mask('column')
        self.CutConfig['pixel_mask'] = self.load_mask('pixel')
        self.CutConfig['fiducial'] = self.load_fiducial()
        self.CutConfig['local_fiducial'] = self.load_fiducial('pixel fiducial')

    def load_mask(self, name):
        data = loads(self.get_config(name, section='MASK'))
        return data[self.DUTName] if self.DUTName in data else []
    # endregion CONFIG
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE
    def generate_dut(self):
        """ Generates the cut strings to apply in the analysis for each of the cuts. """
        self.CutStrings.register(self.generate_trigger_phase(), level=22)
        self.CutStrings.register(self.generate_hit(), 60)
        self.CutStrings.register(self.generate_masks(), 61)
        self.CutStrings.register(self.generate_rhit(), 80)
        self.CutStrings.register(self.generate_fiducial(center=True), 90)

    def generate_aligned(self):
        description = '{:.1f}% of the events excluded'.format(100. * self.find_n_misaligned() / self.Analysis.Run.NEntries) if self.find_n_misaligned() else ''
        return CutString('aligned', 'aligned[{}]'.format(self.Analysis.DUTNumber) if self.find_n_misaligned() else '', description)

    def generate_trigger_phase(self):
        cut_range = self.CutConfig['trigger_phase']
        string = 'trigger_phase[1]>={min}&&trigger_phase[1]<={max}'.format(min=cut_range[0], max=cut_range[1]) if cut_range else ''
        return CutString('trigger_phase', string, '{} <= trigger phase <= {}'.format(*cut_range) if cut_range else '')

    @staticmethod
    def generate_hit():  # TODO implement cut! must change tree structure from tracking telescope (hits in the plane)
        """ Needs to be implemented. Have to change trackingTelescope for this """
        cut_string = ''
        return CutString('hit', cut_string, 'bla')

    def generate_masks(self, col=None, row=None, pixels=None, exclude=True):
        cut_string = TCut('')
        cut_string += self.generate_line_mask('col', col, exclude)
        cut_string += self.generate_line_mask('row', row, exclude)
        cut_string += self.generate_pixel_mask(pixels, exclude)
        return CutString('masks', cut_string, 'masking {} columns, {} rows and {} pixels'.format(*self.find_n_masked()))

    def generate_line_mask(self, var, lines=None, exclude=True):
        lines = self.CutConfig['{}_mask'.format(var)] if lines is None else make_list(lines)
        cut_var = 'cluster_{}[{}]'.format(var, self.DUTNumber)
        cut_string = ' || '.join(('{v}>={} && {v}<={}' if type(line) is list else '{v}=={}').format(v=cut_var, *make_list(line)) for line in lines)
        return '' if not cut_string else invert(cut_string) if exclude else cut_string

    def generate_pixel_mask(self, pixels=None, exclude=True):
        pixels = self.CutConfig['pixel_mask'] if pixels is None else pixels if type(pixels[0]) is list else [pixels]
        cut_string = ' || '.join('cluster_col[{n}]=={} && cluster_row[{n}]=={}'.format(n=self.DUTNumber, *tup) for tup in pixels)
        return '' if not cut_string else invert(cut_string) if exclude else cut_string

    def generate_rhit(self, value=None):
        value = self.CutConfig['rhit'] if value is None else value
        cut_value = self.compute_rhit(value)
        return CutString('rhit', 's_residuals[{d}] < {val}'.format(d=self.DUTNumber, val=cut_value), 'residual < {:1.1f}mm ({}% quantile)'.format(cut_value * 10, value))
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    def compute_rhit(self, value, redo=False):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'RHit', run=self.Analysis.RunNumber, ch=self.DUTNumber)

        def func():
            t = self.Analysis.info('generating rhit cut in for run {run}...'.format(run=self.Analysis.RunNumber), next_line=False)
            values = get_root_vec(self.Analysis.Tree, var='residuals[{}]'.format(self.DUTNumber))
            self.Analysis.add_to_info(t)
            return get_quantiles(values, linspace(0, .2, 2001))

        rhits = do_pickle(pickle_path, func, redo=redo)
        cut_value = rhits[value]
        return cut_value

    def find_n_misaligned(self):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'align', self.RunNumber)

        def f():
            return where(get_root_vec(self.Analysis.Tree, var='aligned[{}]'.format(self.DUTNumber + 3), dtype=bool) == 0)[0].size
        return do_pickle(pickle_path, f)

    def find_n_masked(self):
        cols, rows = [sum(v[1] - v[0] + 1 if type(v) is list else 1 for v in make_list(self.CutConfig['{}_mask'.format(var)])) for var in ['col', 'row']]
        return cols, rows, len(self.CutConfig['pixel_mask'])

    # endregion COMPUTE
    # ----------------------------------------
