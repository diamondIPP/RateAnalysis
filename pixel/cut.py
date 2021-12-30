from helpers.utils import *
from src.cut import Cut, CutString, TCut
from src.binning import Bins
from plotting.draw import fit_fwhm


class PixCut(Cut):
    """ Generates all cut strings which correspond to a single DUT with pixel readout. """
    def __init__(self, analysis):
        Cut.__init__(self, analysis)

        self.N = self.DUT.Number + self.Run.NTelPlanes - 1
        self.generate_dut()
        self.ConsecutiveCuts = self.get_consecutive()

    # ----------------------------------------
    # region GENERATE
    def generate_dut(self):
        """ Generates the cut strings to apply in the analysis for each of the cuts. """
        self.CutStrings.register(self.generate_fiducial(center=True), level=23)
        self.CutStrings.register(self.generate_trigger_phase(), level=30)
        # self.CutStrings.register(self.generate_masks(), 61)
        self.CutStrings.register(self.generate_ncluster(), 80)
        self.CutStrings.register(self.generate_rhit(), 81)

    def generate_trigger_phase(self):
        v = self.load_dut_config('trigger phase')
        return CutString() if v is None else CutString('trigger_phase', f'trigger_phase[1] >= {v[0]} && trigger_phase[1]<={v[1]}', f'{v[0]} <= trigger phase <= {v[1]}')

    def generate_masks(self, col=None, row=None, pixels=None, exclude=True):
        cut_string = TCut('')
        cut_string += self.get_line_mask('col', col, exclude)
        cut_string += self.get_line_mask('row', row, exclude)
        cut_string += self.get_pixel_mask(pixels, exclude)
        return CutString('masks', cut_string, 'masking {} columns, {} rows and {} pixels'.format(*self.find_n_masked(col, row, pixels)) if exclude else '')

    def generate_rhit(self, value=None):
        (sx, sy), n = self.calc_res_sigmas(), choose(value, self.get_config('rhit sigma', dtype=float, required=True))
        cut = f'(sres_x[{self.N}] / {sx}) ** 2 + (sres_y[{self.N}] / {sy}) ** 2 <= {n ** 2}'
        return CutString('rhit', f'({cut} || sres[{self.N}] == -999)', f'(rx/{sx * 1e4:.0f})² + (ry/{sy * 1e4:.0f})² < {n}² [um] ({n} sigma)')

    def generate_ncluster(self, max_n=1):
        return CutString('ncluster', f'n_clusters[{self.N}] <= {max_n}', f'number of clusters <= {max_n}')
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def load_mask(self, option):
        data = self.Config.get_list('MASK', option)
        return data[self.DUT.Name] if self.DUT.Name in data else []

    def get_line_mask(self, var, lines=None, exclude=True):
        string = ' || '.join(f'cluster_{var}[{self.N}] == {line}' for line in (self.load_mask('column') if lines is None else make_list(lines)))
        return '' if lines is None and not exclude else Cut.invert(string) if exclude else string

    def get_pixel_mask(self, pixels=None, exclude=True):
        p = make_list(self.load_mask('pixel') if pixels is None else pixels)
        string = ' || '.join('cluster_col[{n}]=={} && cluster_row[{n}]=={}'.format(n=self.N, *tup) for tup in p.reshape(p.size // 2, 2))
        return '' if pixels is None and not exclude else Cut.invert(string) if exclude else string

    def get_plane(self):
        return TCut(f'plane == {self.N}')

    def get_nhit(self, n=1):
        return TCut(f'n_hits[{self.N}] == {n}')

    def get_ncluster(self, n=1, plane=None):
        return TCut(f'n_clusters[{choose(plane, self.N)}] == {n}')

    def get_cluster_size(self, n=1, plane=None):
        return TCut(f'cluster_size[{choose(plane, self.N)}] == {n}')

    def get_fid_lines(self, cols=None, rows=None, pix=None):
        fid = self.load_pixel_fid()
        cols, rows = make_list(cols) if cols is not None else arange(*fid[0][1:3] + [0, 1]), make_list(rows) if rows is not None else arange(*fid[1][:2] + [0, 1])
        return cols if pix is None else make_list(pix[0]), rows if pix is None else make_list(pix[1])

    def get_track_var(self, n, mode, mm=True, local=False, pixel=False):
        return super(PixCut, self).get_track_var(n, mode, mm, local, pixel)

    def get_track_vars(self, n, mm=True, local=False, pixel=False):
        return super(PixCut, self).get_track_vars(n, mm, local, pixel)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    @save_pickle('RHit', print_dur=True)
    def calc_res_sigmas(self, _redo=False):
        d = self.get_tree_vec([f'residuals_{m}[{self.N}]' for m in ['x', 'y']], self['tracks'] + self.get_ncluster(1, self.N))
        return array([fit_fwhm(self.Draw.distribution(i, show=False))[2].n for i in d])

    def find_n_masked(self, col, row, pixels):
        ncols, nrows = [make_list(self.load_mask(w) if lines is None else lines).size for w, lines in [('column', col), ('row', row)]]
        return ncols, nrows, make_list(self.load_mask('pixel') if pixels is None else pixels).size // 2

    @save_pickle('BeamStops', suf_args='all')
    def find_beam_interruptions(self, bin_width=10, threshold=.4):
        """ Looking for beam interruptions by evaluating the event rate. """
        bin_values, time_bins = histogram(self.Run.Time / 1000, bins=Bins(self.Ana).get_raw_time(bin_width)[1])
        m = mean(bin_values[bin_values.argsort()][-20:-10])  # take the mean of the 20th to the 10th highest bin to get an estimate of the plateau
        deviating_bins = where(abs(1 - bin_values / m) > threshold)[0]
        times = time_bins[deviating_bins] + bin_width / 2 - self.Run.Time[0] / 1000  # shift to the center of the bin
        not_connected = where(concatenate([[False], deviating_bins[:-1] != deviating_bins[1:] - 1]))[0]  # find the bins that are not consecutive
        times = split(times, not_connected)
        interruptions = [[self.get_event_at_time(v, rel=True) for v in [t[0], t[-1]]] for t in times] if len(times[0]) else []
        return interruptions
    # endregion COMPUTE
    # ----------------------------------------
