#!/usr/bin/env python
from ROOT import TFile, gROOT, TGraph, TH2F, gStyle, TCanvas, TCut, TH1F
from sys import argv
from src.run import Run
from src.binning import Bins
from os import chdir, system
from helpers.draw import *
from numpy import genfromtxt, polyfit, polyval, quantile

widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
plot = Draw()


def load_runinfo():
    run_info = {}
    try:
        fi = open(run.InfoFile, 'r')
        data = load(fi)
        fi.close()
    except NoSectionError as err:
        warning('{err}\nCould not load default RunInfo! --> Using default'.format(err=err))
        return None

    if run >= 0:
        run_info = data.get(str(run))
    return run_info


def draw_hitmap(cut=None, plane=None):
    planes = range(4) if plane is None else [int(plane)]
    cut = TCut('') if cut is None else TCut(cut)
    histos = [TH2F('h_hm{0}'.format(i_pl), 'Hitmap Plane {0}'.format(i_pl), 52, 0, 52, 80, 0, 80) for i_pl in planes]
    for plane in planes:
        cut_string = cut + TCut('plane == {0}'.format(plane))
        t.Draw('row:col>>h_hm{0}'.format(plane), cut_string, 'goff')
        format_histo(histos[plane], x_tit='col', y_tit='row')
    c = TCanvas('c_hm', 'Hitmaps', 2000, 2000)
    c.Divide(2, 2)
    for i, h in enumerate(histos, 1):
        h.SetStats(0)
        pad = c.cd(i)
        pad.SetBottomMargin(.15)
        h.Draw('colz')
    stuff.append([histos, c])


def get_track_vars(roc=0, local=False):
    return ['cluster_{}pos_{}[{}]'.format(n, 'local' if local else 'tel', roc) for n in ['x', 'y']]


def draw_occupancies():
    c = plot.canvas('c', w=1.5, h=1.5, divide=(2, 2))
    for i in range(z.NPlanes):
        x, y = z.get_tree_vec(['cluster_col[{}]'.format(i), 'cluster_row[{}]'.format(i)])
        plot.histo_2d(x, y, Bins.get_pixel(), stats=0, show=0, canvas=c.cd(i + 1))


def trig_edges(nwf=None):
    nwf = entries if nwf is None else nwf
    pbar = ProgressBar(widgets=widgets, maxval=nwf).start()
    h = TH1F('h_te', 'Trigger Edges', 1024, 0, 1024)
    t.Draw('wf8', '', 'goff', nwf, 0)
    buf = t.GetV1()
    for i in range(nwf):
        pbar.update(i + 1)
        for k in range(1023):
            # print i * 1204 + j, int(buf[i * 1204 + j])
            if abs(buf[i * 1024 + k] - buf[i * 1024 + k + 1]) > 50:
                h.Fill(k)
    format_histo(h, x_tit='Bin Number', y_tit='Number of Entries', y_off=1.4, stats=0, fill_color=407)
    plot.histo(h, lm=.12)


def draw_waveforms(n=1000, start_event=0, cut_string='', show=True, fixed_range=None, ch=0):
    channel = ch
    global count
    start = start_event + count
    print('starting at event', start)
    if not wf_exists(channel):
        warning('This waveform is not stored in the tree!')
        return
    cut = cut_string
    n_events = find_n_events(n, cut, start)
    h = TH2F('wf', 'Waveform', 1024, 0, 1024, 1204, -512, 512)
    gStyle.SetPalette(55)
    t.Draw('wf{ch}:Iteration$>>wf'.format(ch=channel), cut, 'goff', n_events, start)
    h = TGraph(t.GetSelectedRows(), t.GetV2(), t.GetV1()) if n == 1 else h
    h.SetMarkerStyle(20)
    if fixed_range is None and n > 1:
        ymin, ymax = h.GetYaxis().GetBinCenter(h.FindFirstBinAbove(0, 2)), h.GetYaxis().GetBinCenter(h.FindLastBinAbove(0, 2))
        delta = ymax - ymin
        h.GetYaxis().SetRangeUser(ymin - delta * .1, ymax + delta * .1)
    elif fixed_range:
        assert type(fixed_range) is list, 'Range has to be a list!'
        h.GetYaxis().SetRangeUser(fixed_range[0], fixed_range[1])
    if show:
        gROOT.SetBatch(0)
    format_histo(h, title='Waveform', name='wf', x_tit='Time [ns]', y_tit='Signal [mV]', markersize=.4, y_off=.4, stats=0, tit_size=.05)
    plot.histo(h, show, lm=.06, rm=.045, draw_opt='alp' if n == 1 else 'col', w=1.5, h=.5)
    count += n_events
    return h, n_events


def show_single_waveforms(n=1, cut='', start_event=0):
    global count
    start = start_event + count
    activated_wfs = [wf for wf in range(4) if wf_exists(wf)]
    print('activated wafeforms:', activated_wfs)
    print('Start at event number:', start)
    wfs = [draw_waveforms(n=n, start_event=start, cut_string=cut, show=False, ch=wf) for wf in activated_wfs]
    n_wfs = len(activated_wfs)
    cs = {c.GetName(): c for c in gROOT.GetListOfCanvases()}
    if 'c_wfs' not in cs:
        c = TCanvas('c_wfs', 'Waveforms', 2000, n_wfs * 500)
        c.Divide(1, n_wfs)
    else:
        c = cs['c_wfs']
    for i, wf in enumerate(wfs, 1):
        wf[0].SetTitle('{nam} WaveForm'.format(nam=channels[activated_wfs[i - 1]]))
        wf[0].GetXaxis().SetTitleSize(.06)
        wf[0].GetYaxis().SetTitleSize(.06)
        wf[0].GetXaxis().SetLabelSize(.06)
        wf[0].GetYaxis().SetLabelSize(.06)
        pad = c.cd(i)
        pad.SetMargin(.05, .05, .1, .1)
        wf[0].Draw('alp')
    stuff.append([c, wfs])
    cnt = wfs[0][1]
    # if cnt is None:
    #     return
    count += cnt


def find_n_events(n, cut, start):
    total_events = t.Draw('event_number', cut, 'goff', entries, start)
    evt_numbers = [t.GetV1()[i] for i in range(total_events)]
    return int(evt_numbers[:n][-1] + 1 - start)


def wf_exists(channel):
    wf_exist = True if t.FindBranch('wf{ch}'.format(ch=channel)) else False
    if not wf_exist:
        print('The waveform for Channel {ch} is not stored in the tree'.format(ch=channel))
    return wf_exist


def get_branch(n, name='IntegralPeaks[13]'):
    print(get_tree_vec(t, name, nentries=1, firstentry=n))


def load_diamond_name(ch):
    p = ConfigParser()
    p.read('config/DiamondAliases.ini')
    dia = runinfo['dia{0}'.format(ch)]
    return p.get('ALIASES', dia)


def get_bias(ch):
    return runinfo['dia{0}hv'.format(ch)]


def draw_runinfo_legend(ch):
    testcamp = datetime.strptime(tc, '%Y%m')
    leg = Draw.make_legend(.005, .156, y1=.003, nentries=2, clean=False, scale=1, margin=.05)
    leg.SetTextSize(.05)
    leg.AddEntry(0, 'Test Campaign: {tc}, Run {run} @ {rate:2.1f} kHz/cm^{{2}}'.format(tc=testcamp.strftime('%b %Y'), run=run, rate=run.calculate_flux()), '')
    leg.AddEntry(0, 'Diamond: {diamond} @ {bias:+}V'.format(diamond=load_diamond_name(ch), bias=get_bias(ch)), '')
    leg.Draw()
    stuff.append(leg)


def draw_peak_timings(xmin=100, xmax=400, ch=0, corr=False):
    h = TH1F('h_pt', 'PeakTimings', xmax - xmin, xmin, xmax)
    t.Draw('max_peak_{p}[{c}]>>h_pt'.format(c=ch, p='position' if not corr else 'time'), '', 'goff')
    format_histo(h, x_tit='Digitiser Bin', y_tit='Number of Entries')
    plot.histo(h)


def draw_both_wf(ch, show=True):
    pulser = draw_waveforms(100, cut_string='pulser', show=False, ch=ch)[0]
    signal = draw_waveforms(100, cut_string='!pulser', show=False, ch=ch)[0]
    p_proj = pulser.ProjectionY()
    p = [p_proj.GetBinCenter(p_proj.FindFirstBinAbove(0)), p_proj.GetBinCenter(p_proj.FindLastBinAbove(0))]
    s_proj = signal.ProjectionY()
    s = [s_proj.GetBinCenter(p_proj.FindFirstBinAbove(0)), s_proj.GetBinCenter(p_proj.FindLastBinAbove(0))]
    y = [min(p + s) - 10, max(p + s) + 10]
    format_histo(pulser, y_range=y, stats=0, markersize=.4, y_off=.5, tit_size=.05)
    plot.histo(pulser, show, lm=.06, rm=.045, bm=.2, draw_opt='col', w=1.5, h=.5)
    signal.Draw('samecol')
    leg_ch = 1 if not ch else 2
    draw_runinfo_legend(leg_ch)


def read_macro(f):
    macro = f.Get('region_information')
    chan = None
    try:
        regions = [str(line) for line in macro.GetListOfLines()]
        chan = []
        for j, line in enumerate(regions):
            if 'Sensor Names' in line:
                chan = regions[j + 1].strip(' ').split(',')
        print(chan)
    except AttributeError:
        pass
    return chan


def draw_peak_values(ch=0, show=True):
    h = TH1F('h_pv', 'PH', 100, -450, -350)
    t.Draw('max_peak_position[{c}]>>h_pv'.format(c=ch), '', 'goff')
    format_histo(h, x_tit='Peak Value [mV]', y_tit='Number of Entries', y_off=1.4)
    plot.histo(h, show=show)
    return h


def fit_peak_values(ch=0, show=True):
    h = draw_peak_values(ch, show)
    format_statbox(h, fit=True)
    fit = h.Fit('gaus', 'sq')
    return fit


def draw_pedestal(ch=0, show=True):
    h = TH1F('h_pd', 'PH', 120, -30, 30)
    t.Draw('pedestals[{c}]>>h_pd'.format(c=ch), '', 'goff')
    format_histo(h, x_tit='Peak Value [mV]', y_tit='Number of Entries', y_off=1.4)
    plot.histo(h, show=show)
    return h


def fit_pedestal(ch=0, show=True):
    h = draw_pedestal(ch, show)
    format_statbox(h, fit=True)
    fit = h.Fit('gaus', 'sq')
    return fit


def draw_pulser_pulse_height(ch=0, show=True):
    h = TH1F('h_ph', 'PH', 2000, -500, 500)
    t.Draw('peak_values[{c}] - pedestals[{c}]>>h_ph'.format(c=ch), '', 'goff')
    format_histo(h, x_tit='Pulser Pulse Height [mV]', y_tit='Number of Entries', y_off=1.4, x_range=[h.GetMean() - 10, h.GetMean() + 20])
    plot.histo(h, show=show)
    return h


def fit_pulser_ph(ch=0, show=True):
    h = draw_pulser_pulse_height(ch, show)
    format_statbox(h, fit=True, entries=4, w=.3)
    fit = h.Fit('gaus', 'sq')
    print(fit.Parameter(1), fit.Parameter(2))
    return fit


def draw_integral(ch=0, show=True):
    h = TH1F('h_int', 'PH', 2000, -500, 500)
    t.Draw('peak_integrals[{c}] - pedestals[{c}]>>h_int'.format(c=ch), '', 'goff')
    format_histo(h, x_tit='Pulser Pulse Height [mV]', y_tit='Number of Entries', y_off=1.4, x_range=[h.GetMean() - 10, h.GetMean() + 20])
    plot.histo(h, show=show)
    return h


def fit_integral(ch=0, show=True):
    h = draw_integral(ch, show)
    format_statbox(h, fit=True, entries=4, w=.3)
    fit = h.Fit('gaus', 'sq')
    print(fit.Parameter(1), fit.Parameter(2))
    return fit


def calc_integral():
    n = t.Draw('peak_positions[0]:pedestals[0]', '', 'goff', 100)
    # ped = [t.GetV1()[i] for i in range(n)]
    pp = [t.GetV2()[i] for i in range(n)]
    n = t.Draw('wf0', '', 'goff', 100)
    wf = [t.GetV1()[i] for i in range(n)]
    wf = [wf[i * 1024: (i + 1) * 1024] for i in range(len(wf) - 1)]
    print(wf[0])
    integral = mean([sum(w[pp[i] - 20:pp[i] + 20]) for i, w in enumerate(wf)])
    print(integral)


def get_real_zero(nwf=0, channel=0):
    print(nwf)
    try:
        n = t.Draw('wf{ch}:Iteration$'.format(ch=channel), '', 'goff', 1, nwf)
        wf = [t.GetV1()[i] for i in range(n)]
        delta = [wf[i + 1] - wf[i] for i in range(len(wf) - 2)]
        m, s = mean_sigma(delta[5:200])
        start = next(delta.index(d) for d in delta[5:] if d > 6 * s)
        return wf[start]
    except StopIteration:
        return 0


def get_real_zeros(ch=0):
    z_ = [get_real_zero(i, ch) for i in range(entries)]
    # peaks = t.Draw('peak_values[{c}]:pedestals[{c}]'.format(c=ch), '', 'goff')
    m, s = mean_sigma([t.GetV1()[i] - t.GetV2()[i] - z_[i] for i in range(entries)])
    return m, s


def convert(file_name):
    if not file_exists(file_name):
        file_dir = dirname(file_name)
        chdir(file_dir)
        raw_file = file_name.replace('test', 'run').replace('root', 'raw')
        cmd = '~/scripts/converter.py -t waveformtree -p {r}'.format(r=raw_file)
        print(cmd)
        system(cmd)


def get_z_positions(e=0):
    x = genfromtxt(join(z.Converter.TrackingDir, 'data', 'alignments.txt'), usecols=[0, 2, 6])
    return array([ufloat(ix, e) if e else ix for ix in x[(x[:, 0] == z.Converter.TelescopeID) & (x[:, 1] > -1)][:, 2]])  # [cm]


def get_res_cut(res, n, cut, max_res):
    x = array([min(i) for i in split(res, cumsum(n)[:-1])])[cut]
    h = plot.distribution(x, make_bins(-1000, 1000, n=100), show=False)
    m = h.GetMean()
    return array((m - max_res < x) & (x < m + max_res))
    # f = Draw.make_f('f', 'gaus', -1000, 1000)
    # h.Fit(f, 'qs', '', *get_fwhm(h, ret_edges=True, err=False))
    # f.Draw('same')
    # return f.Integral(-1e5, 1e5) / h.GetBinWidth(1)
    # m, s = f.Parameter(1), f.Parameter(2)
    # return (x > m - nsig * s) & (x < m + nsig * s)


# noinspection PyTupleAssignmentBalance
def get_inner_efficiency(roc=1, q=.2):
    i_pl = 2 if roc == 1 else 1
    cut = TCut('n_clusters[0] == 1 & n_clusters[3] == 1 & n_clusters[{}] == 1'.format(i_pl))  # select only events with exactly one cluster in the other three planes
    n_tracks = t.GetEntries(cut.GetTitle())
    info('Found {} tracks'.format(n_tracks))
    cut += 'n_clusters[{}] > 0'.format(roc)
    z_ = get_z_positions()
    x0, y0 = z.get_tree_vec(get_track_vars(roc=0), cut=cut)
    xi, yi = z.get_tree_vec(get_track_vars(roc=i_pl), cut=cut)
    x3, y3 = z.get_tree_vec(get_track_vars(roc=3), cut=cut)
    (x, y), n = z.get_tree_vec(get_track_vars(roc=roc), cut=cut), z.get_tree_vec('n_clusters[{}]'.format(roc), cut, dtype='i2')
    xfits, chi_x, _, _, _ = polyfit(z_[[0, i_pl, 3]], [x0, xi, x3], deg=1, full=True)
    yfits, chi_y, _, _, _ = polyfit(z_[[0, i_pl, 3]], [y0, yi, y3], deg=1, full=True)
    chi_cut = (chi_x < quantile(chi_x, q)) & (chi_y < quantile(chi_y, q))
    dx, dy = array([polyval(xfits.repeat(n, axis=1), z_[roc]) - x, polyval(yfits.repeat(n, axis=1), z_[roc]) - y]) * 10000  # to um
    rcut = get_res_cut(dx, n, chi_cut, 2 * z.Plane.PX * 1000) & get_res_cut(dy, n, chi_cut, 2 * z.Plane.PY * 1000)
    # n_good = round(mean([get_res_cut(dx, n, chi_cut), get_res_cut(dy, n, chi_cut)]))
    n_max = chi_cut[chi_cut].size
    n_good = rcut[rcut].size
    eff = calc_eff(n_good, n_max)
    info('Efficiency for plane {}: ({:1.1f}+{:1.1f}-{:1.1f})%'.format(roc, *eff))
    return eff


def get_p_miss(q=.2):
    e1, e2 = (get_inner_efficiency(i, q) for i in [1, 2])
    p_miss = (100 - ufloat(e1[0], mean(e1[1:]))) * (100 - ufloat(e2[0], mean(e2[1:]))) / 100
    info('p-miss: {:1.3f}%'.format(p_miss))
    return p_miss


def get_raw_efficiency(roc=0):
    x = get_tree_vec(t, 'n_clusters[{}]'.format(roc))
    eff = calc_eff(values=x > 0)
    print('Efficiency of Plane {}: ({:1.1f}+{:1.1f}-{:1.1f})%'.format(roc, *eff))
    return eff


def draw_eff_evo(roc=0, n=50):
    e = get_tree_vec(t, 'n_clusters[{}]'.format(roc)) > 0
    plot.efficiency(arange(e.size), e, make_bins(e.size, n=n))


def get_efficiencies():
    return [get_raw_efficiency(i) for i in range(z.NPlanes)]


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('run')
    parser.add_argument('-tc', nargs='?', default=None)
    args = parser.parse_args()

    if isint(args.run):
        z = Run(int(args.run), args.tc)
        t = z.Tree
    else:
        rootfile = TFile(args.run)
        try:
            run = int(argv[1].split('/')[-1].strip('.root').split('00')[-1])
        except (IndexError, ValueError):
            if argv[1].endswith('Tracks.root'):
                run = int(argv[1].split('/')[-1].strip('_withTracks.roottest'))
            elif 'Tracked' in argv[1]:
                run = int(argv[1].split('/')[-1].strip('.root').strip('TrackedRun'))
                tc = remove_letters(args.run.split('/')[3]).replace('_', '')
            else:
                run = None

        try:
            run = Run(run, testcampaign=args.tc, load_tree=False)
            runinfo = load_runinfo()
        except ValueError:
            run = Run(2, testcampaign='201807', load_tree=False)

        channels = read_macro(rootfile)
        t = rootfile.Get('tree')
        entries = t.GetEntries()
        t.SetEstimate(entries * 1024)
        stuff = []
        count = 0
        save_dir = 'WaveForms/'

        if len(argv) > 2:
            draw_both_wf(int(argv[2]), show=False)
