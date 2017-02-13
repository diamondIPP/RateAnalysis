#!/usr/bin/python
from ROOT import TFile, gROOT, TGraph, TH2F, gStyle, TCanvas, TCut, TH1F
from sys import argv
from AbstractClasses.Utils import *
from json import load
from AbstractClasses.RunClass import Run
from datetime import datetime
from ConfigParser import ConfigParser, NoSectionError
from numpy import mean
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar

widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]


def load_runinfo():
    info = {}
    try:
        fi = open(run.runinfofile, 'r')
        data = load(fi)
        fi.close()
    except NoSectionError as err:
        log_warning('{err}\nCould not load default RunInfo! --> Using default'.format(err=err))
        return None

    if run >= 0:
        info = data.get(str(run))
    return info


def draw_hitmap(show=True, cut=None, plane=None):
    planes = range(4) if plane is None else [int(plane)]
    cut = TCut('') if cut is None else TCut(cut)
    histos = [TH2F('h_hm{0}'.format(i_pl), 'Hitmap Plane {0}'.format(i_pl), 52, 0, 52, 80, 0, 80) for i_pl in planes]
    for plane in planes:
        cut_string = cut + TCut('plane == {0}'.format(plane))
        t.Draw('row:col>>h_hm{0}'.format(plane), cut_string, 'goff')
        run.format_histo(histos[plane], x_tit='col', y_tit='row')
    run.set_root_output(show)
    c = TCanvas('c_hm', 'Hitmaps', 2000, 2000)
    c.Divide(2, 2)
    for i, h in enumerate(histos, 1):
        h.SetStats(0)
        pad = c.cd(i)
        pad.SetBottomMargin(.15)
        h.Draw('colz')
    stuff.append([histos, c])
    run.set_root_output(True)


def trig_edges(nwf=None):
    nwf = entries if nwf is None else nwf
    pbar = ProgressBar(widgets=widgets, maxval=nwf).start()
    h = TH1F('h_te', 'Trigger Edges', 1024, 0, 1024)
    t.Draw('wf8', '', 'goff', nwf, 0)
    buf = t.GetV1()
    for i in xrange(nwf):
        pbar.update(i + 1)
        for k in xrange(1023):
            # print i * 1204 + j, int(buf[i * 1204 + j])
            if abs(buf[i * 1024 + k] - buf[i * 1024 + k + 1]) > 50:
                h.Fill(k)
    run.format_histo(h, x_tit='Bin Number', y_tit='Number of Entries', y_off=1.4, stats=0, fill_color=407)
    run.draw_histo(h, lm=.12)


def draw_waveforms(n=1000, start_event=0, cut_string='', show=True, fixed_range=None, ch=0):
    channel = ch
    global count
    start = start_event + count
    print 'starting at event', start
    if not wf_exists(channel):
        log_warning('This waveform is not stored in the tree!')
        return
    cut = cut_string
    n_events = find_n_events(n, cut, start)
    h = TH2F('wf', 'Waveform', 1024, 0, 1024, 1204, -512, 512)
    gStyle.SetPalette(55)
    t.Draw('wf{ch}:Iteration$>>wf'.format(ch=channel), cut, 'goff', n_events, start)
    h = TGraph(t.GetSelectedRows(), t.GetV2(), t.GetV1()) if n == 1 else h
    if fixed_range is None and n > 1:
        ymin, ymax = h.GetYaxis().GetBinCenter(h.FindFirstBinAbove(0, 2)), h.GetYaxis().GetBinCenter(h.FindLastBinAbove(0, 2))
        diff = ymax - ymin
        h.GetYaxis().SetRangeUser(ymin - diff * .1, ymax + diff * .1)
    elif fixed_range:
        assert type(fixed_range) is list, 'Range has to be a list!'
        h.GetYaxis().SetRangeUser(fixed_range[0], fixed_range[1])
    if show:
        gROOT.SetBatch(0)
    run.format_histo(h, title='Waveform', name='wf', x_tit='Time [ns]', y_tit='Signal [mV]', markersize=.4, y_off=.4, stats=0, tit_size=.05)
    run.draw_histo(h, '', show, lm=.06, rm=.045, draw_opt='scat' if n == 1 else 'col', x=1.5, y=.5)
    count += n_events
    return h, n_events


def show_single_waveforms(n=1, cut='', start_event=0):
    global count
    start = start_event + count
    activated_wfs = [wf for wf in xrange(4) if wf_exists(wf)]
    print 'activated wafeforms:', activated_wfs
    print 'Start at event number:', start
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
    evt_numbers = [t.GetV1()[i] for i in xrange(total_events)]
    return int(evt_numbers[:n][-1] + 1 - start)


def wf_exists(channel):
    wf_exist = True if t.FindBranch('wf{ch}'.format(ch=channel)) else False
    if not wf_exist:
        print 'The waveform for channel {ch} is not stored in the tree'.format(ch=channel)
    return wf_exist


def get_branch(n, name='IntegralPeaks[13]'):
    t.GetEntry(n)
    info = None
    exec 'info = t.{br}'.format(br=name)
    print info


def load_diamond_name(ch):
    parser = ConfigParser()
    parser.read('Configuration/DiamondAliases.cfg')
    dia = runinfo['dia{0}'.format(ch)]
    return parser.get('ALIASES', dia)


def calc_flux():
    f = open(joinpath(run.load_mask_file_dir(), runinfo['maskfile']), 'r')
    data = []
    for line in f:
        if len(line) > 3:
            line = line.split()
            data.append([int(line[2])] + [int(line[3])])
    f.close()
    pixel_size = 0.01 * 0.015
    area = [(data[1][0] - data[0][0]) * (data[1][1] - data[0][1]) * pixel_size, (data[3][0] - data[2][0]) * (data[3][1] - data[2][1]) * pixel_size]
    flux = [runinfo['for{0}'.format(i + 1)] / area[i] / 1000. for i in xrange(2)]
    return mean(flux)


def get_bias(ch):
    return runinfo['dia{0}hv'.format(ch)]


def draw_runinfo_legend(ch):
    testcamp = datetime.strptime(tc, '%Y%m')
    l = run.make_legend(.005, .156, y1=.003, x2=.39, nentries=2, felix=False, scale=1, margin=.05)
    l.SetTextSize(.05)
    l.AddEntry(0, 'Test Campaign: {tc}, Run {run} @ {rate:2.1f} kHz/cm^{{2}}'.format(tc=testcamp.strftime('%b %Y'), run=run, rate=run.calc_flux()), '')
    l.AddEntry(0, 'Diamond: {diamond} @ {bias:+}V'.format(diamond=load_diamond_name(ch), bias=get_bias(ch)), '')
    l.Draw()
    stuff.append(l)


def draw_peak_timings(xmin=0, xmax=100, ch=0):
    h = TH1F('h_pt', 'PeakTimings', 100, xmin, xmax)
    t.Draw('peak_positions{0}>>h_pt'.format(ch), '', 'goff')
    run.format_histo(h, x_tit='Digitiser Bin', y_tit='Number of Entries')
    run.draw_histo(h)


def draw_both_wf(ch, show=True):
    pulser = draw_waveforms(100, cut_string='pulser', show=False, ch=ch)[0]
    signal = draw_waveforms(100, cut_string='!pulser', show=False, ch=ch)[0]
    p_proj = pulser.ProjectionY()
    p = [p_proj.GetBinCenter(p_proj.FindFirstBinAbove(0)), p_proj.GetBinCenter(p_proj.FindLastBinAbove(0))]
    s_proj = signal.ProjectionY()
    s = [s_proj.GetBinCenter(p_proj.FindFirstBinAbove(0)), s_proj.GetBinCenter(p_proj.FindLastBinAbove(0))]
    y = [min(p + s) - 10, max(p + s) + 10]
    run.format_histo(pulser, y_range=y, stats=0, markersize=.4, y_off=.5, tit_size=.05)
    run.draw_histo(pulser, '', show, lm=.06, rm=.045, bm=.2, draw_opt='col', x=1.5, y=.5)
    signal.Draw('samecol')
    leg_ch = 1 if not ch else 2
    draw_runinfo_legend(leg_ch)
    run.save_plots('WaveForm_{run}_{ch}'.format(run=run, ch=ch, tc=tc), save_dir)


def read_macro(f):
    macro = f.Get('region_information')
    chan = None
    try:
        regions = [str(line) for line in macro.GetListOfLines()]
        chan = []
        for j, line in enumerate(regions):
            if 'Sensor Names' in line:
                chan = regions[j + 1].strip(' ').split(',')
        print chan
    except AttributeError:
        pass
    return chan


if __name__ == '__main__':
    rootfile = TFile(argv[1])
    tc = '20' + argv[1].split('/')[-1].strip('test').split('00')[0]
    tc += '0' if len(tc) < 6 else ''
    tc = '201610'

    try:
        run = int(argv[1].split('/')[-1].strip('.root').split('00')[-1])
    except (IndexError, ValueError):
        if argv[1].endswith('Tracks.root'):
            run = int(argv[1].split('/')[-1].strip('_withTracks.roottest'))
        elif 'Tracked' in argv[1]:
            run = int(argv[1].split('/')[-1].strip('.root').strip('TrackedRun'))
        else:
            run = None

    run = Run(run, test_campaign=tc, load_tree=False)

    channels = read_macro(rootfile)
    t = rootfile.Get('tree')
    entries = t.GetEntries()
    t.SetEstimate(entries * 1024)
    stuff = []
    count = 0
    bla = 5
    save_dir = 'WaveForms/'

    runinfo = load_runinfo()
    if len(argv) > 2:
        draw_both_wf(int(argv[2]), show=False)
