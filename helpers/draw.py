#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TCanvas, gStyle, TLegend, TArrow, TPad, TCutG, TLine, TPaveText, TPaveStats, TH1F, TEllipse, TColor, TProfile
from ROOT import TProfile2D, TH2F, THStack, TMultiGraph, TPie, gROOT
from numpy import sign, linspace, ones, ceil, append, tile, absolute, rot90, flip, argsort
from helpers.utils import *


def get_color_gradient():
    stops = array([0., .5, 1], 'd')
    green = array([0. / 255., 200. / 255., 80. / 255.], 'd')
    blue = array([0. / 255., 0. / 255., 0. / 255.], 'd')
    red = array([180. / 255., 200. / 255., 0. / 255.], 'd')
    color_gradient = TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 255)
    return array([color_gradient + ij for ij in range(255)])


def load_resolution(default=800):
    try:
        from screeninfo import get_monitors
        return int(round_up_to(get_monitors()[0].height // 2, 100))
    except Exception as err:
        warning(err)
        return default


class Draw(object):

    Verbose = False
    Config = None
    Count = {}
    Res = load_resolution()
    Colors = get_color_gradient()
    Objects = []
    Dir = get_base_dir()

    Title = True
    Legend = False
    FillColor = 871
    Font = 42

    DefaultStats = {'x2': None, 'y2': None, 'h': None, 'w': .3, 'entries': False, 'm': False, 'rms': False, 'all_stat': True, 'fit': False, 'center_x': False, 'center_y': False, 'form': None}
    Stats = {}

    def __init__(self, config='', verbose=True):

        # Basics
        Draw.Verbose = verbose
        Draw.Config = Config(choose(config, default=join(Draw.Dir, 'config', 'main.ini')))

        # Settings
        self.IColor = 0  # color index
        Draw.Title = Draw.Config.get_value('SAVE', 'activate title', default=True)
        Draw.Legend = Draw.Config.get_value('SAVE', 'info legend', default=False)
        Draw.FillColor = Draw.Config.get_value('PLOTS', 'fill color', default=821)
        Draw.Font = Draw.Config.get_value('PLOTS', 'legend font', default=42)

        Draw.setup()

        self.Dic = {'TH1F': self.distribution, 'TH1I': self.distribution, 'TH1D': self.distribution,
                    'TH1': self.function,
                    'TGraph': self.graph, 'TGraphErrors': self.graph, 'TGraphAsymmErrors': self.graph,
                    'TProfile': self.profile,
                    'TH2I': self.histo_2d, 'TH2D': self.histo_2d, 'TH2F': self.histo_2d,
                    'TProfile2D': self.prof2d,
                    'TMultiGraph': self.multigraph}

    def __call__(self, th, *args, **kwargs):
        if th.ClassName() in self.Dic:
            return self.Dic[th.ClassName()](th, *args, **kwargs)
        return Draw.histo(th, *args, **kwargs)

    @staticmethod
    def add(*args):
        for obj in args:
            Draw.Objects.append(obj)
        Draw.clean()
        return args[0] if len(args) == 1 else args

    @staticmethod
    def clean():
        for obj in Draw.Objects:
            # if '0x(nil)' in str(obj) or obj is None:
            if obj is None:
                Draw.Objects.remove(obj)

    # ----------------------------------------
    # region SET
    @staticmethod
    def setup():
        gStyle.SetLegendFont(Draw.Font)
        gStyle.SetOptTitle(Draw.Title)
        gStyle.SetPalette(Draw.Config.get_value('PLOTS', 'palette', default=1))
        gStyle.SetNumberContours(Draw.Config.get_value('PLOTS', 'contours', default=20))

    @staticmethod
    def set_margin(c, side, value=None, default=.1, off=0):
        do(getattr(c, f'Set{side}Margin'), None if round(getattr(c, f'Get{side}Margin')(), 2) != .1 and value is None else max(choose(value, default) + off, 0))

    @staticmethod
    def set_pad_margins(c=None, l_=None, r=None, b=None, t=None):
        Draw.set_margin(c, 'Left', l_, default=.13)
        Draw.set_margin(c, 'Right', r, default=.02)
        Draw.set_margin(c, 'Bottom', b, default=.115, off=.06 if Draw.Legend else 0)
        Draw.set_margin(c, 'Top', t, default=.02, off=.08 if Draw.Title else 0)
    # endregion SET
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    @staticmethod
    def color(n, i):
        return Draw.get_colors(n)[i]

    def get_color(self, n, i=None):
        color = Draw.get_colors(n)[choose(i, self.IColor)]
        if i is None:
            self.IColor = self.IColor + 1 if self.IColor < n - 1 else 0
        return color

    @staticmethod
    def get_colors(n):
        return Draw.Colors[linspace(0, Draw.Colors.size - 1, n).round().astype(int)].tolist()

    @staticmethod
    def get_count(name='a'):
        if name not in Draw.Count:
            Draw.Count[name] = -1
        Draw.Count[name] += 1
        return Draw.Count[name]

    @staticmethod
    def reset_count(name='a'):
        Draw.Count[name] = 0
    
    @staticmethod
    def get_name(string='a'):
        return '{}{}'.format(string, Draw.get_count(string))
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAWING
    @staticmethod
    def canvas(title='c', x=None, y=None, w=1., h=1., logx=None, logy=None, logz=None, gridx=None, gridy=None, transp=None, divide=None, show=True):
        set_root_output(show)
        c0 = get_last_canvas(warn=False)
        x = x if x is not None else 0 if c0 is None else c0.GetWindowTopX() + 50
        y = y if y is not None else 0 if c0 is None else c0.GetWindowTopY() + 20
        c = TCanvas(Draw.get_name('c'), title, int(x), int(y), int(w * Draw.Res), int(h * Draw.Res))
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx, gridy])
        do(make_transparent, c, transp)
        if divide is not None:
            c.Divide(*(divide if type(divide) in [list, tuple] else [divide]))
        return Draw.add(c)

    @staticmethod
    def axis(x1, x2, y1, y2, title, limits=None, name='ax', col=1, width=1, off=.15, tit_size=.035, lab_size=0.035, tick_size=0.03, line=False, opt='+SU', l_off=.01, log=False, center=None):
        limits = ([y1, y2] if x1 == x2 else [x1, x2]) if limits is None else limits
        a = TGaxis(x1, y1, x2, y2, limits[0], limits[1], 510, opt + ('G' if log else ''))
        a.SetName(name)
        a.SetLineColor(col)
        a.SetLineWidth(width)
        a.SetLabelSize(lab_size if not line else 0)
        a.SetTitleSize(tit_size)
        a.SetTitleOffset(off)
        a.SetTitle(title)
        a.SetTitleColor(col)
        a.SetLabelColor(col)
        a.SetLabelFont(Draw.Font)
        a.SetTitleFont(Draw.Font)
        do(a.CenterTitle, center)
        a.SetTickSize(tick_size if not line else 0)
        a.SetTickLength(tick_size if not line else 0)
        a.SetNdivisions(0) if line else do_nothing()
        a.SetLabelOffset(l_off)
        a.Draw()
        return Draw.add(a)

    @staticmethod
    def y_axis(x, ymin, ymax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False, center=None):
        return Draw.axis(x, x, ymin, ymax, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log, center)

    @staticmethod
    def x_axis(y, xmin, xmax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False, center=None):
        return Draw.axis(xmin, xmax, y, y, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log, center)

    @staticmethod
    def line(x1, x2, y1, y2, color=1, width=1, style=1):
        line = TCutG(Draw.get_name('l'), 2, array([x1, x2], 'd'), array([y1, y2], 'd'))
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw('same')
        return Draw.add(line)

    @staticmethod
    def tline(x1, x2, y1, y2, color=1, width=1, style=1, ndc=None):
        line = TLine(x1, y1, x2, y2)
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        do(line.SetNDC, ndc)
        line.Draw()
        return Draw.add(line)

    @staticmethod
    def vertical_line(x, ymin=-1e9, ymax=1e9, color=1, w=1, style=1, tline=False):
        return Draw.line(x, x, ymin, ymax, color, w, style) if not tline else Draw.tline(x, x, ymin, ymax, color, w, style)

    @staticmethod
    def horizontal_line(y, xmin=-1e9, xmax=1e9, color=1, w=1, style=1, tline=False, ndc=None):
        return Draw.line(xmin, xmax, y, y, color, w, style) if not tline else Draw.tline(xmin, xmax, y, y, color, w, style, ndc)

    @staticmethod
    def polygon(x, y, line_color=1, width=1, style=1, name=None, fillstyle=None, fill_color=None, opacity=None, show=True):
        if get_object(name) is not None:  # check if name already exists
            get_object(name).Clear()
        line = TCutG(choose(name, Draw.get_name('poly')), len(x) + 1, append(x, x[0]).astype('d'), append(y, y[0]).astype('d'))
        format_histo(line, line_color=line_color, lw=width, line_style=style, fill_color=fill_color, fill_style=fillstyle, opacity=opacity)
        if show:
            line.Draw('l')
            line.Draw('f') if fill_color is not None or fillstyle is not None and fillstyle < 4000 else do_nothing()
        return Draw.add(line)

    @staticmethod
    def box(x1, y1, x2, y2, line_color=1, width=1, style=1, name=None, fillstyle=None, fillcolor=None, opacity=None, show=True):
        return Draw.polygon(*make_box_args(x1, y1, x2, y2), line_color, width, style, name, fillstyle, fillcolor, opacity, show)

    @staticmethod
    def fypolygon(f, x1, x2, y, name=None, n=100, **kwargs):
        x, y = array([(x, f(x)) for x in linspace(x1, x2, n)] + [(x2, y), (x1, y)]).T
        Draw.polygon(x, y, name=name, **kwargs)

    @staticmethod
    def tlatex(x, y, text, name=None, align=20, color=1, size=.05, angle=None, ndc=None, font=None, show=True):
        tlatex = TLatex(x, y, text)
        format_text(tlatex, choose(name, Draw.get_name('t')), align, color, size, angle, ndc, font)
        tlatex.Draw() if show else do_nothing()
        return Draw.add(tlatex)

    @staticmethod
    def date(x, y, align=20, color=1, size=.05, angle=None, font=42, **kwargs):
        return Draw.tlatex(x, y, datetime.now().strftime('%Y-%m-%d %H:%M'), None, align, color, size, angle, ndc=True, font=font, **kwargs)

    @staticmethod
    def arrow(x1, x2, y1, y2, col=1, width=1, opt='<|', size=.005):
        ar = TArrow(x1, y1, x2, y2, size, opt)
        ar.SetLineWidth(width)
        ar.SetLineColor(col)
        ar.SetFillColor(col)
        ar.Draw()
        return Draw.add(ar)

    @staticmethod
    def tpad(tit='', pos=None, fill_col=0, gridx=None, gridy=None, margins=None, transparent=False, logy=None, logx=None, logz=None, lm=None, rm=None, bm=None, tm=None, c=None):
        c.cd() if c is not None else do_nothing()
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(Draw.get_name('pd'), tit, *pos)
        p.SetFillColor(fill_col)
        margins = margins if all(m is None for m in [lm, rm, bm, tm]) else [lm, rm, bm, tm]
        Draw.set_pad_margins(p, *full(4, .1) if margins is None else margins)
        do([p.SetLogx, p.SetLogy, p.SetLogz], [logx, logy, logz])
        do([p.SetGridx, p.SetGridy], [gridx, gridy])
        make_transparent(p) if transparent else do_nothing()
        p.Draw()
        p.cd()
        return Draw.add(p)

    @staticmethod
    def tpavetext(text, x1, x2, y1, y2, font=42, align=0, size=0, angle=0, margin=.05, color=1):
        p = TPaveText(x1, y1, x2, y2, 'ndc')
        p.SetFillColor(0)
        p.SetFillStyle(0)
        p.SetBorderSize(0)
        p.SetMargin(margin)
        t = p.AddText(text)
        format_text(t, 'pave', align, color, size, angle, ndc=True, font=font)
        p.Draw()
        return Draw.add(p)

    @staticmethod
    def stats(fit, y2=None, width=.3, prec='5.1f', names=None):
        names = fit.Names if names is None else names
        c = get_last_canvas()
        tm = .98 - .05 - c.GetTopMargin() if y2 is None else y2
        rm = .98 - c.GetRightMargin()
        p = TPaveStats(rm - width, tm - .06 * (fit.NPars + 1), rm, tm, 'ndc')
        p.SetBorderSize(1)
        p.SetFillColor(0)
        p.SetFillStyle(0)
        leg = p.AddText('Fit Result')
        leg.SetTextFont(42)
        ls = p.GetListOfLines()
        ls.Add(Draw.tlatex(0, 0, '#chi^{{2}} / ndf  = {chi2:{p}} / {ndf}'.format(ndf=fit.Ndf(), chi2=fit.Chi2(), p=prec), size=0, align=0, font=42))
        for i in range(fit.NPars):
            ls.Add(Draw.tlatex(0, 0, '{n}  = {v:{p}} #pm {e:{p}}'.format(n=names[i], v=fit.Parameter(i), e=fit.ParError(i), p=prec), size=0, align=0, font=42))
        p.Draw()
        return Draw.add(p)

    @staticmethod
    def add_stats_entry(h, key, value, form='.2f', line=None):
        s = h.GetListOfFunctions()[0]
        s.SetName(Draw.get_name('st'))
        value = f'{value.n:{form}} #pm {value.s:{form}}' if is_ufloat(value) else f'{value:{form}}'
        text = f'{key} = {value}'
        h.SetStats(0)
        y2, hl = s.GetY2NDC(), (s.GetY2NDC() - s.GetY1NDC()) / s.GetSize()
        s.AddText(text)
        s.SetY1NDC(s.GetY1NDC() - hl)
        [Draw.horizontal_line(y2 - i * hl, s.GetX1NDC(), s.GetX2NDC(), tline=True, ndc=True) for i in make_list(line)] if line is not None else do_nothing()

    @staticmethod
    def frame(pad, xmin, xmax, ymin, ymax, tit, div=None, y_cent=None):
        pad.cd()
        fr = pad.DrawFrame(xmin, ymin, xmax, ymax)
        pad.Modified()
        fr.GetYaxis().SetTitle(tit)
        do(fr.GetYaxis().CenterTitle, y_cent)
        fr.GetYaxis().SetNdivisions(div) if div is not None else do_nothing()
        format_frame(fr)
        Draw.add(fr)

    @staticmethod
    def grid(x_vals, y_vals, width=1, color=1):
        for x in x_vals:
            Draw.line(x, x, min(y_vals), max(y_vals), width=width, color=color)
        for y in y_vals:
            Draw.line(min(x_vals), max(x_vals), y, y, width=width, color=color)

    @staticmethod
    def ellipse(a=1, b=1, x_off=0, y_off=0, color=2, w=2):
        e = TEllipse(x_off, y_off, a, b)
        do(e.SetLineColor, color)
        do(e.SetLineWidth, w)
        e.SetFillStyle(4000)
        e.Draw()
        return Draw.add(e)

    @staticmethod
    def circle(r, x_off=0, y_off=0, color=None, w=None):
        return Draw.ellipse(r, r, x_off, y_off, color, w)

    @staticmethod
    def preliminary(canvas=None, height=.06):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        return Draw.tpavetext('#font[62]{RD42} Preliminary', c.GetLeftMargin(), .5, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=72, align=12, margin=0.04)

    @staticmethod
    def irradiation(irr, canvas=None, height=.06, left=True):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        x1, x2 = (c.GetLeftMargin(), .5) if left else (.5, 1 - c.GetRightMargin())
        return Draw.tpavetext('Irradiation: {}'.format(irr), x1, x2, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=42, align=12, margin=0.04)

    @staticmethod
    def legend(histos, titles, styles=None, x2=None, y2=None, *args, **kwargs):
        leg = Draw.make_legend(x2, y2, *args, **prep_kw(kwargs, nentries=len(histos)))
        for i in range(len(histos)):
            leg.AddEntry(histos[i], titles[i], 'lpf' if styles is None else styles[i] if type(styles) is list else styles)
        leg.Draw('same')
        return leg

    @staticmethod
    def histo(th, show=True, lm=None, rm=None, bm=None, tm=None, m=None, draw_opt=None, w=1, h=1, logx=None, logy=None, logz=None, grid=None, gridy=None, gridx=None, phi=None, theta=None,
              leg=None, canvas=None, sumw2=None, stats=False, **kwargs):
        w += .16 if not Draw.Title and w == 1 else 0  # rectify if there is no title
        th.Sumw2(sumw2) if hasattr(th, 'Sumw2') and sumw2 is not None else do_nothing()
        set_root_output(show)
        c = Draw.canvas(th.GetTitle().split(';')[0], None, None, w, h, logx, logy, logz, gridx or grid, gridy or grid, show=show) if canvas is None else canvas
        Draw.set_pad_margins(c, *[lm, rm, bm, tm] if m is None else m)
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx or grid, gridy or grid])
        do([c.SetPhi, c.SetTheta], [phi, theta])
        th.Draw(draw_opt if draw_opt is not None else 'ap' if is_graph(th) else 'hist' if 'TH' in th.ClassName() else '')
        if leg is not None:
            update_canvas()
            for i_leg in make_list(leg):
                i_leg.Draw('same')
        set_root_output(True)
        if stats or stats is None:
            for i in (th.GetListOfGraphs() if 'Multi' in th.ClassName() else [th]):
                format_statbox(i, **Draw.Stats if stats else Draw.DefaultStats)
        return Draw.add(c, th, leg)[0]

    @staticmethod
    def mode(m=1, **kwargs):
        d = {1: {'tit_size': .05, 'lab_size': .045, 'y_off': 1.35},
             2: {'w': 1.5, 'h': .75, 'tit_size': .06, 'lab_size': .05, 'y_off': .7, 'lm': .08, 'bm': .15},
             3: {'w': 1.5, 'h': .5, 'tit_size': .07, 'lab_size': .06, 'y_off': .5, 'lm': .073, 'bm': .225, 'rm': .03, 'x_tit': 'Time [ns]', 'y_tit': 'Signal [mV]', 'markersize': .5}
             }[m]
        return prep_kw(kwargs, **d)

    def distribution(self, x, binning=None, title='', **kwargs):
        if hasattr(x, 'GetName'):
            th = x
        else:
            th = TH1F(Draw.get_name('h'), title, *choose(binning, find_bins, values=x))
            fill_hist(th, x)
        format_histo(th, **prep_kw(kwargs, **Draw.mode(), fill_color=Draw.FillColor, y_tit='Number of Entries'))
        self.histo(th, **prep_kw(kwargs, stats=None))
        return th

    def function(self, f, title='', c=None, bm=None, lm=None, rm=None, show=True, logy=None, w=1, h=1, stats=None, draw_opt=None, grid=None, **kwargs):
        format_histo(f, title=title, **prep_kw(kwargs, **Draw.mode()))
        c = get_last_canvas() if draw_opt and 'same' in draw_opt and c is None else c
        self.histo(f, show=show, bm=bm, lm=lm, rm=rm, logy=logy, w=w, h=h, stats=stats, draw_opt=draw_opt, canvas=c, grid=grid)
        return f

    def graph(self, x, y=None, title='', c=None, bm=None, show=True, bin_labels=None, **kwargs):
        g = x if y is None else Draw.make_tgrapherrors(x, y)
        format_histo(g, title=title, **prep_kw(kwargs, **Draw.mode(), fill_color=Draw.FillColor))
        set_bin_labels(g, bin_labels)
        self.histo(g, show=show, bm=choose(bm, .24 if bin_labels else None), canvas=c, **kwargs)
        return g

    def profile(self, x, y=None, binning=None, title='', thresh=.02, show=True, graph=False, **kwargs):
        if y is None:
            p = x
        else:
            x, y = array(x, dtype='d'), array(y, dtype='d')
            p = TProfile(Draw.get_name('p'), title, *choose(binning, find_bins, values=x, thresh=thresh))
            fill_hist(p, x, y)
        p = self.make_graph_from_profile(p) if graph else p
        format_histo(p, **prep_kw(kwargs, **Draw.mode(), fill_color=Draw.FillColor, stats=None))
        self.histo(p, show=show, **kwargs)
        return p

    def prof2d(self, x, y=None, zz=None, binning=None, title='', lm=None, rm=.17, bm=None, show=True, phi=None, theta=None, draw_opt='colz', stats=None,
               rot=None, mirror=None, centre=None, **kwargs):
        if y is None:
            p = x
        else:
            dflt_bins = make_bins(min(x), max(x), sqrt(len(x))) + make_bins(min(y), max(y), sqrt(len(y))) if binning is None else None
            p = TProfile2D(Draw.get_name('p2'), title, *choose(binning, dflt_bins))
            fill_hist(p, x, y, uarr2n(zz))
        p = self.rotate_2d(p, rot)
        p = self.flip_2d(p, mirror)
        rx, ry = get_2d_centre_ranges(p, centre)
        format_histo(p, **prep_kw(kwargs, **Draw.mode(), z_off=1.2, pal=55, stats=stats, x_range=rx, y_range=ry))
        set_statbox(entries=True, w=.25) if stats is None else do_nothing()
        self.histo(p, show=show, lm=lm, rm=rm, bm=bm, phi=phi, theta=theta, draw_opt=draw_opt, stats=choose(stats, True), **kwargs)
        return p

    def histo_2d(self, x, y=None, binning=None, title='', lm=None, rm=.17, bm=None, tm=None, show=True, stats=None, canvas=None,
                 rot=None, mirror=None, centre=None, **kwargs):
        if y is None:
            th = x
        else:
            x, y = array(x, dtype='d'), array(y, dtype='d')
            dflt_bins = make_bins(min(x), max(x), sqrt(x.size)) + make_bins(min(y), max(y), sqrt(x.size)) if binning is None else None
            th = TH2F(Draw.get_name('h2'), title, *choose(binning, dflt_bins))
            fill_hist(th, x, y)
        th = self.rotate_2d(th, rot)
        th = self.flip_2d(th, mirror)
        rx, ry = get_2d_centre_ranges(th, centre)
        format_histo(th, **prep_kw(kwargs, **Draw.mode(), z_off=1.2, z_tit='Number of Entries', pal=55, stats=stats, x_range=rx, y_range=ry))
        set_statbox(entries=True, w=.25) if stats is None else do_nothing()
        self.histo(th, show=show, lm=lm, rm=rm, bm=bm, tm=tm, stats=choose(stats, True), canvas=canvas, **prep_kw(kwargs, draw_opt='colz'))
        return th

    def efficiency(self, x, e, binning=None, title='Efficiency', lm=None, show=True, **kwargs):
        binning = choose(binning, make_bins, min(x), max(x), (max(x) - min(x)) / sqrt(x.size))
        p = self.profile(x, e, binning, show=False)
        x = get_hist_args(p)
        y = array([calc_eff(p0 * n, n) for p0, n in [[p.GetBinContent(ibin), p.GetBinEntries(ibin)] for ibin in range(1, p.GetNbinsX() + 1)]])
        return self.graph(x, y, title, lm=lm, show=show, y_tit='Efficiency [%]', **kwargs)

    def pull(self, h, binning=None, **kwargs):
        x = get_hist_vec(h)
        self.distribution(x, choose(binning, find_bins, values=x[x != 0], lfac=.2, rfac=.2, n=2), x_tit=h.GetYaxis().GetTitle(), **kwargs)
        return mean_sigma(x[x != 0])

    def stack(self, histos, title, leg_titles, scale=False, draw_opt='nostack', show=True, fill=None, w=.2, *args, **kwargs):
        s = THStack(Draw.get_name('s'), title)
        leg = Draw.make_legend(nentries=len(histos), w=w)
        for h, tit in zip(histos, leg_titles):
            s.Add(h, 'hist')
            leg.AddEntry(h, tit, choose('lf', 'l', fill))
            color = self.get_color(len(histos))
            format_histo(h, color=color, fill_color=choose(color, 0, fill), fill_style=choose(1001, 4000, fill), stats=0, *args, **kwargs)
            if scale:
                h.Scale(1 / h.GetMaximum())
        h0 = histos[0]
        format_histo(s, draw_first=True, x_tit=h0.GetXaxis().GetTitle(), y_tit=h0.GetYaxis().GetTitle(), y_off=h0.GetYaxis().GetTitleOffset())
        self.histo(s, draw_opt=draw_opt, leg=leg, lm=get_last_canvas().GetLeftMargin(), show=show)
        return s

    def multigraph(self, graphs, title='', leg_titles=None, bin_labels=None, draw_opt='p', wleg=.2, **kwargs):
        if hasattr(graphs, 'GetName'):
            m, g = graphs, graphs.GetListOfGraphs()[0]
        else:
            g = graphs[0]
            m = TMultiGraph(Draw.get_name('mg'), ';'.join([title, g.GetXaxis().GetTitle(), g.GetYaxis().GetTitle()]))
            for i, g in enumerate(graphs):
                m.Add(g, draw_opt)
                format_histo(g, **prep_kw(kwargs, color=self.get_color(len(graphs)), stats=False))
        y_range = ax_range(get_graph_y(graphs, err=False), 0, .3, .6)
        format_histo(m, draw_first=True, **prep_kw(kwargs, **Draw.mode(1, y_off=g.GetYaxis().GetTitleOffset()), y_range=y_range, x_tit=choose('', None, bin_labels)))
        set_bin_labels(m, bin_labels)
        leg = self.legend(graphs, leg_titles, 'p', w=wleg) if leg_titles else None
        self.histo(m, leg=leg, **prep_kw(kwargs, bm=choose(.26, None, bin_labels), draw_opt='ap'))
        return m

    def pie(self, labels, values=None, colors=None, title='', offset=0, show=True, flat=False, draw_opt=None, **kwargs):
        labels, (values, colors) = (labels.keys(), array(list(labels.values())).T) if values is None else (labels, values, colors)
        pie = TPie(self.get_name('pie'), title, len(labels), array(values, 'f'), array(colors, 'i'))
        for i, label in enumerate(labels):
            pie.SetEntryRadiusOffset(i, offset)
            pie.SetEntryLabel(i, label)
        format_pie(pie, **kwargs)
        draw_opt = choose(draw_opt, f'{"" if flat else "3d"}rsc')
        self.histo(pie, draw_opt=draw_opt, show=show)
        return pie

    def prof2hist(self, p):
        bins = get_2d_bins(p)
        h, nx, ny = self.histo_2d([], [], bins, show=False), bins[0], bins[2]
        xax, yax = p.GetXaxis(), p.GetYaxis()
        [h.Fill(xax.GetBinCenter(ix), yax.GetBinCenter(iy)) for ix in range(1, nx + 2) for iy in range(1, ny + 2) for _ in range(_get_2d_bin_entries(p, ix, iy, nx))]
        return h

    @staticmethod
    def info(txt, c=None, size=.04):
        c = (get_last_canvas() if c is None else c).cd()
        Draw.tlatex(c.GetLeftMargin() + .02, 1 - (c.GetTopMargin() + .02), txt, align=13, ndc=True, size=size)

    @staticmethod
    def bin_numbers(h, show=True):
        if show:
            x, y = get_2d_bins(h, arr=True)
            dx, dy = diff(x)[0] / 2, diff(y)[0] / 2
            [Draw.tlatex(x[m] + dx, y[n] + dy, str((x.size - 1) * n + m)) for n in range(y.size - 1) for m in range(x.size - 1)]
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region OPERATIONS
    def operate(self, h, f, *args, **kwargs):
        h0, h = h, self.prof2d([], [], [], get_2d_bins(h), show=False)
        prof = 'Profile' in h0.ClassName()
        x, n = get_2d_hist_vec(h0, err=False, flat=False), get_2d_bin_entries(h0) if prof else 1
        set_2d_values(h, f(x * n, *args, **kwargs))
        set_2d_entries(h, f(n, *args, **kwargs)) if prof else do_nothing()
        h.SetEntries(int(h0.GetEntries()))
        format_histo(h, z_range=[h0.GetMinimum(), h0.GetMaximum()], **{f'{i}_tit': getattr(h0, f'Get{i.title()}axis')().GetTitle() for i in ['x', 'y', 'z']}, ncont=h0.GetContour())
        return h

    def rotate_2d(self, h, n=2):
        return self.operate(h, rot90, n) if n is not None else h

    def flip_2d(self, h, axis=0):
        return self.operate(h, flip, axis=axis) if axis is not None else h
    # endregion OPERATIONS
    # ----------------------------------------

    # ----------------------------------------
    # region CREATE
    @staticmethod
    def make_histo(title, bins):
        h = TH1F(Draw.get_name('h'), title, *bins)
        return Draw.add(h)

    @staticmethod
    def make_f(name, function, xmin=0, xmax=1, pars=None, limits=None, fix=None, npx=None, **kwargs):
        f = TF1(choose(name, Draw.get_name('f')), function, xmin, xmax)
        f.SetParameters(*pars) if pars is not None else do_nothing()
        [f.SetParLimits(i, *lim) for i, lim in enumerate(limits)] if limits else do_nothing()
        [f.FixParameter(i, value) for i, value in enumerate(make_list(fix))] if fix is not None else do_nothing()
        do(f.SetNpx, npx)
        format_histo(f, **kwargs)
        return Draw.add(f)

    @staticmethod
    def make_tf1(name, f, xmin=0, xmax=1, color=None, w=None, style=None, title='', npx=None, *args, **kwargs):
        def tmp(x, pars):
            _ = pars
            return f(x[0], pars, *args, **kwargs) if 'pars' in signature(f).parameters else f(x[0], *args, **kwargs)

        Draw.add(tmp)
        f0 = TF1(choose(name, Draw.get_name('f')), tmp, xmin, xmax)
        do(f0.SetNpx, npx)
        format_histo(f0, title, line_color=color, line_style=style, lw=w)
        return Draw.add(f0)

    @staticmethod
    def make_tgrapherrors(x=None, y=None, **kwargs):
        if len(list(x)) != len(list(y)) or not len(x):
            return warning('Arrays have different size!')
        x, y = array(x), array(y)
        asym = len(x.shape) == 2 or len(y.shape) == 2
        s, utypes, has_ers = len(x), [type(v[0]) in [Variable, AffineScalarFunc] for v in [x, y]], [len(v.shape) > 1 for v in [x, y]]
        ex, ey = [array([[v.s for v in vals]] if is_u and not asym else vals[:, 1:3].T if has_e else zeros((2, s)) if asym else [zeros(s)], 'd') for vals, is_u, has_e in zip([x, y], utypes, has_ers)]
        x, y = [array([v.n for v in vals] if utype else vals[:, 0] if has_e else vals, 'd') for vals, utype, has_e in zip([x, y], utypes, has_ers)]
        g = (TGraphAsymmErrors if asym else TGraphErrors)(s, x, y, *array(ex.tolist()), *array(ey.tolist()))  # doesn't work without double conversion...
        format_histo(g, Draw.get_name('g'), **prep_kw(kwargs, marker=20, markersize=1))
        return Draw.add(g)

    @staticmethod
    def make_graph_from_profile(p):
        x, y = get_hist_vecs(p)
        return Draw.make_tgrapherrors(x[y != 0], y[y != 0], x_tit=p.GetXaxis().GetTitle(), y_tit=p.GetYaxis().GetTitle())

    @staticmethod
    def make_legend(x2=None, y2=None, w=.25, nentries=2, scale=1, d=.01, y1=None, x1=None, clean=False, margin=.25, cols=None, fix=False, bottom=False, left=False):
        use_margins = y2 is None
        h = nentries * .05 * scale
        x2, y2 = get_stat_margins(None, x2, y2, d, bottom, left, h, w)
        x1 = choose(x1, x2 - w)
        y1 = choose(y1, y2 - h)
        if not use_margins:
            y1 += .07 if not Draw.Title and y1 + h > .8 and not fix else 0
            y1 -= .07 if not Draw.Legend and y1 < .3 and not fix else 0
        leg = TLegend(x1, max(y1, 0), x1 + w, min(y1 + h, 1))
        leg.SetName(Draw.get_name('l'))
        leg.SetTextFont(Draw.Font)
        leg.SetMargin(margin)
        do(leg.SetNColumns, cols)
        if clean:
            leg.SetLineWidth(2)
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetTextAlign(12)
        Draw.add(leg)
        return leg
    # endregion CREATE
    # ----------------------------------------

# END OF CLASS ---------------------------


# ----------------------------------------
# region FORMATTING
def format_histo(histo, name=None, title=None, x_tit=None, y_tit=None, z_tit=None, marker=None, color=None, line_color=None, line_style=None, markersize=None, x_off=None, y_off=None, z_off=None,
                 lw=None, fill_color=None, fill_style=None, stats=None, tit_size=None, lab_size=None, l_off_y=None, l_off_x=None, draw_first=False, x_range=None, y_range=None, z_range=None,
                 sumw2=None, do_marker=True, style=None, ndivx=None, ndivy=None, ncont=None, tick_size=None, t_ax_off=None, center_y=False, center_x=False, yax_col=None, normalise=None, pal=None,
                 rebin=None, y_ticks=None, x_ticks=None, z_ticks=None, opacity=None, center_tit=None, **kwargs):
    _ = kwargs
    h = histo
    if draw_first:
        set_root_output(False)
        h.Draw('nostack' if h.ClassName() == 'THStack' else 'a')
        set_root_output(True)
    do(h.SetTitle, title)
    do(h.SetName, name)
    do(set_palette, pal)
    if normalise is not None:
        y_tit = 'Frequency' if 'Number' in choose(y_tit, h.GetYaxis().GetTitle()) else choose(y_tit, h.GetYaxis().GetTitle())
        normalise_histo(h)
    try:
        do(h.SetStats, stats)
    except AttributeError or ReferenceError:
        pass
    do(h.Rebin, rebin) if hasattr(h, 'Rebin') else do_nothing()
    # markers
    try:
        if do_marker:
            do(h.SetMarkerStyle, marker)
            do(h.SetMarkerColor, color)
            do(h.SetMarkerSize, markersize)
    except AttributeError or ReferenceError:
        pass
    # lines/fill
    try:
        h.SetLineColor(line_color) if line_color is not None else h.SetLineColor(color) if color is not None else do_nothing()
        do(h.SetLineWidth, lw)
        do(h.SetLineStyle, line_style)
        h.SetFillColor(fill_color) if fill_color is not None and opacity is None else do_nothing()
        h.SetFillColorAlpha(fill_color, opacity) if fill_color is not None and opacity is not None else do_nothing()
        h.SetFillStyle(fill_style) if fill_style is not None else do_nothing()
        h.SetFillStyle(style) if style is not None else do_nothing()
        h.SetContour(ncont) if ncont is not None else do_nothing()
    except AttributeError or ReferenceError:
        pass
    # axes
    try:
        x_args = [x_tit, x_off, tit_size, choose(center_tit, center_x), lab_size, l_off_x, x_range, ndivx, choose(x_ticks, tick_size), ]
        y_args = [y_tit, y_off, tit_size, choose(center_tit, center_y), lab_size, l_off_y, y_range, ndivy, choose(y_ticks, tick_size), yax_col]
        z_args = [z_tit, z_off, tit_size, False, lab_size, None, z_range, None, choose(z_ticks, tick_size)]
        for i, name in enumerate(['X', 'Y', 'Z']):
            format_axis(getattr(h, 'Get{}axis'.format(name))(), h, *[x_args, y_args, z_args][i])
    except AttributeError or ReferenceError:
        pass
    set_time_axis(h, off=t_ax_off) if t_ax_off is not None else do_nothing()
    do(h.Sumw2, sumw2) if hasattr(h, 'Sumw2') else do_nothing()
    update_canvas()
    return h


def set_statbox(x2=None, y2=None, h=None, w=.3, entries=False, m=False, rms=False, all_stat=False, fit=False, center_x=False, center_y=False, form=None, stats=True):
    Draw.Stats = {'x2': x2, 'y2': y2, 'h': h, 'w': w, 'entries': entries, 'm': m, 'rms': rms, 'all_stat': all_stat, 'fit': fit, 'center_x': center_x, 'center_y': center_y, 'form': form}
    return stats


def set_entries():
    Draw.Stats = {'entries': True, 'w': .2}
    return True


def get_stat_margins(c=None, x2=None, y2=None, d=.01, bottom=False, left=False, h=0, w=0):
    c = choose(c, get_last_canvas(warn=False))
    r = c.GetWindowHeight() / c.GetWindowWidth()
    x2 = choose(x2, c.GetLeftMargin() + w + d * r if left else 1 - c.GetRightMargin() - d * r)
    y2 = choose(y2, c.GetBottomMargin() + h + d if bottom else 1 - c.GetTopMargin() - d)
    return x2, y2


def format_statbox(th, x2=None, y2=None, d=.01, h=None, w=.3, entries=False, m=False, rms=False, all_stat=False, fit=False, center_x=False,
                   center_y=False, bottom=False, left=False, form=None, c=None):
    c = choose(c, get_last_canvas(warn=False))
    update_canvas(c)
    f = None if 'TF1' in th.ClassName() else next((o for o in th.GetListOfFunctions() if 'TF1' in o.ClassName()), None)
    if 'TGraph' in th.ClassName() and fit and f:
        gStyle.SetOptFit(True)
    p = None if 'TF1' in th.ClassName() else next((o for o in th.GetListOfFunctions() if 'Pave' in o.ClassName()), None)
    if p is not None:
        r = c.GetWindowHeight() / c.GetWindowWidth()
        stats = ones(3, 'i') if all_stat else array([rms, m, entries], 'i')
        fit_pars = f.GetNpar() + 1 if fit and f is not None else 0
        h = choose(h, .05 / r * (stats.nonzero()[0].size + fit_pars))
        x2, y2 = get_stat_margins(c, x2, y2, d, bottom, left, h, w)
        cx, cy = mean([1 - c.GetRightMargin(), c.GetLeftMargin()]), mean([1 - c.GetTopMargin(), c.GetBottomMargin()])
        p.SetX1NDC(cx - w / 2 if center_x else x2 - w)
        p.SetX2NDC(cx + w / 2 if center_x else x2)
        p.SetY1NDC(cy - h / 2 if center_y else y2 - h)
        p.SetY2NDC(cy + h / 2 if center_y else y2)
        p.SetOptStat(int('100000{}{}{}0'.format(*stats * 2)))
        p.SetOptFit(int(fit))
        p.SetStatFormat(choose(form, '1.1f'))
        p.SetFitFormat(form) if form is not None else do_nothing()
        p.Draw()
        update_canvas(c)


def format_axis(axis, h, title, tit_offset, tit_size, centre_title, lab_size, label_offset, limits, ndiv, tick_size, color=None):
    do(axis.SetTitle, title)
    do(axis.SetTitleOffset, tit_offset)
    do(axis.SetTitleSize, tit_size)
    axis.CenterTitle(centre_title)
    do(axis.SetLabelSize, lab_size)
    do(axis.SetLabelOffset, label_offset)
    if limits is not None:
        axis.SetLimits(*limits) if is_graph(h) and 'xaxis' in axis.GetName() else axis.SetRangeUser(*limits)
    do(axis.SetNdivisions, ndiv)
    do(axis.SetTickSize, tick_size)
    do(axis.SetTitleColor, color)
    do(axis.SetLabelColor, color)
    do(axis.SetAxisColor, color)


def format_pie(pie, h=None, r=None, text_size=None, angle3d=None, angle_off=None, label_format=None):
    do([pie.SetHeight, pie.SetRadius], [h, r])
    do(pie.SetTextSize, text_size)
    do(pie.SetAngle3D, angle3d)
    do(pie.SetLabelFormat, label_format)
    do(pie.SetAngularOffset, angle_off)


def format_text(t, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None):
    t.SetName(name)
    t.SetTextAlign(align)
    t.SetTextColor(color)
    t.SetTextSize(size)
    do(t.SetTextAngle, angle)
    do(t.SetTextFont, font)
    do(t.SetNDC, ndc)
    return t


def format_frame(frame):
    fr = frame
    fr.GetYaxis().SetTitleSize(.06)
    fr.GetYaxis().SetTitleOffset(.6)
    fr.GetYaxis().SetLabelSize(.06)
    fr.SetTitleSize(.05)
    fr.GetXaxis().SetTickLength(0)
    fr.GetXaxis().SetLabelOffset(99)
    fr.SetLineColor(0)
    fr.GetXaxis().SetTimeDisplay(1)
# endregion FORMATTING
# ----------------------------------------


def fill_hist(h, x, y=None, zz=None):
    x, y, zz = array(x).astype('d'), array(y).astype('d'), array(zz).astype('d')
    if len(x.shape) > 1:
        y = array(x[:, 1])
        x = array(x[:, 0])
    if not x.size:  # return if there are no entries
        return
    if h.ClassName() == 'TProfile2D':
        for i in range(x.size):
            h.Fill(x[i], y[i], zz[i])
    elif 'TH1' in h.ClassName():
        h.FillN(x.size, x, ones(x.size))
    elif any(name in h.ClassName() for name in ['TH2', 'TProfile']):
        h.FillN(x.size, x, y, ones(x.size))
    else:
        h.FillN(x.size, x, y, zz, ones(x.size))


def set_2d_ranges(h, dx, dy):
    # find centers in x and y
    xmid, ymid = [(p.GetBinCenter(p.FindFirstBinAbove(0)) + p.GetBinCenter(p.FindLastBinAbove(0))) / 2 for p in [h.ProjectionX(), h.ProjectionY()]]
    format_histo(h, x_range=[xmid - dx, xmid + dx], y_range=[ymid - dy, ymid + dx])


def find_bins(values, lfac=.2, rfac=.2, q=.02, n=1):
    width, (xmin, xmax) = freedman_diaconis(values) * n, find_range(values, lfac, rfac, q)
    bins = arange(xmin, xmax + width, width)
    return [bins.size - 1, bins]


def find_range(values, lfac=.2, rfac=.2, q=.02):
    return ax_range(*quantile(values, [q, 1 - q]), lfac, rfac)


def fix_chi2(g, prec=.01, show=True):
    it = 0
    error = 2
    chi2 = 0
    fit = None
    while abs(chi2 - 1) > prec and it < 20:
        for i in range(g.GetN()):
            g.SetPointError(i, g.GetErrorX(i), error)
        fit = g.Fit('pol0', 'qs{}'.format('' if show else 0))
        chi2 = fit.Chi2() / fit.Ndf()
        error += .5 ** it * sign(chi2 - 1)
        it += 1
    return None if fit is None else FitRes(fit)


def make_darray(values):
    return array([v.n for v in values] if is_ufloat(values[0]) else values, dtype='d')


def set_bin_labels(g, labels):
    if labels is not None:
        for i, label in enumerate(labels):
            g.GetXaxis().SetBinLabel(g.GetXaxis().FindBin(i), label)
        update_canvas()


def make_box_args(x1, y1, x2, y2):
    return array([[x1, x1, x2, x2], [y1, y2, y2, y1]])


def make_star(cx=0, cy=0, r=1, n=5):
    coods = pol2cart(tile([r, r / (2 * cos(pi / n) + 1)], n), linspace(0, 2 * pi, 2 * n, endpoint=False) - pi / 2)
    return (coods.T + array([cx, cy])).T


def set_titles(status=True):
    gStyle.SetOptTitle(status)


def get_graph_vecs(g, err=True):
    return get_graph_x(g, err), get_graph_y(g, err)


def get_graph_x(g, err=True):
    return make_ufloat(frombuffer(g.GetX()), frombuffer(g.GetEX())) if err else frombuffer(g.GetX())


def get_graph_y(g, err=True):
    if is_iter(g):
        return array([v for ig in g for v in get_graph_y(ig, err)])
    return make_ufloat(frombuffer(g.GetY()), frombuffer(g.GetEY())) if err else frombuffer(g.GetY())


def get_hist_vec(p, err=True):
    return array([ufloat(p.GetBinContent(ibin), p.GetBinError(ibin)) if err else p.GetBinContent(ibin) for ibin in range(1, p.GetNbinsX() + 1)])


def get_hist_args(h, err=True, raw=False, axis='X'):
    ax = getattr(h, f'Get{axis.title()}axis')()
    if raw:
        return array([ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins() + 2)], 'd')
    return array([ufloat(ax.GetBinCenter(ibin), ax.GetBinWidth(ibin) / 2) if err else ax.GetBinCenter(ibin) for ibin in range(1, ax.GetNbins() + 1)])


def get_hist_vecs(p, err=True, raw=False):
    if type(p) in [ndarray, list]:
        return array([tup for ip in p for tup in array(get_hist_vecs(ip, err, raw)).T]).T
    return get_hist_args(p, err, raw), get_hist_vec(p, err)


def get_h_values(h):
    return get_graph_y(h) if 'Graph' in h.ClassName() else get_hist_vec(h)


def get_h_args(h):
    return get_graph_x(h) if 'Graph' in h.ClassName() else get_hist_args(h)


def get_2d_hist_vec(h, err=True, flat=True, zero_supp=True):
    xbins, ybins = range(1, h.GetNbinsX() + 1), range(1, h.GetNbinsY() + 1)
    values = array([ufloat(h.GetBinContent(xbin, ybin), h.GetBinError(xbin, ybin)) for ybin in ybins for xbin in xbins])
    values = values if err else array([v.n for v in values])
    return (values[values != 0] if zero_supp else values) if flat else values.reshape(len(ybins), len(xbins))


def get_2d_bins(h, arr=False):
    x, y = [get_hist_args(h, raw=True, axis=ax) for ax in ['X', 'Y']]
    return [x, y] if arr else [x.size - 1, x, y.size - 1, y]


def set_2d_values(h, arr):
    [h.SetBinContent(ix + 1, iy + 1, arr[iy, ix]) for ix in range(arr.shape[1]) for iy in range(arr.shape[0])]


def set_2d_entries(h, arr):
    ny, nx = arr.shape
    [h.SetBinEntries((nx + 2) * (iy + 1) + (ix + 1), arr[iy, ix]) for ix in range(nx) for iy in range(ny)]


def _get_2d_bin_entries(h, ix, iy, nx):
    return int(h.GetBinEntries((nx + 2) * iy + ix))


def get_2d_bin_entries(h, flat=False):
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    entries = array([[_get_2d_bin_entries(h, ix, iy, nx) for ix in range(1, nx + 1)] for iy in range(1, ny + 1)])
    return entries.flatten() if flat else entries


def get_2d_args(h):
    return array([[getattr(h, 'Get{}axis'.format(ax))().GetBinCenter(ibin) for ibin in range(1, getattr(h, 'GetNbins{}'.format(ax))() + 1)] for ax in ['X', 'Y']])


def get_2d_vecs(h, err=True, flat=False):
    return get_2d_args(h), get_2d_hist_vec(h, err, flat)


def get_h_entries(h):
    return array([h.GetBinEntries(ibin) for ibin in range(1, h.GetNbinsX() + 1)])


def scale_multigraph(mg, val=1, to_low_flux=False):
    graphs = mg.GetListOfGraphs()
    first_graph = graphs[0]
    scale = scale_graph(first_graph, val=val, to_low_flux=to_low_flux)
    for gr in graphs[1:]:
        scale_graph(gr, scale)
    for i, l in enumerate(first_graph.GetListOfFunctions()):
        y, ey = first_graph.GetY()[i], first_graph.GetErrorY(i)
        l.SetY(y - ey - .003)


def scale_graph(gr, scale=None, val=1, to_low_flux=False):
    x, y = get_graph_vecs(gr)
    if scale is None:
        m, s = mean_sigma(y, err=False)
        scale = val / (y[where(x == min(x))[0]] if to_low_flux else m)
    for i in range(x.size):
        gr.SetPoint(i, gr.GetX()[i], gr.GetY()[i] * scale)
        gr.SetPointError(i, gr.GetErrorX(i), gr.GetErrorY(i) * scale) if 'Error' in gr.ClassName() else do_nothing()
    return scale


def get_quantile(h, q):
    quantiles = make_list(q)
    v = zeros(quantiles.size)
    h.GetQuantiles(v.size, v, quantiles)
    return v[0] if v.size == 1 else v


def markers(i):
    return ((list(range(20, 24)) + [29, 33, 34]) * 2)[i]


def set_palette(pal):
    gStyle.SetPalette(pal)


def is_graph(h):
    return 'Graph' in h.ClassName()


def update_canvas(c=None):
    c = choose(c, get_last_canvas(warn=False))
    if c is not None:
        c.Modified()
        c.Update()
    return c


def show_colors(colors):
    n = len(colors)
    c = Draw.canvas(divide=(int(ceil(sqrt(n))), int(ceil(sqrt(n)))))
    for i, col in enumerate(colors, 1):
        c.cd(i)
        Draw.box(0, 0, 1, 1, fillstyle=1001, fillcolor=col)
        Draw.tlatex(.5, .5, str(i - 1), align=22, size=.2)


def show_wheel():
    from ROOT import TColorWheel
    t = TColorWheel()
    t.Draw()
    Draw.add(t)


def show_line_styles():
    Draw.canvas(w=1.5, title='Line Styles')
    for i in range(1, 11):
        Draw.horizontal_line(1 / 11 * i, 0.1, .95, w=2, style=i)
        Draw.tlatex(.07, 1 / 11 * i, str(i), align=32)


def ax_range(low=None, high=None, fl=0, fh=0, h=None, rnd=False, thresh=None):
    if type(low) in [list, ndarray]:
        utypes = [Variable, AffineScalarFunc]
        if len(low) == 2:
            return ax_range(low[0], low[1], fl, fh)
        m, s = mean_sigma(low, err=0)
        v = low[absolute(low - m) < thresh * s] if thresh is not None else low
        return ax_range(min(v).n if type(v[0]) in utypes else min(v), max(v).n if type(v[0]) in utypes else max(v), fl, fh, rnd=rnd)
    if h is not None:
        lo, hi = choose(thresh, low), choose(thresh, high)
        if 'TH2' in h.ClassName() or '2D' in h.ClassName():
            axes = enumerate([h.GetXaxis(), h.GetYaxis()], 1)
            return [ax_range(ax.GetBinCenter(h.FindFirstBinAbove(lo, i)), ax.GetBinCenter(h.FindLastBinAbove(hi, i)), fl, fh, rnd=rnd) for i, ax in axes]
        return ax_range(h.GetBinCenter(h.FindFirstBinAbove(lo)), h.GetBinCenter(h.FindLastBinAbove(hi)), fl, fh, rnd=rnd)
    d = abs(high - low)
    l, h = low - d * fl, high + d * fh
    return [int(l), int(ceil(h))] if rnd else [l, h]


def set_drawing_range(h, legend=True, lfac=None, rfac=None, thresh=None):
    for i in range(1, 4):
        h.SetBinContent(i, 0)
    thresh = choose(thresh, .05 * h.GetMaximum())
    range_ = [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(thresh), h.FindLastBinAbove(thresh)]]
    lfac = lfac if lfac is not None else .2
    rfac = rfac if rfac is not None else .55 if legend else .1
    h.GetXaxis().SetRangeUser(*ax_range(range_, lfac, rfac))


def normalise_histo(histo, x_range=None, from_min=False):
    h = histo
    x_axis = h.GetXaxis()
    x_axis.SetRangeUser(*x_range) if x_range is not None else do_nothing()
    min_bin = h.GetMinimumBin() if from_min else 1
    integral = h.Integral(min_bin, h.GetNbinsX() - 2)
    return scale_histo(h, integral)


def normalise_bins(h):
    px = h.ProjectionX()
    for xbin in range(h.GetNbinsX()):
        for ybin in range(h.GetNbinsY()):
            h.SetBinContent(xbin, ybin, h.GetBinContent(xbin, ybin) / (px.GetBinContent(xbin) if px.GetBinContent(xbin) else 1))
    update_canvas()


def make_bins(min_val, max_val=None, bin_width=1, last=False, n=None, off=0):
    bins = array(min_val, 'd')
    if type(min_val) not in [ndarray, list]:
        min_val, max_val = choose(min_val, 0, decider=max_val), choose(max_val, min_val)
        bins = append(arange(min_val, max_val, bin_width, dtype='d'), max_val if last else []) if n is None else linspace(min_val, max_val, int(n) + 1, endpoint=True)
    return [bins.size - 1, bins + off]


def set_z_range(zmin, zmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetZaxis().SetRangeUser(zmin, zmax)


def set_axes_range(xmin, xmax, ymin, ymax, c=None):
    set_x_range(xmin, xmax, c)
    set_y_range(ymin, ymax, c)
    update_canvas()


def set_x_range(xmin, xmax, c=None):
    c = choose(c, get_last_canvas())
    h = c.GetListOfPrimitives()[1]
    h.GetXaxis().SetRangeUser(xmin, xmax)


def set_y_range(ymin, ymax, c=None):
    c = choose(c, get_last_canvas())
    h = c.GetListOfPrimitives()[1]
    h.GetYaxis().SetRangeUser(ymin, ymax)


def get_last_canvas(warn=True):
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        warning('There is no canvas is in the list...', prnt=warn)


def close_last_canvas():
    get_last_canvas().Close()


def get_object(name):
    if name is not None:
        o = gROOT.FindObject(name)
        return None if o.__class__.Class_Name() == 'TObject' else o


def set_time_axis(histo, form='%H:%M', off=0):
    histo.GetXaxis().SetTimeFormat(form)
    histo.GetXaxis().SetTimeOffset(-off - 3600 if off else 0)
    histo.GetXaxis().SetTimeDisplay(1)


def find_mpv_fwhm(histo, bins=15):
    max_bin = histo.GetMaximumBin()
    fit = TF1('fit', 'gaus', 0, 500)
    histo.Fit('fit', 'qs0', '', histo.GetBinCenter(max_bin - bins), histo.GetBinCenter(max_bin + bins))
    mpv = ufloat(fit.GetParameter(1), fit.GetParError(1))
    fwhm = histo.FindLastBinAbove(fit(mpv.n) / 2) - histo.FindFirstBinAbove(fit(mpv.n) / 2)
    return mpv, fwhm, mpv / fwhm


def get_fw_center(h):
    return mean(get_fwhm(h, ret_edges=True))  # center of FWHM as MPV


def get_fwhm(h, fit_range=.8, ret_edges=False, err=True):
    x, y = [f(get_hist_vec(h, err=False)) for f in [argsort, sorted]]
    x, ymax = (x[-1] + 1, y[-1]) if y[-1] < 2 * y[-2] else (x[-2] + 1, y[-2])
    fit_range = [f(ymax * fit_range) for f in [h.FindFirstBinAbove, h.FindLastBinAbove]]
    fit_range = fit_range if diff(fit_range)[0] > 5 else (x + array([-5, 5])).tolist()
    half_max = FitRes(h.Fit('gaus', 'qs0', '', *[h.GetBinCenter(i) for i in fit_range]))[0] * .5  # fit the top with a gaussian to get better maxvalue
    half_max = ufloat(1, .05) * ymax / 2 if half_max > .9 * ymax else ufloat(1, .02) * half_max  # half max must be lower than max ...
    low, high = [ufloat(v.n, v.s + abs(v.n - i.n)) for v, i in zip(_get_fwhm(h, half_max), _get_fwhm(h, half_max - half_max.s))]
    return ((low, high) if err else (low.n, high.n)) if ret_edges else high - low


def _get_fwhm(h, half_max):
    blow, bhigh, w = h.FindFirstBinAbove(half_max.n), h.FindLastBinAbove(half_max.n), h.GetBinWidth(1)
    low = get_x(h.GetBinCenter(blow - 1), h.GetBinCenter(blow), h.GetBinContent(blow - 1), h.GetBinContent(blow), half_max)
    high = get_x(h.GetBinCenter(bhigh), h.GetBinCenter(bhigh + 1), h.GetBinContent(bhigh), h.GetBinContent(bhigh + 1), half_max)
    return low, high


def fit_fwhm(h, fitfunc='gaus', show=False, fit_range=.8):
    low, high = get_fwhm(h, fit_range, ret_edges=True)
    return FitRes(h.Fit(fitfunc, 'qs{}'.format('' if show else 0), '', low.n, high.n))


def get_f_fwhm(f: TF1):
    half_max = f.GetMaximum() / 2
    return f.GetX(half_max, f.GetMaximumX(), 1e9) - f.GetX(half_max)


def scale_histo(histo, value=None, to_max=False, x_range=None):
    h = histo
    maximum = h.GetBinContent(h.GetMaximumBin())
    if x_range is not None:
        h.GetXaxis().SetRangeUser(*x_range) if x_range is not None else do_nothing()
        maximum = h.GetBinContent(h.GetMaximumBin())
        h.GetXaxis().UnZoom()
    value = maximum if to_max else value
    if value:
        h.Scale(1 / value)
    return h


def find_2d_centre(h, thresh=.5):
    px, py = h.ProjectionX(), h.ProjectionY()
    return array([mean([p.GetBinCenter(b) for b in [p.FindFirstBinAbove(p.GetMaximum() * thresh), p.FindLastBinAbove(p.GetMaximum() * thresh)]]) for p in [px, py]])


def get_2d_centre_ranges(h, dx, dy=None, thresh=.5):
    if dx is None:
        return None, None
    cx, cy = find_2d_centre(h, thresh)
    dx, dy = dx / 2, choose(dy, dx) / 2
    return [cx - dx, cx + dx], [cy - dy, cy + dy]


def centre_2d(h, dx, dy=None, thresh=.5):
    c = get_last_canvas()
    set_axes_range(*concatenate(get_2d_centre_ranges(h, dx, dy, thresh)), c)


def make_transparent(pad):
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.SetFrameFillStyle(4000)


def hide_axis(axis):
    axis.SetTickLength(0)
    axis.SetLabelOffset(99)
    axis.SetTitleOffset(99)


def set_root_warnings(status):
    gROOT.ProcessLine('gErrorIgnoreLevel = {e};'.format(e='0' if status else 'kError'))


def set_root_output(status=True):
    gROOT.SetBatch(not status)
    set_root_warnings(status)


if __name__ == '__main__':
    z = Draw(join(Draw.Dir, 'config', 'main.ini'))
