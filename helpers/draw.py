#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from os.path import join
from ROOT import TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TCanvas, gStyle, TLegend, TArrow, TPad, TCutG, TLine, TPaveText, TPaveStats, TH1F, TSpectrum, TEllipse, TColor, TProfile
from ROOT import TProfile2D, TH2F
from numpy import sign, linspace, ones, ceil
from src.binning import Bins
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
    Count = 0
    Res = load_resolution()
    Colors = get_color_gradient()
    Objects = []
    Dir = get_base_dir()

    Title = True
    Legend = False
    FillColor = 871
    Font = 42

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

    def __call__(self, *args, **kwargs):
        return Draw.histo(**kwargs)

    @staticmethod
    def add(*args):
        for obj in args:
            if obj not in Draw.Objects:
                Draw.Objects.append(obj)
        Draw.clean()
        return args[0] if len(args) == 1 else args

    @staticmethod
    def clean():
        n_none = sum(str(obj) == 'None' for obj in Draw.Objects)
        for _ in range(n_none):
            Draw.Objects.remove(None)

    # ----------------------------------------
    # region SET
    @staticmethod
    def setup():
        gStyle.SetLegendFont(Draw.Font)
        gStyle.SetOptTitle(Draw.Title)
        gStyle.SetPalette(Draw.Config.get_value('PLOTS', 'palette', default=1))
        gStyle.SetNumberContours(Draw.Config.get_value('PLOTS', 'contours', default=20))

    @staticmethod
    def set_pad_margins(c=None, l_=None, r=None, b=None, t=None):
        do(c.SetLeftMargin, l_)
        do(c.SetRightMargin, r if r is not None else None if round(c.GetRightMargin(), 1) != .1 else .03)
        do(c.SetBottomMargin, None if round(c.GetBottomMargin(), 1) != .1 else (.17 if b is None else b) - (.07 if not Draw.Legend else 0))
        do(c.SetTopMargin, None if round(c.GetTopMargin(), 1) != .1 else (.1 if t is None else t) - (0 if Draw.Title else .07))
    # endregion SET
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    @staticmethod
    def color(n, i):
        return Draw.get_colors(n)[i]

    def get_color(self, n, i=None):
        color = Draw.get_colors(n)[choose(i, self.IColor)]
        self.IColor = self.IColor + 1 if self.IColor < n - 1 else 0
        return color

    @staticmethod
    def get_colors(n):
        return Draw.Colors[linspace(0, Draw.Colors.size - 1, n).round().astype(int)].tolist()

    @staticmethod
    def get_count():
        Draw.Count += 1
        return Draw.Count
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
        c = TCanvas('c{}'.format(Draw.get_count()), title, int(x), int(y), int(w * Draw.Res), int(h * Draw.Res))
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx, gridy])
        do(make_transparent, c, transp)
        if divide is not None:
            c.Divide(*(divide if type(divide) in [list, tuple] else [divide]))
        return Draw.add(c)

    @staticmethod
    def axis(x1, x2, y1, y2, title, limits=None, name='ax', col=1, width=1, off=.15, tit_size=.035, lab_size=0.035, tick_size=0.03, line=False, opt='+SU', l_off=.01, log=False):
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
        a.SetTickSize(tick_size if not line else 0)
        a.SetTickLength(tick_size if not line else 0)
        a.SetNdivisions(0) if line else do_nothing()
        a.SetLabelOffset(l_off)
        a.Draw()
        return Draw.add(a)

    @staticmethod
    def y_axis(x, ymin, ymax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return Draw.axis(x, x, ymin, ymax, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    @staticmethod
    def x_axis(y, xmin, xmax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return Draw.axis(xmin, xmax, y, y, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    @staticmethod
    def line(x1, x2, y1, y2, color=1, width=1, style=1, name='li'):
        line = TCutG(name, 2, array([x1, x2], 'd'), array([y1, y2], 'd'))
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw('same')
        return Draw.add(line)

    @staticmethod
    def tline(x1, x2, y1, y2, color=1, width=1, style=1):
        line = TLine(x1, y1, x2, y2)
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw()
        return Draw.add(line)

    @staticmethod
    def vertical_line(x, ymin, ymax, color=1, w=1, style=1, name='li', tline=False):
        return Draw.line(x, x, ymin, ymax, color, w, style, name) if not tline else Draw.tline(x, x, ymin, ymax, color, w, style)

    @staticmethod
    def horizontal_line(y, xmin, xmax, color=1, w=1, style=1, name='li', tline=False):
        return Draw.line(xmin, xmax, y, y, color, w, style, name) if not tline else Draw.tline(xmin, xmax, y, y, color, w, style)

    @staticmethod
    def polygon(x, y, line_color=1, width=1, style=1, name=None, fillstyle=None, fill_color=None, show=True):
        line = TCutG(choose(name, 'poly{}'.format(Draw.get_count())), len(x) + 1, array(x + [x[0]], 'd'), array(y + [y[0]], 'd'))
        line.SetLineColor(line_color)
        do(line.SetFillColor, fill_color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        do(line.SetFillStyle, fillstyle)
        if show:
            line.Draw('l')
            line.Draw('f') if fill_color is not None or fillstyle is not None and fillstyle < 4000 else do_nothing()
        return Draw.add(line)

    @staticmethod
    def box(x1, y1, x2, y2, line_color=1, width=1, style=1, name=None, fillstyle=None, fillcolor=None, show=True):
        return Draw.polygon([x1, x1, x2, x2], [y1, y2, y2, y1], line_color, width, style, name, fillstyle, fillcolor, show)

    @staticmethod
    def tlatex(x, y, text, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None, show=True):
        tlatex = TLatex(x, y, text)
        format_text(tlatex, name, align, color, size, angle, ndc, font)
        tlatex.Draw() if show else do_nothing()
        return Draw.add(tlatex)

    @staticmethod
    def arrow(x1, x2, y1, y2, col=1, width=1, opt='<|', size=.005):
        ar = TArrow(x1, y1, x2, y2, size, opt)
        ar.SetLineWidth(width)
        ar.SetLineColor(col)
        ar.SetFillColor(col)
        ar.Draw()
        return Draw.add(ar)

    @staticmethod
    def tpad(name, tit='', pos=None, fill_col=0, gridx=None, gridy=None, margins=None, transparent=False, logy=None, logx=None, logz=None, lm=None, rm=None, bm=None, tm=None):
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(name, tit, *pos)
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
            Draw.line(x, x, min(y_vals), max(y_vals), name='x{}'.format(x), width=width, color=color)
        for y in y_vals:
            Draw.line(min(x_vals), max(x_vals), y, y, name='y{}'.format(y), width=width, color=color)

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
    def histo(th, show=True, lm=None, rm=None, bm=None, tm=None, m=None, draw_opt=None, w=1, h=1, logy=False, logx=False, logz=False, grid=False, gridy=False, gridx=False, phi=None, theta=None,
              leg=None, canvas=None, sumw2=None):
        w += .16 if not Draw.Title and w == 1 else 0  # rectify if there is no title
        th.Sumw2(sumw2) if hasattr(th, 'Sumw2') and sumw2 is not None else do_nothing()
        set_root_output(show)
        c = Draw.canvas(th.GetTitle().split(';')[0], None, None, w, h, logx, logy, logz, gridx or grid, gridy or grid, show=show) if canvas is None else canvas
        Draw.set_pad_margins(c, *[lm, rm, bm, tm] if m is None else m)
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx or grid, gridy or grid])
        do([c.SetPhi, c.SetTheta], [phi, theta])
        th.Draw(draw_opt if draw_opt is not None else 'ap' if is_graph(th) else '')
        if leg is not None:
            for i_leg in make_list(leg):
                i_leg.Draw()
        set_root_output(True)
        return Draw.add(c, th, leg)[0]

    # TODO: revise and put at correct place
    def save_combined_pulse_heights(self, mg, mg1, mg_y, show=True, name=None, pulser_leg=None,
                                    x_range=None, y_range=None, rel_y_range=None, draw_objects=None, prnt=True):
        set_root_output(show)
        c = TCanvas('c', 'c', int(Draw.Res * 10 / 11.), Draw.Res)
        make_transparent(c)
        bm = .11
        scale = 1.5
        pm = bm + (1 - bm - .1) / 5.

        # set unified x-range:
        mgn = mg1.Clone()
        mgn.GetXaxis().SetLimits(1, 3e4) if x_range is None else do_nothing()
        mg.GetXaxis().SetLimits(1, 3e4) if x_range is None else do_nothing()

        # bottom pad with 20%
        p0 = self.tpad('p0', 'p0', pos=[0, 0, 1, pm], margins=[.14, .03, bm / pm, 0], transparent=True, logx=True, gridy=True)
        scale_multigraph(mgn)
        rel_y_range = [.7, 1.3] if rel_y_range is None else rel_y_range
        format_histo(mgn, title='', y_range=rel_y_range, y_tit='Rel. ph [au]' if not scale > 1 else ' ', y_off=66, tit_size=.1 * scale, x_off=99, lab_size=.1 * scale)
        mgn.GetYaxis().SetNdivisions(3)
        hide_axis(mgn.GetXaxis())
        mgn.Draw('alp')
        x_range = [mgn.GetXaxis().GetXmin(), mgn.GetXaxis().GetXmax()] if x_range is None else x_range
        self.x_axis(1.3, x_range[0], x_range[1], mgn.GetXaxis().GetTitle() + ' ', opt='SG+-=', tit_size=.1, lab_size=.1 * scale, off=99, tick_size=.1, l_off=0)
        c.cd()

        # top pad with zero suppression
        self.tpad('p1', 'p1', pos=[0, pm, 1, 1], margins=[.14, .03, 0, .1], transparent=True, logx=True)
        mg.Draw('alp')
        hide_axis(mg.GetXaxis())
        if pulser_leg:
            pulser_leg()
        if y_range:
            mg.SetMinimum(y_range[0])
            mg.SetMaximum(y_range[1])
        format_histo(mg, tit_size=.04 * scale, y_off=1.75 / scale, lab_size=.04 * scale)
        self.x_axis(mg_y, x_range[0], x_range[1], mgn.GetXaxis().GetTitle() + ' ', opt='SG=', tit_size=.035 * scale, lab_size=0, off=1, l_off=99)
        leg = mg.GetListOfFunctions()[0]
        move_legend(leg, .17, .03)
        leg.Draw()
        if draw_objects is not None:
            for obj, opt in draw_objects:
                obj.Draw(opt)

        if hasattr(self, 'InfoLegend'):
            run_info = self.InfoLegend.draw(p0, all_pads=False)
            scale_legend(run_info[0], txt_size=.09, height=0.098 / pm)
            run_info[1].SetTextSize(.05)

        for obj in p0.GetListOfPrimitives():
            if obj.GetName() == 'title':
                obj.SetTextColor(0)
        self.save_canvas(c, name='CombinedPulseHeights' if name is None else name, show=show, print_names=prnt)

        Draw.add(c, *draw_objects)
        set_root_output(True)

    def distribution(self, values, binning=None, title='', thresh=.02, lm=None, rm=None, show=True, logy=None, **kwargs):
        values = array(values, dtype='d')
        kwargs['fill_color'] = Draw.FillColor if 'fill_color' not in kwargs else kwargs['fill_color']
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        kwargs['y_tit'] = 'Number of Entries' if 'y_tit' not in kwargs else kwargs['y_tit']
        h = TH1F('h{}'.format(Draw.get_count()), title, *choose(binning, make_bins, values=values, thresh=thresh))
        fill_hist(h, values)
        format_histo(h, **kwargs)
        self.histo(h, show=show, lm=lm, rm=rm, logy=logy)
        return h

    def graph(self, x, y, ex=None, ey=None, asym_errors=False, lm=None, rm=None, tm=None, w=1, h=1, show=True, draw_opt=None, logy=False, **kwargs):
        g = Draw.make_tgrapherrors(x, y, ex, ey, asym_errors)
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        format_histo(g, **kwargs)
        self.histo(g, show=show, lm=lm, rm=rm, tm=tm, w=w, h=h, draw_opt=draw_opt, logy=logy)
        return g

    def profile(self, x, y, binning=None, title='', thresh=.02, lm=None, rm=None, w=1, h=1, show=True, draw_opt=None, **kwargs):
        x, y = array(x, dtype='d'), array(y, dtype='d')
        kwargs['fill_color'] = Draw.FillColor if 'fill_color' not in kwargs else kwargs['fill_color']
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        p = TProfile('p{}'.format(Draw.get_count()), title, *choose(binning, make_bins, values=x, thresh=thresh))
        fill_hist(p, x, y)
        format_histo(p, **kwargs)
        self.histo(p, show=show, lm=lm, rm=rm, w=w, h=h, draw_opt=draw_opt)
        return p

    def prof2d(self, x, y, zz, binning=None, title='', lm=None, rm=.15, w=1, h=1, show=True, draw_opt='colz', **kwargs):
        x, y, zz = array(x, dtype='d'), array(y, dtype='d'), array(zz, dtype='d')
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        kwargs['z_off'] = 1.2 if 'z_off' not in kwargs else kwargs['z_off']
        dflt_bins = Bins.make(min(x), max(x), sqrt(x.size)) + Bins.make(min(y), max(y), sqrt(x.size))
        p = TProfile2D('p{}'.format(self.get_count()), title, *choose(binning, dflt_bins))
        fill_hist(p, x, y, zz)
        format_histo(p, pal=55, **kwargs)
        self.histo(p, show=show, lm=lm, rm=rm, w=w, h=h, draw_opt=draw_opt)
        return p

    def histo_2d(self, x, y, binning=None, title='', lm=None, rm=.15, show=True, draw_opt='colz', **kwargs):
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        kwargs['z_off'] = 1.2 if 'z_off' not in kwargs else kwargs['z_off']
        kwargs['z_tit'] = 'Number of Entries' if 'z_tit' not in kwargs else kwargs['z_tit']
        x, y = array(x, dtype='d'), array(y, dtype='d')
        dflt_bins = Bins.make(min(x), max(x), sqrt(x.size)) + Bins.make(min(y), max(y), sqrt(x.size))
        h = TH2F('h{}'.format(self.get_count()), title, *choose(binning, dflt_bins))
        fill_hist(h, x, y)
        format_histo(h, **kwargs)
        self.histo(h, show=show, lm=lm, rm=rm, draw_opt=draw_opt)
        return h

    def efficiency(self, x, e, binning=None, title='', lm=None, show=True, **kwargs):
        binning = choose(binning, Bins.make, min(x), max(x), (max(x) - min(x)) / sqrt(x.size))
        p = self.profile(x, e, binning, show=False)
        x = get_hist_args(p, err=False)
        values = [[p.GetBinContent(ibin), p.GetBinEntries(ibin)] for ibin in range(1, p.GetNbinsX() + 1)]
        e = array([calc_eff(p0 / 100 * n, n) for p0, n in values])
        ey = array([e[:, 1], e[:, 2]])
        g = Draw.make_tgrapherrors(x=x, y=e[:, 0], ey=ey, asym_err=True)
        format_histo(g, title=title, **kwargs)
        self.histo(g, show=show, lm=lm, draw_opt='ap')
        return g
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region CREATE
    @staticmethod
    def make_histo(title, *bins):
        h = TH1F('h{}'.format(Draw.get_count()), title, *bins)
        return Draw.add(h)

    @staticmethod
    def make_tgrapherrors(x=None, y=None, ex=None, ey=None, asym_err=False, **kwargs):
        g = (TGraphAsymmErrors if asym_err else TGraphErrors)(*make_graph_args(x, y, ex, ey, asym_err))
        kwargs['marker'] = 20 if 'marker' not in kwargs else kwargs['marker']
        kwargs['markersize'] = 1 if 'markersize' not in kwargs else kwargs['markersize']
        format_histo(g, 'g{}'.format(Draw.get_count()), **kwargs)
        return Draw.add(g)

    @staticmethod
    def make_graph_from_profile(p):
        x_range = [i for i in range(p.GetNbinsX()) if p.GetBinContent(i)]
        x = [make_ufloat([p.GetBinCenter(i), p.GetBinWidth(i) / 2]) for i in x_range]
        y = [make_ufloat([p.GetBinContent(i), p.GetBinError(i)]) for i in x_range]
        return Draw.make_tgrapherrors(x, y)

    @staticmethod
    def make_legend(x1=.65, y2=.88, nentries=2, scale=1, name='l', y1=None, clean=False, margin=.25, x2=None, w=None, cols=None, fix=False):
        x2 = .95 if x2 is None else x2
        x1 = x2 - w if w is not None else x1
        h = nentries * .05 * scale
        y = array([y2 - h if y1 is None else y1, y1 + h if y1 is not None else y2])
        y += .07 if not Draw.Title and y[1] > .7 and not fix else 0
        y -= .07 if not Draw.Legend and y[1] < .7 and not fix else 0
        leg = TLegend(x1, max(y[0], 0), x2, min(y[1], 1))
        leg.SetName(name)
        leg.SetTextFont(42)
        leg.SetTextSize(0.03 * scale)
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
def format_histo(histo, name=None, title=None, x_tit=None, y_tit=None, z_tit=None, marker=20, color=None, line_color=None, markersize=None, x_off=None, y_off=None, z_off=None, lw=1,
                 fill_color=None, fill_style=None, stats=True, tit_size=None, lab_size=None, l_off_y=None, l_off_x=None, draw_first=False, x_range=None, y_range=None, z_range=None, sumw2=None,
                 do_marker=True, style=None, ndivx=None, ndivy=None, ncont=None, tick_size=None, t_ax_off=None, center_y=False, center_x=False, yax_col=None, normalise=None, pal=None, rebin=None,
                 y_ticks=None, x_ticks=None, z_ticks=None, opacity=None):
    h = histo
    if draw_first:
        set_root_output(False)
        h.Draw('nostack' if h.ClassName() == 'THStack' else 'a')
        set_root_output(True)
    do(h.SetTitle, title)
    do(h.SetName, name)
    do(set_palette, pal)
    if normalise is not None:
        y_tit = y_tit.replace('Number', 'Percentage') if y_tit is not None else y_tit
        h.Sumw2(True)
        normalise_histo(h)
    try:
        h.SetStats(stats)
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
        h.SetLineWidth(lw)
        h.SetFillColor(fill_color) if fill_color is not None and opacity is None else do_nothing()
        h.SetFillColorAlpha(fill_color, opacity) if fill_color is not None and opacity is not None else do_nothing()
        h.SetFillStyle(fill_style) if fill_style is not None else do_nothing()
        h.SetFillStyle(style) if style is not None else do_nothing()
        h.SetContour(ncont) if ncont is not None else do_nothing()
    except AttributeError or ReferenceError:
        pass
    # axes
    try:
        x_args = [x_tit, x_off, tit_size, center_x, lab_size, l_off_x, x_range, ndivx, choose(x_ticks, tick_size), ]
        y_args = [y_tit, y_off, tit_size, center_y, lab_size, l_off_y, y_range, ndivy, choose(y_ticks, tick_size), yax_col]
        z_args = [z_tit, z_off, tit_size, False, lab_size, None, z_range, None, choose(z_ticks, tick_size)]
        for i, name in enumerate(['X', 'Y', 'Z']):
            format_axis(getattr(h, 'Get{}axis'.format(name))(), h, *[x_args, y_args, z_args][i])
    except AttributeError or ReferenceError:
        pass
    set_time_axis(h, off=t_ax_off) if t_ax_off is not None else do_nothing()
    do(h.Sumw2, sumw2) if hasattr(h, 'Sumw2') else do_nothing()
    update_canvas()


def format_statbox(x=.95, y=None, w=.2, h=.15, only_fit=False, fit=False, entries=False, form=None, m=False, rms=False, all_stat=False):
    gStyle.SetOptFit(int(only_fit or fit))
    opt_stat = '100000{}{}{}0'.format(*ones(3, 'i') if all_stat else array([rms, m, entries]).astype('i'))
    if only_fit:
        opt_stat = '0011'
    y = (.88 if Draw.Title else .95) if y is None else y
    gStyle.SetOptStat(int(opt_stat))
    gStyle.SetFitFormat(form) if form is not None else do_nothing()
    gStyle.SetStatX(x)
    gStyle.SetStatY(y)
    gStyle.SetStatW(w)
    gStyle.SetStatH(h)


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


def adapt_z_range(h, n_sigma=2):
    values = get_2d_hist_vec(h)
    m, s = mean_sigma(values[5:-5], err=False)
    z_range = [min(values).n, .8 * max(values).n] if s > m else [m - n_sigma * s, m + n_sigma * s]
    format_histo(h, z_range=z_range)


def make_bins(values, thresh=.02):
    binning = linspace(*(find_range(values, thresh=thresh) + [int(sqrt(values.size))]))
    return [binning.size - 1, binning]


def find_range(values, lfac=.2, rfac=.2, thresh=.02):
    v = array(sorted(values))
    xmin, xmax = v[int(thresh * v.size)], v[int(v.size - thresh * v.size)]
    return ax_range(xmin, xmax, lfac, rfac)


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
    return FitRes(fit) if fit is not None else FitRes()


def make_graph_args(x, y, ex=None, ey=None, asym_errors=False):
    if x is None:
        return []
    if len(list(x)) != len(list(y)):
        warning('Arrays have different size!')
        return []
    s = len(x)
    utypes = [Variable, AffineScalarFunc]
    ex, ey = [v if type(v) in [list, ndarray] or v is None else full(s, v) for v in [ex, ey]]
    ex, ey = [array([v.s for v in vals], 'd') if type(vals[0]) in utypes else zeros((2, s) if asym_errors else s) if ers is None else array(ers, 'd') for vals, ers in zip([x, y], [ex, ey])]
    x, y = [array([v.n for v in vals] if type(vals[0]) in utypes else vals, 'd') for vals in [x, y]]
    return [s, array(x, 'd'), array(y, 'd')] + ([ex[0], ex[1], ey[0], ey[1]] if asym_errors else [ex, ey])


def set_titles(status=True):
    gStyle.SetOptTitle(status)


def get_graph_vecs(g, err=True):
    return get_graph_x(g, err), get_graph_y(g, err)


def get_graph_x(g, err=True):
    values = array([make_ufloat([g.GetX()[i], g.GetEX()[i]]) for i in range(g.GetN())]) if 'Error' in g.ClassName() else array([make_ufloat(g.GetX()[i]) for i in range(g.GetN())])
    return values if err else array([v.n for v in values])


def get_graph_y(g, err=True):
    values = array([make_ufloat([g.GetY()[i], g.GetEY()[i]]) for i in range(g.GetN())]) if 'Error' in g.ClassName() else array([make_ufloat(g.GetY()[i]) for i in range(g.GetN())])
    return values if err else array([v.n for v in values])


def get_hist_vec(p, err=True):
    return array([make_ufloat([p.GetBinContent(ibin), p.GetBinError(ibin)]) if err else p.GetBinContent(ibin) for ibin in range(1, p.GetNbinsX() + 1)])


def get_hist_args(p, err=True):
    return array([make_ufloat([p.GetBinCenter(ibin), p.GetBinWidth(ibin) / 2]) if err else p.GetBinCenter(ibin) for ibin in range(1, p.GetNbinsX() + 1)])


def get_hist_vecs(p, err=True):
    return get_hist_args(p, err), get_hist_vec(p, err)


def get_h_values(h):
    return get_graph_y(h) if 'Graph' in h.ClassName() else get_hist_vec(h)


def get_h_args(h):
    return get_graph_x(h) if 'Graph' in h.ClassName() else get_hist_args(h)


def get_2d_hist_vec(h, err=True, flat=True):
    xbins, ybins = range(1, h.GetNbinsX() + 1), range(1, h.GetNbinsY() + 1)
    values = array([ufloat(h.GetBinContent(xbin, ybin), h.GetBinError(xbin, ybin)) for xbin in xbins for ybin in ybins])
    values = values if err else array([v.n for v in values])
    return values[values != 0] if flat else values.reshape(len(xbins), len(ybins))


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


def get_pull(h, name, binning, fit=True):
    set_root_output(False)
    h_out = TH1F('hp{}'.format(name[:3]), name, *binning)
    values = array([h.GetBinContent(ibin + 1) for ibin in range(h.GetNbinsX())], 'd')
    h_out.FillN(values.size, values, full(values.size, 1, 'd'))
    h_out.Fit('gaus', 'q') if fit else do_nothing()
    format_histo(h_out, x_range=ax_range(values.min(), values.max(), .1, .3))
    return h_out


def get_quantile(h, q):
    quantiles = make_list(q)
    v = zeros(quantiles.size)
    h.GetQuantiles(v.size, v, quantiles)
    return v[0] if v.size == 1 else v


def markers(i):
    return ((list(range(20, 24)) + [29, 33, 34]) * 2)[i]


def fit_bucket(histo, show=True):
    # TODO move to appropriate spot and revise
    set_root_output(False)
    h = histo
    format_histo(h, rebin=int(h.GetBinCenter(h.FindLastBinAbove(h.GetMaximum() * .02))) // 40)
    fit = TF1('fit', 'gaus(0) + gaus(3) + gaus(6)', h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    s = TSpectrum(3)
    n = s.Search(h, 2.5)
    points = [(s.GetPositionX()[i], s.GetPositionY()[i]) for i in [0, 1 if n == 2 else 2]]
    x1, x2 = (p[0] for p in sorted(points))
    y1, y2 = (p[1] for p in sorted(points))
    if y1 < 20 or y1 > 1e10:
        return  # didn't find pedestal peak!
    diff = x2 - x1
    fit.SetParameters(*[y2, x2, 10, y1, x1, 3, min(y1, y2) / 4, x1 + diff / 4, 5])
    # signal
    fit.SetParLimits(1, x2 - 5, x2 + 5)
    # pedestal
    fit.SetParLimits(3, 1, y1 * 2)
    fit.SetParLimits(4, x1 - 10, x1 + 10)
    # middle ped
    fit.SetParLimits(6, 1, min(y1, y2) / 2)
    fit.SetParLimits(7, x1, x1 + diff / 2)
    for i in range(1):
        h.Fit(fit, 'qs{0}'.format('' if show else '0'), '', -50, x2 + 5)
    set_root_warnings(1)
    return fit


def set_palette(pal):
    gStyle.SetPalette(pal)


def is_graph(h):
    return 'Graph' in h.ClassName()


def update_canvas(c=None):
    c = choose(c, get_last_canvas(warn=False))
    if c is not None:
        c.Modified()
        c.Update()


def show_colors(colors):
    n = len(colors)
    c = Draw.canvas(divide=(int(ceil(sqrt(n))), int(ceil(sqrt(n)))))
    for i, col in enumerate(colors, 1):
        c.cd(i)
        Draw.box(0, 0, 1, 1, fillstyle=1001, fillcolor=col)
        Draw.tlatex(.5, .5, str(i - 1), align=22, size=.2)


def ax_range(low, high=None, fl=0, fh=0, h=None):
    if h is not None:
        if 'TH2' in h.ClassName() or '2D' in h.ClassName():
            return [ax_range(axis.GetBinCenter(h.FindFirstBinAbove(low, i)), axis.GetBinCenter(h.FindLastBinAbove(high, i)), fl, fh) for i, axis in enumerate([h.GetXaxis(), h.GetYaxis()], 1)]
        return ax_range(h.GetBinCenter(h.FindFirstBinAbove(low)), h.GetBinCenter(h.FindLastBinAbove(high)), fl, fh)
    low, high = low if high is None else (low, high)
    d = abs(high - low)
    return [low - d * fl, high + d * fh]


if __name__ == '__main__':
    z = Draw(join(Draw.Dir, 'config', 'main.ini'))
