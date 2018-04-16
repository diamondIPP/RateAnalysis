#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gROOT, TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TCanvas, gStyle, TLegend, TArrow, TPad, TCutG, TLine, kGreen, kOrange, kViolet, kYellow, kRed, kBlue, kMagenta, kAzure, \
    kCyan, kTeal
from Utils import *
from os.path import dirname
from numpy import ndarray, zeros, sign


class Draw:

    def __init__(self, elementary):

        self.Config = elementary.MainConfigParser
        self.Dir = elementary.Dir
        self.TESTCAMPAIGN = elementary.TESTCAMPAIGN
        self.Verbose = elementary.verbose
        self.TCString = elementary.TCString
        self.Res = elementary.Res
        self.ResultsDir = self.generate_results_directory()

        self.Objects = []

        # colors
        self.count = 0
        self.colors = create_colorlist()
        self.FillColor = 821
        gStyle.SetLegendFont(42)

        self.set_titles()

    def set_titles(self, on=None):
        gStyle.SetOptTitle(self.Config.getboolean('SAVE', 'activate_title') if on is None else on)

    def set_save_directory(self, name):
        self.ResultsDir = join(self.Dir, name)

    def generate_results_directory(self):
        return join(self.Dir, 'Results{tc}'.format(tc=self.TCString))

    def make_bias_string(self, bias=None):
        if bias is None:
            return self.make_bias_string(self.bias) if hasattr(self, 'bias') else ''
        pol = 'm' if bias < 0 else 'p'
        return '_{pol}{bias:04d}'.format(pol=pol, bias=int(abs(bias)))

    def make_info_string(self):
        info = ''
        if not self.Config.getboolean('SAVE', 'short_name'):
            info = '_{dia}'.format(dia=self.diamond_name) if hasattr(self, 'diamond_name') else ''
            info += self.make_bias_string()
            info += '_{tc}'.format(tc=self.TESTCAMPAIGN)
            info = info.replace('-', '')
        return info

    def get_color(self):
        self.count %= 20
        color = self.colors[self.count]
        self.count += 1
        return color

    def reset_colors(self):
        self.count = 0
    
    # region DRAWING

    def draw_axis(self, x1, x2, y1, y2, title, name='ax', col=1, width=1, off=.15, tit_size=.035, lab_size=0.035, line=False, opt='+SU', tick_size=0.03, l_off=.01):
        range_ = [y1, y2] if x1 == x2 else [x1, x2]
        a = TGaxis(x1, y1, x2, y2, range_[0], range_[1], 510, opt)
        a.SetName(name)
        a.SetLineColor(col)
        a.SetLineWidth(width)
        a.SetLabelSize(lab_size if not line else 0)
        a.SetTitleSize(tit_size)
        a.SetTitleOffset(off)
        a.SetTitle(title)
        a.SetTitleColor(col)
        a.SetLabelColor(col)
        a.SetLabelFont(42)
        a.SetTitleFont(42)
        a.SetTickSize(tick_size if not line else 0)
        a.SetTickLength(tick_size if not line else 0)
        a.SetNdivisions(0) if line else do_nothing()
        a.SetLabelOffset(l_off)
        a.Draw()
        self.Objects.append(a)
        return a

    def draw_y_axis(self, x, ymin, ymax, tit, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False):
        return self.draw_axis(x, x, ymin, ymax, tit, name=name, col=col, off=off, opt=opt, width=w, tit_size=tit_size, lab_size=lab_size, tick_size=tick_size, l_off=l_off, line=line)

    def draw_x_axis(self, y, xmin, xmax, tit, col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False):
        return self.draw_axis(xmin, xmax, y, y, tit, col=col, off=off, opt=opt, width=w, tit_size=tit_size, lab_size=lab_size, tick_size=tick_size, l_off=l_off, line=line)

    def draw_line(self, x1, x2, y1, y2, color=1, width=1, style=1, name='li'):
        l = TCutG(name, 2, array([x1, x2], 'd'), array([y1, y2], 'd'))
        l.SetLineColor(color)
        l.SetLineWidth(width)
        l.SetLineStyle(style)
        l.Draw('same')
        self.Objects.append(l)
        return l

    def draw_tline(self, x1, x2, y1, y2, color=1, width=1, style=1):
        l = TLine(x1, y1, x2, y2)
        l.SetLineColor(color)
        l.SetLineWidth(width)
        l.SetLineStyle(style)
        l.Draw()
        self.Objects.append(l)
        return l

    def draw_box(self, x1, y1, x2, y2, color=1, width=1, style=1, fillstyle=None, name='box', show=True):
        l = TCutG(name, 5, array([x1, x1, x2, x2, x1], 'd'), array([y1, y2, y2, y1, y1], 'd'))
        l.SetLineColor(color)
        l.SetFillColor(color)
        l.SetLineWidth(width)
        l.SetLineStyle(style)
        l.SetFillStyle(fillstyle) if fillstyle is not None else do_nothing()
        if show:
            l.Draw('same')
        self.Objects.append(l)
        return l

    def draw_vertical_line(self, x, ymin, ymax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(x, x, ymin, ymax, color, w, style, name) if not tline else self.draw_tline(x, x, ymin, ymax, color, w, style)

    def draw_horizontal_line(self, y, xmin, xmax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(xmin, xmax, y, y, color, w, style, name) if not tline else self.draw_tline(xmin, xmax, y, y, color, w, style)

    def draw_tlatex(self, x, y, text, align=20, color=1, size=.05, angle=None, ndc=False):
        l = TLatex(x, y, text)
        l.SetName(text)
        l.SetTextAlign(align)
        l.SetTextColor(color)
        l.SetTextSize(size)
        do(l.SetTextAngle, angle)
        l.SetNDC() if ndc else do_nothing()
        l.Draw()
        self.Objects.append(l)
        return l

    def draw_arrow(self, x1, x2, y1, y2, col=1, width=1, opt='<|', size=.005):
        ar = TArrow(x1, y1, x2, y2, size, opt)
        ar.SetLineWidth(width)
        ar.SetLineColor(col)
        ar.SetFillColor(col)
        ar.Draw()
        self.Objects.append(ar)

    def draw_tpad(self, name, tit='', pos=None, fill_col=0, gridx=None, gridy=None, margins=None, transparent=False, logy=None, logx=None, logz=None):
        margins = [.1, .1, .1, .1] if margins is None else margins
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(name, tit, *pos)
        p.SetFillColor(fill_col)
        p.SetMargin(*margins)
        do([p.SetLogx, p.SetLogy, p.SetLogz], [logx, logy, logz])
        do([p.SetGridx, p.SetGridy], [gridx, gridy])
        make_transparent(p) if transparent else do_nothing()
        p.Draw()
        p.cd()
        self.Objects.append(p)
        return p

    def draw_preliminary(self, canvas=None):
        c = get_last_canvas() if canvas is None else canvas
        self.draw_tlatex((1 - c.GetRightMargin() + c.GetLeftMargin()) / 2, .5, 'Preliminary', align=22, color=19, size=.17, angle=30, ndc=True)
        make_transparent(c)
        for obj in c.GetListOfPrimitives():
            if obj.IsA().GetName() in ['TH1F', 'TH2F', 'TGraph', 'TGraphErrors']:
                obj.Draw('same')
        c.RedrawAxis()

    def draw_histo(self, histo, save_name='', show=True, sub_dir=None, lm=.1, rm=.03, bm=.15, tm=.1, draw_opt='', x=None, y=None, all_pads=True,
                   l=None, logy=False, logx=False, logz=False, canvas=None, gridy=False, gridx=False, both_dias=False, prnt=True, phi=None, theta=None, ind=None):
        return self.save_histo(histo, save_name, show, sub_dir, lm, rm, bm, tm, draw_opt, x, y, all_pads, l, logy, logx, logz, canvas, gridx, gridy, False, both_dias, ind, prnt, phi, theta)
    
    # endregion

    # region SAVE
    def save_on_server(self, canvas, file_name):
        if server_is_mounted():
            if hasattr(self, 'DiamondName'):
                if hasattr(self, 'RunPlan'):
                    rp = self.RunPlan
                    run_string = 'RunPlan{r}'.format(r=rp[1:] if rp[0] == '0' else rp)
                elif hasattr(self, 'RunNumber'):
                    run_string = str(self.RunNumber)
                else:
                    return
                path = join(get_base_dir(), 'mounts/psi/Diamonds', self.DiamondName, 'BeamTests', make_tc_str(self.TCString, long_=False), run_string, file_name)
                canvas.SaveAs('{p}.pdf'.format(p=path))
                canvas.SaveAs('{p}.png'.format(p=path))

    def save_plots(self, savename, sub_dir=None, canvas=None, all_pads=True, both_dias=False, ind=None, prnt=True, save=True, show=True):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        canvas = get_last_canvas() if canvas is None else canvas
        if canvas is None:
            return
        if ind is None:
            self.InfoLegend.draw(canvas, all_pads, both_dias, show) if hasattr(self, 'InfoLegend') else log_warning('Did not find InfoLegend class...') \
                if not any(hasattr(self, att) for att in ['RunSelections', 'CurrentGraph']) else do_nothing()
        else:
            self.collection.values()[ind].InfoLegend.draw(canvas, all_pads, both_dias, show) if hasattr(self, 'collection') else log_critical('sth went wrong...')
        canvas.Modified()
        canvas.Update()
        if save:
            try:
                self.save_canvas(canvas, sub_dir=sub_dir, name=savename, print_names=prnt, show=show)
                self.Objects.append(canvas)
            except Exception as inst:
                log_warning('Error in save_canvas:\n{0}'.format(inst))

    def save_canvas(self, canvas, sub_dir=None, name=None, print_names=True, show=True):
        sub_dir = self.save_dir if hasattr(self, 'save_dir') and sub_dir is None else sub_dir
        canvas.Update()
        file_name = canvas.GetName() if name is None else name
        file_path = join(self.ResultsDir, sub_dir, '{typ}', file_name)
        ftypes = ['root', 'png', 'pdf', 'eps']
        out = 'Saving plots: {nam}'.format(nam=name)
        run_number = self.run_number if hasattr(self, 'run_number') else None
        run_number = 'rp{nr}'.format(nr=self.run_plan) if hasattr(self, 'run_plan') else run_number
        set_root_output(show)
        gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
        info = self.make_info_string()
        for f in ftypes:
            ext = '.{typ}'.format(typ=f)
            if not f == 'png' and run_number is not None:
                ext = '{str}_{run}.{typ}'.format(str=info, run=run_number, typ=f)
            ensure_dir(dirname(file_path.format(typ=f)))
            out_file = '{fname}{ext}'.format(fname=file_path, ext=ext)
            out_file = out_file.format(typ=f)
            canvas.SaveAs(out_file)
        self.save_on_server(canvas, file_name)
        if print_names:
            log_message(out, prnt=self.Verbose)
        set_root_output(True)

    def save_histo(self, histo, save_name='test', show=True, sub_dir=None, lm=.1, rm=.03, bm=.15, tm=None, draw_opt='', x_fac=None, y_fac=None, all_pads=True,
                   l=None, logy=False, logx=False, logz=False, canvas=None, gridx=False, gridy=False, save=True, both_dias=False, ind=None, prnt=True, phi=None, theta=None):
        if tm is None:
            tm = .1 if self.Config.getboolean('SAVE', 'activate_title') else .03
        x = self.Res if x_fac is None else int(x_fac * self.Res)
        y = self.Res if y_fac is None else int(y_fac * self.Res)
        h = histo
        set_root_output(show)
        c = TCanvas('c_{0}'.format(h.GetName()), h.GetTitle().split(';')[0], x, y) if canvas is None else canvas
        c.SetMargin(lm, rm, bm, tm)
        c.SetLogx() if logx else do_nothing()
        c.SetLogy() if logy else do_nothing()
        c.SetLogz() if logz else do_nothing()
        c.SetGridx() if gridx else do_nothing()
        c.SetGridy() if gridy else do_nothing()
        c.SetPhi(phi) if phi is not None else do_nothing()
        c.SetTheta(theta) if theta is not None else do_nothing()
        h.Draw(draw_opt)
        if l is not None:
            l = [l] if type(l) is not list else l
            for i in l:
                i.Draw()
        self.save_plots(save_name, sub_dir=sub_dir, both_dias=both_dias, all_pads=all_pads, ind=ind, prnt=prnt, save=save, show=show)
        set_root_output(True)
        lst = [c, h, l] if l is not None else [c, h]
        self.Objects.append(lst)
        return c

    # endregion

    @staticmethod
    def format_histo(histo, name='', title='', x_tit='', y_tit='', z_tit='', marker=20, color=1, markersize=1, x_off=None, y_off=None, z_off=None, lw=1,
                     fill_color=None, fill_style=None, stats=True, tit_size=.04, lab_size=.04, l_off_y=None, draw_first=False, x_range=None, y_range=None, z_range=None,
                     do_marker=True, style=None, ndivx=None, ndivy=None, ncont=None, tick_size=None, t_ax_off=None):
        h = histo
        if draw_first:
            set_root_output(False)
            h.Draw('a')
            set_root_output(True)
        h.SetTitle(title) if title else h.SetTitle(h.GetTitle())
        h.SetName(name) if name else h.SetName(h.GetName())
        try:
            h.SetStats(stats)
        except AttributeError or ReferenceError:
            pass
        # markers
        try:
            if do_marker:
                h.SetMarkerStyle(marker) if marker is not None else do_nothing()
                h.SetMarkerColor(color) if color is not None else do_nothing()
                h.SetMarkerSize(markersize) if markersize is not None else do_nothing()
        except AttributeError or ReferenceError:
            pass
        # lines/fill
        try:
            h.SetLineColor(color) if color is not None else h.SetLineColor(h.GetLineColor())
            h.SetFillColor(fill_color) if fill_color is not None else do_nothing()
            h.SetFillStyle(fill_style) if fill_style is not None else do_nothing()
            h.SetLineWidth(lw)
            h.SetFillStyle(style) if style is not None else do_nothing()
            h.SetContour(ncont) if ncont is not None else do_nothing()
        except AttributeError or ReferenceError:
            pass
        # axis titles
        try:
            x_axis = h.GetXaxis()
            if x_axis:
                x_axis.SetTitle(x_tit) if x_tit else h.GetXaxis().GetTitle()
                x_axis.SetTitleOffset(x_off) if x_off is not None else do_nothing()
                x_axis.SetTitleSize(tit_size)
                x_axis.SetLabelSize(lab_size)
                x_axis.SetRangeUser(x_range[0], x_range[1]) if x_range is not None else do_nothing()
                x_axis.SetNdivisions(ndivx) if ndivx is not None else do_nothing()
                do(x_axis.SetTickSize, tick_size)
            y_axis = h.GetYaxis()
            if y_axis:
                y_axis.SetTitle(y_tit) if y_tit else y_axis.GetTitle()
                y_axis.SetTitleOffset(y_off) if y_off is not None else do_nothing()
                y_axis.SetTitleSize(tit_size)
                y_axis.SetLabelSize(lab_size)
                do(y_axis.SetLabelOffset, l_off_y)
                y_axis.SetRangeUser(y_range[0], y_range[1]) if y_range is not None else do_nothing()
                do(y_axis.SetNdivisions, ndivy)
            z_axis = h.GetZaxis()
            if z_axis:
                z_axis.SetTitle(z_tit) if z_tit else h.GetZaxis().GetTitle()
                z_axis.SetTitleOffset(z_off) if z_off is not None else do_nothing()
                z_axis.SetTitleSize(tit_size)
                z_axis.SetLabelSize(lab_size)
                z_axis.SetRangeUser(z_range[0], z_range[1]) if z_range is not None else do_nothing()
        except AttributeError or ReferenceError:
            pass
        set_time_axis(h, off=t_ax_off) if t_ax_off is not None else do_nothing()

    @staticmethod
    def make_tgrapherrors(name, title, color=1, marker=20, marker_size=1, width=1, asym_err=False, style=1, x=None, y=None, ex=None, ey=None):
        x = list(x) if type(x) == ndarray else x
        if x is None:
            gr = TGraphErrors() if not asym_err else TGraphAsymmErrors()
        else:
            gr = TGraphErrors(*make_graph_args(x, y, ex, ey)) if not asym_err else TGraphAsymmErrors(*make_graph_args(x, y, ex, ey))
        gr.SetTitle(title)
        gr.SetName(name)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
        gr.SetMarkerSize(marker_size)
        gr.SetLineWidth(width)
        gr.SetLineStyle(style)
        return gr

    @staticmethod
    def make_legend(x1=.65, y2=.88, nentries=2, scale=1, name='l', y1=None, clean=False, margin=.25, x2=None):
        x2 = .95 if x2 is None else x2
        y1 = y2 - nentries * .05 * scale if y1 is None else y1
        l = TLegend(x1, y1, x2, y2)
        l.SetName(name)
        l.SetTextFont(42)
        l.SetTextSize(0.03 * scale)
        l.SetMargin(margin)
        if clean:
            l.SetLineWidth(2)
            l.SetBorderSize(0)
            l.SetFillColor(0)
            l.SetFillStyle(0)
            l.SetTextAlign(12)
        return l

    def make_canvas(self, name='c', title='c', x=1., y=1., show=True, logx=None, logy=None, logz=None, gridx=None, gridy=None, transp=None):
        set_root_output(show)
        c = TCanvas(name, title, int(x * self.Res), int(y * self.Res))
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx, gridy])
        do(make_transparent, c, transp)
        self.Objects.append(c)
        return c

    def make_graph_from_profile(self, p):
        x_range = [i for i in xrange(p.GetNbinsX()) if p.GetBinContent(i)]
        x, y = [p.GetBinCenter(i) for i in x_range], [p.GetBinContent(i) for i in x_range]
        ex, ey = [p.GetBinWidth(i) / 2 for i in x_range], [p.GetBinError(i) for i in x_range]
        return self.make_tgrapherrors('g{n}'.format(n=p.GetName()[1:]), p.GetTitle(), x=x, y=y, ex=ex, ey=ey)

    @staticmethod
    def fix_chi2(g, prec=.01, show=True):
        it = 0
        error = 2
        chi2 = 0
        fit = None
        while abs(chi2 - 1) > prec and it < 20:
            for i in xrange(g.GetN()):
                g.SetPointError(i, g.GetErrorX(i), error)
            fit = g.Fit('pol0', 'qs{}'.format('' if show else 0))
            chi2 = fit.Chi2() / fit.Ndf()
            error += .5 ** it * sign(chi2 - 1)
            it += 1
        return fit if fit is not None else FitRes()


def create_colorlist():
    col_names = [kGreen, kOrange, kViolet, kYellow, kRed, kBlue, kMagenta, kAzure, kCyan, kTeal]
    colors = []
    for color in col_names:
        colors.append(color + (1 if color != 632 else -7))
    for color in col_names:
        colors.append(color + (3 if color != 800 else 9))
    return colors


def make_graph_args(x, y, ex=None, ey=None):
    if len(list(x)) != len(list(y)):
        log_warning('Arrays have different size!')
        return []
    lx = len(x)
    return [lx, array(x, 'd'), array(y, 'd'), array(ex, 'd') if ex is not None else zeros(lx), array(ey, 'd') if ey is not None else zeros(lx)]