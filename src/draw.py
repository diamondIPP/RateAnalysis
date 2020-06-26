#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from __future__ import print_function
from utils import *
from ROOT import TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TCanvas, gStyle, TLegend, TArrow, TPad, TCutG, TLine, TPaveText, TPaveStats, TH1F, TSpectrum, TEllipse, TColor, TProfile
from numpy import ndarray, zeros, sign, linspace, ones
from os.path import expanduser, join, basename


class Draw:

    Count = 0
    Res = None

    def __init__(self, tc_string='', verbose=True, config=None):

        # Basics
        self.load_resolution(default=800)
        self.Verbose = verbose
        self.Dir = get_base_dir()
        self.TCString = tc_string
        self.MainConfig = self.init_main_config(config)
        self.ResultsDir = self.get_results_dir()
        self.SubDir = ''
        self.ServerDir = self.load_server_dir()

        # Colors
        self.Color = 0
        self.Colors = self.get_colors(10)
        self.FillColor = 821

        # Settings
        gStyle.SetLegendFont(42)
        self.Title = self.get_config('SAVE', 'activate title')
        set_titles(self.Title)
        self.Legend = self.get_config('SAVE', 'info legend')
        self.Save = self.get_config('SAVE', 'save')

        self.Objects = []

    def init_main_config(self, config):
        return load_parser(join(self.Dir, 'config', 'main.ini')) if config is None else config

    @staticmethod
    def load_resolution(default=800):
        if Draw.Res is None:
            Draw.Res = load_resolution(default)

    # ----------------------------------------
    # region BASIC
    def set_save_directory(self, name):
        self.SubDir = name

    def set_results_dir(self, name):
        self.ResultsDir = join(self.Dir, name)

    def get_results_dir(self):
        return join(self.Dir, 'Results', self.TCString)

    def get_color(self, n, i=None):
        color = get_color_gradient(n)[choose(i, self.Color)]
        self.Color = self.Color + 1 if self.Color < n - 1 else 0
        return color

    @staticmethod
    def get_colors(n):
        return get_color_gradient(n)

    def get_config(self, section, option):
        return True if self.MainConfig is None else self.MainConfig.getboolean(section, option)

    def load_server_dir(self):
        return expanduser(self.MainConfig.get('SAVE', 'server mount directory')) if self.MainConfig is not None else None

    def add(self, *args):
        for obj in args:
            if obj not in self.Objects:
                self.Objects.append(obj)
        self.clean()

    def clean(self):
        n_none = sum(str(obj) == 'None' for obj in self.Objects)
        for _ in range(n_none):
            self.Objects.remove(None)

    def set_pad_margins(self, c=None, l=None, r=None, b=None, t=None):
        do(c.SetLeftMargin, l)
        do(c.SetRightMargin, r if r is not None else None if round(c.GetRightMargin(), 1) != .1 else .03 )
        do(c.SetBottomMargin, None if round(c.GetBottomMargin(), 1) != .1 else (.17 if b is None else b) - (.07 if not self.Legend else 0))
        do(c.SetTopMargin, None if round(c.GetTopMargin(), 1) != .1 else (.1 if t is None else t) - (0 if self.Title else .07))
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region DRAWING
    def draw_axis(self, x1, x2, y1, y2, title, limits=None, name='ax', col=1, width=1, off=.15, tit_size=.035, lab_size=0.035, tick_size=0.03, line=False, opt='+SU', l_off=.01, log=False):
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
        a.SetLabelFont(42)
        a.SetTitleFont(42)
        a.SetTickSize(tick_size if not line else 0)
        a.SetTickLength(tick_size if not line else 0)
        a.SetNdivisions(0) if line else do_nothing()
        a.SetLabelOffset(l_off)
        a.Draw()
        self.add(a)
        return a

    def draw_y_axis(self, x, ymin, ymax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return self.draw_axis(x, x, ymin, ymax, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    def draw_x_axis(self, y, xmin, xmax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return self.draw_axis(xmin, xmax, y, y, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    def draw_line(self, x1, x2, y1, y2, color=1, width=1, style=1, name='li'):
        line = TCutG(name, 2, array([x1, x2], 'd'), array([y1, y2], 'd'))
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw('same')
        self.add(line)
        return line

    def draw_tline(self, x1, x2, y1, y2, color=1, width=1, style=1):
        line = TLine(x1, y1, x2, y2)
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw()
        self.add(line)
        return line

    def draw_box(self, x1, y1, x2, y2, line_color=1, width=1, style=1, fillstyle=None, name='box', show=True):
        return self.draw_n_box([x1, x1, x2, x2], [y1, y2, y2, y1], line_color, width, style, fillstyle, name, show)

    def draw_n_box(self, x, y, line_color=1, width=1, style=1, fillstyle=None, name='box', show=True, fill_color=None):
        line = TCutG(name, len(x) + 1, array(x + [x[0]], 'd'), array(y + [y[0]], 'd'))
        line.SetLineColor(line_color)
        do(line.SetFillColor, fill_color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.SetFillStyle(fillstyle) if fillstyle is not None else do_nothing()
        if show:
            line.Draw('l')
            line.Draw('f') if fillstyle < 4000 else do_nothing()
        self.add(line)
        return line

    def draw_vertical_line(self, x, ymin, ymax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(x, x, ymin, ymax, color, w, style, name) if not tline else self.draw_tline(x, x, ymin, ymax, color, w, style)

    def draw_horizontal_line(self, y, xmin, xmax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(xmin, xmax, y, y, color, w, style, name) if not tline else self.draw_tline(xmin, xmax, y, y, color, w, style)

    def draw_tlatex(self, x, y, text, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None, show=True):
        tlatex = TLatex(x, y, text)
        format_text(tlatex, name, align, color, size, angle, ndc, font)
        tlatex.Draw() if show else do_nothing()
        self.add(tlatex)
        return tlatex

    def draw_arrow(self, x1, x2, y1, y2, col=1, width=1, opt='<|', size=.005):
        ar = TArrow(x1, y1, x2, y2, size, opt)
        ar.SetLineWidth(width)
        ar.SetLineColor(col)
        ar.SetFillColor(col)
        ar.Draw()
        self.add(ar)

    def draw_tpad(self, name, tit='', pos=None, fill_col=0, gridx=None, gridy=None, margins=None, transparent=False, logy=None, logx=None, logz=None, lm=None, rm=None, bm=None, tm=None):
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(name, tit, *pos)
        p.SetFillColor(fill_col)
        self.set_pad_margins(p, *(margins if all(m is None for m in [lm, rm, bm, tm]) else [lm, rm, bm, tm]))
        do([p.SetLogx, p.SetLogy, p.SetLogz], [logx, logy, logz])
        do([p.SetGridx, p.SetGridy], [gridx, gridy])
        make_transparent(p) if transparent else do_nothing()
        p.Draw()
        p.cd()
        self.add(p)
        return p

    def draw_tpavetext(self, text, x1, x2, y1, y2, font=42, align=0, size=0, angle=0, margin=.05, color=1):
        p = TPaveText(x1, y1, x2, y2, 'ndc')
        p.SetFillColor(0)
        p.SetFillStyle(0)
        p.SetBorderSize(0)
        p.SetMargin(margin)
        t = p.AddText(text)
        format_text(t, 'pave', align, color, size, angle, ndc=True, font=font)
        p.Draw()
        self.add(p)
        return p

    def draw_stats(self, fit, y2=None, width=.3, prec='5.1f', names=None):
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
        ls.Add(self.draw_tlatex(0, 0, '#chi^{{2}} / ndf  = {chi2:{p}} / {ndf}'.format(ndf=fit.Ndf(), chi2=fit.Chi2(), p=prec), size=0, align=0, font=42))
        for i in range(fit.NPars):
            ls.Add(self.draw_tlatex(0, 0, '{n}  = {v:{p}} #pm {e:{p}}'.format(n=names[i], v=fit.Parameter(i), e=fit.ParError(i), p=prec), size=0, align=0, font=42))
        p.Draw()
        self.add(p)
        return p

    def draw_frame(self, pad, xmin, xmax, ymin, ymax, tit, div=None, y_cent=None):
        pad.cd()
        fr = pad.DrawFrame(xmin, ymin, xmax, ymax)
        pad.Modified()
        fr.GetYaxis().SetTitle(tit)
        do(fr.GetYaxis().CenterTitle, y_cent)
        fr.GetYaxis().SetNdivisions(div) if div is not None else do_nothing()
        format_frame(fr)
        self.add(fr)

    def draw_grid(self, x_vals, y_vals, width=1, color=1):
        for x in x_vals:
            self.draw_line(x, x, min(y_vals), max(y_vals), name='x{}'.format(x), width=width, color=color)
        for y in y_vals:
            self.draw_line(min(x_vals), max(x_vals), y, y, name='y{}'.format(y), width=width, color=color)

    def draw_ellipse(self, a, b=0, x_off=0, y_off=0, color=2, w=2):
        e = TEllipse(x_off, y_off, a, b)
        do(e.SetLineColor, color)
        do(e.SetLineWidth, w)
        e.SetFillStyle(4000)
        e.Draw()
        self.add(e)

    def draw_circle(self, r, x_off=0, y_off=0, color=None, w=None):
        self.draw_ellipse(r, 0, x_off, y_off, color, w)

    def draw_preliminary(self, canvas=None, height=.06):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        return self.draw_tpavetext('#font[62]{RD42} Preliminary', c.GetLeftMargin(), .5, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=72, align=12, margin=0.04)

    def draw_irradiation(self, irr, canvas=None, height=.06, left=True):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        x1, x2 = (c.GetLeftMargin(), .5) if left else (.5, 1 - c.GetRightMargin())
        return self.draw_tpavetext('Irradiation: {}'.format(irr), x1, x2, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=42, align=12, margin=0.04)

    def draw_histo(self, histo, show=True, lm=None, rm=None, bm=None, tm=None, draw_opt=None, x=None, y=None, all_pads=True,
                   leg=None, logy=False, logx=False, logz=False, canvas=None, grid=False, gridy=False, gridx=False, both_dias=False, prnt=True, phi=None, theta=None, ind=None):
        return self.save_histo(histo, '', show, '', lm, rm, bm, tm, draw_opt, x, y, all_pads, leg, logy, logx, logz, canvas, grid, gridx, gridy, False, both_dias, ind,
                               prnt, phi, theta)

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
        p0 = self.draw_tpad('p0', 'p0', pos=[0, 0, 1, pm], margins=[.14, .03, bm / pm, 0], transparent=True, logx=True, gridy=True)
        scale_multigraph(mgn)
        rel_y_range = [.7, 1.3] if rel_y_range is None else rel_y_range
        format_histo(mgn, title='', y_range=rel_y_range, y_tit='Rel. ph [au]' if not scale > 1 else ' ', y_off=66, tit_size=.1 * scale, x_off=99, lab_size=.1 * scale)
        mgn.GetYaxis().SetNdivisions(3)
        hide_axis(mgn.GetXaxis())
        mgn.Draw('alp')
        x_range = [mgn.GetXaxis().GetXmin(), mgn.GetXaxis().GetXmax()] if x_range is None else x_range
        self.draw_x_axis(1.3, x_range[0], x_range[1], mgn.GetXaxis().GetTitle() + ' ', opt='SG+-=', tit_size=.1, lab_size=.1 * scale, off=99, tick_size=.1, l_off=0)
        c.cd()

        # top pad with zero suppression
        self.draw_tpad('p1', 'p1', pos=[0, pm, 1, 1], margins=[.14, .03, 0, .1], transparent=True, logx=True)
        mg.Draw('alp')
        hide_axis(mg.GetXaxis())
        if pulser_leg:
            pulser_leg()
        if y_range:
            mg.SetMinimum(y_range[0])
            mg.SetMaximum(y_range[1])
        format_histo(mg, tit_size=.04 * scale, y_off=1.75 / scale, lab_size=.04 * scale)
        self.draw_x_axis(mg_y, x_range[0], x_range[1], mgn.GetXaxis().GetTitle() + ' ', opt='SG=', tit_size=.035 * scale, lab_size=0, off=1, l_off=99)
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

        self.add(c, *draw_objects)
        set_root_output(True)

    @staticmethod
    def make_bins(values, thresh=.02):
        bins = linspace(*(find_range(values, thresh=thresh) + [int(sqrt(values.size))]))
        return [bins.size - 1, bins]

    @staticmethod
    def get_count():
        Draw.Count += 1
        return Draw.Count

    def draw_disto(self, values, title='', bins=None, thresh=.02, lm=None, rm=None, show=True, **kwargs):
        values = array(values, dtype='d')
        kwargs['fill_color'] = self.FillColor if 'fill_color' not in kwargs else kwargs['fill_color']
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        kwargs['y_tit'] = 'Number of Entries' if 'y_tit' not in kwargs else kwargs['y_tit']
        h = TH1F('h{}'.format(self.get_count()), title, *choose(bins, self.make_bins, values=values, thresh=thresh))
        fill_hist(h, values)
        format_histo(h, **kwargs)
        self.draw_histo(h, show, lm, rm)
        return h

    def draw_profile(self, x, y, bins=None, title='', thresh=.02, lm=None, rm=None, show=True, **kwargs):
        x, y = array(x, dtype='d'), array(y, dtype='d')
        kwargs['fill_color'] = self.FillColor if 'fill_color' not in kwargs else kwargs['fill_color']
        kwargs['y_off'] = 1.4 if 'y_off' not in kwargs else kwargs['y_off']
        p = TProfile('p{}'.format(self.get_count()), title, *choose(bins, self.make_bins, values=x, thresh=thresh))
        fill_hist(p, x, y)
        format_histo(p, **kwargs)
        self.draw_histo(p, show, lm, rm)
        return p

    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region SAVE
    def make_bias_string(self, bias=None):
        if bias is None:
            return self.make_bias_string(self.bias) if hasattr(self, 'bias') else ''
        pol = 'm' if bias < 0 else 'p'
        return '_{pol}{bias:04d}'.format(pol=pol, bias=int(abs(bias))) if bias else ''

    def make_info_string(self):
        string = '_{dia}'.format(dia=self.DUT.Name) if hasattr(self, 'DUT') else ''
        string += self.make_bias_string()
        string += '_{tc}'.format(tc=self.TCString)
        string = string.replace('-', '')
        return string

    def server_is_mounted(self):
        return dir_exists(join(self.ServerDir, 'Diamonds'))

    def save_on_server(self, canvas, file_name, ftype=None):
        if self.ServerDir is None:
            return
        if not self.server_is_mounted():
            warning('Diamond server is not mounted in {}'.format(self.ServerDir))
            return
        if hasattr(self, 'DUT'):
            if hasattr(self, 'RunPlan') and self.RunPlan is not None:
                run_string = 'RunPlan{}'.format(self.RunPlan.lstrip('0'))
            elif hasattr(self, 'RunNumber'):
                run_string = str(self.RunNumber)
            else:
                return
            path = join(self.ServerDir, 'Diamonds', self.DUT.Name, 'BeamTests', make_tc_str(self.TCString, long_=False), run_string, file_name)
            for ft in ['pdf', 'png'] if ftype is None else [ftype]:
                canvas.SaveAs('{}.{}'.format(path, ft))

    def server_pickle(self, old_path, value):
        if self.server_is_mounted():
            picklepath = join(self.ServerDir, 'Pickles', basename(dirname(old_path)), basename(old_path))
            do_pickle(picklepath, do_nothing, value)

    def save_plots(self, savename, sub_dir=None, canvas=None, all_pads=True, both_dias=False, ind=None, prnt=True, save=True, show=True):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        canvas = get_last_canvas() if canvas is None else canvas
        if canvas is None:
            return
        if ind is None:
            try:
                self.InfoLegend.draw(canvas, all_pads, both_dias) if hasattr(self, 'InfoLegend') else log_warning('Did not find InfoLegend class...') \
                    if not any(hasattr(self, att) for att in ['RunSelections', 'CurrentGraph']) else do_nothing()
            except Exception as err:
                warning(err)
        else:
            self.Analyses.values()[ind].InfoLegend.draw(canvas, all_pads, both_dias) if hasattr(self, 'Analyses') else log_critical('sth went wrong...')
        canvas.Modified()
        canvas.Update()
        if save and self.Save:
            try:
                if both_dias and sub_dir is None:
                    sub_dir = self.TelSaveDir if hasattr(self, 'TelSaveDir') else sub_dir
                self.save_canvas(canvas, sub_dir=sub_dir, name=savename, print_names=prnt, show=show)
                self.add(canvas)
            except Exception as inst:
                log_warning('Error in save_canvas:\n{0}'.format(inst))

    def save_tel_plots(self, savename, sub_dir=None, canvas=None, all_pads=True, ind=None, prnt=True, save=True, show=True):
        self.save_plots(savename, sub_dir, canvas, all_pads, True, ind, prnt, save, show)

    def save_canvas(self, canvas, sub_dir=None, name=None, print_names=True, show=True, ftype=None, res_dir=None):
        """should not be used in analysis methods..."""
        sub_dir = self.SubDir if sub_dir is None else sub_dir
        canvas.Update()
        file_name = canvas.GetName() if name is None else name
        ftypes = ['root', 'png', 'pdf', 'eps'] if ftype is None else [ftype]
        file_path = join(self.ResultsDir if res_dir is None else res_dir, sub_dir, '{typ}' if len(ftypes) > 1 else '', file_name)
        out = 'saving plot: {nam}'.format(nam=name)
        run_number = self.RunNumber if hasattr(self, 'RunNumber') else None
        run_number = 'rp{nr}'.format(nr=self.run_plan) if hasattr(self, 'run_plan') else run_number
        set_root_output(show)  # needs to be in the same batch so that the pictures are created, takes forever...
        set_root_warnings(False)
        info_str = self.make_info_string()
        for f in ftypes:
            ext = '.{typ}'.format(typ=f)
            if not f == 'png' and run_number is not None:
                ext = '{str}_{run}.{typ}'.format(str=info_str, run=run_number, typ=f)
            ensure_dir(dirname(file_path.format(typ=f)))
            out_file = '{fname}{ext}'.format(fname=file_path, ext=ext)
            out_file = out_file.format(typ=f)
            canvas.SaveAs(out_file)
        self.save_on_server(canvas, file_name, ftype)
        if print_names:
            log_info(out, prnt=self.Verbose)
        set_root_output(True)

    def save_histo(self, histo, save_name='test', show=True, sub_dir=None, lm=None, rm=None, bm=None, tm=None, draw_opt=None, x=None, y=None, all_pads=True,
                   leg=None, logy=False, logx=False, logz=False, canvas=None, grid=False, gridx=False, gridy=False, save=True, both_dias=False, ind=None, prnt=True, phi=None, theta=None, sumw2=False):
        fac = 1 if self.Title else 1.16
        x = int(Draw.Res * fac) if x is None else int(x * Draw.Res)
        y = Draw.Res if y is None else int(y * Draw.Res)
        h = histo
        h.Sumw2(sumw2) if hasattr(h, 'Sumw2') and sumw2 is not None else do_nothing()
        set_root_output(show)
        c = TCanvas('c_{0}'.format(h.GetName()), h.GetTitle().split(';')[0], x, y) if canvas is None else canvas
        self.set_pad_margins(c, lm, rm, bm, tm)
        c.SetLogx() if logx else do_nothing()
        c.SetLogy() if logy else do_nothing()
        c.SetLogz() if logz else do_nothing()
        c.SetGridx() if gridx or grid else do_nothing()
        c.SetGridy() if gridy or grid else do_nothing()
        c.SetPhi(phi) if phi is not None else do_nothing()
        c.SetTheta(theta) if theta is not None else do_nothing()
        h.Draw(draw_opt if draw_opt is not None else 'ap' if 'Graph' in h.ClassName() else '')
        if leg is not None:
            leg = [leg] if type(leg) is not list else leg
            for i in leg:
                i.Draw()
        self.save_plots(save_name, sub_dir=sub_dir, both_dias=both_dias, all_pads=all_pads, ind=ind, prnt=prnt, save=save, show=show)
        set_root_output(True)
        self.add(c, h, leg)
        return c

    def save_tel_histo(self, histo, save_name='test', show=True, sub_dir=None, lm=.1, rm=.03, bm=.15, tm=None, draw_opt=None, x_fac=None, y_fac=None, all_pads=True,
                       leg=None, logy=False, logx=False, logz=False, canvas=None, grid=False, gridx=False, gridy=False, save=True, ind=None, prnt=True, phi=None, theta=None):
        return self.save_histo(histo, save_name, show, sub_dir, lm, rm, bm, tm, draw_opt, x_fac, y_fac, all_pads, leg, logy, logx, logz, canvas, grid, gridx, gridy, save, True, ind, prnt, phi, theta)
    # endregion SAVE
    # ----------------------------------------

    # ----------------------------------------
    # region CREATE
    def make_tgrapherrors(self, name, title, color=1, marker=20, marker_size=1, width=1, asym_err=False, style=1, x=None, y=None, ex=None, ey=None):
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
        self.add(gr)
        return gr

    def make_legend(self, x1=.65, y2=.88, nentries=2, scale=1, name='l', y1=None, clean=False, margin=.25, x2=None, w=None, cols=None):
        x2 = .95 if x2 is None else x2
        x1 = x2 - w if w is not None else x1
        h = nentries * .05 * scale
        y = array([y2 - h if y1 is None else y1, y1 + h if y1 is not None else y2])
        y += .07 if not self.Title and y[1] > .7 else 0
        y -= .07 if not self.Legend and y[1] < .7 else 0
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
        self.add(leg)
        return leg

    def make_canvas(self, name='c', title='c', x=1., y=1., logx=None, logy=None, logz=None, gridx=None, gridy=None, transp=None, divide=None, show=True):
        set_root_output(show)
        c = TCanvas(name, title, int(x * Draw.Res), int(y * Draw.Res))
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx, gridy])
        do(make_transparent, c, transp)
        if divide is not None:
            c.Divide(*(divide if type(divide) in [list, tuple] else [divide]))
        self.add(c)
        return c

    def make_graph_from_profile(self, p):
        x_range = [i for i in range(p.GetNbinsX()) if p.GetBinContent(i)]
        x = [make_ufloat([p.GetBinCenter(i), p.GetBinWidth(i) / 2]) for i in x_range]
        y = [make_ufloat([p.GetBinContent(i), p.GetBinError(i)]) for i in x_range]
        return self.make_tgrapherrors('g{n}'.format(n=p.GetName()[1:]), p.GetTitle(), x=x, y=y)
    # endregion CREATE
    # ----------------------------------------

    def format_statbox(self, x=.95, y=None, w=.2, h=.15, only_fit=False, fit=False, entries=False, form=None, m=False, rms=False, all_stat=False):
        gStyle.SetOptFit(int(only_fit or fit))
        opt_stat = '100000{}{}{}0'.format(*[int(v) for v in [rms, m, entries]] if not all_stat else [1, 1, 1])
        if only_fit:
            opt_stat = '0011'
        y = (.88 if self.Title else .95) if y is None else y
        gStyle.SetOptStat(int(opt_stat))
        gStyle.SetFitFormat(form) if form is not None else do_nothing()
        gStyle.SetStatX(x)
        gStyle.SetStatY(y)
        gStyle.SetStatW(w)
        gStyle.SetStatH(h)
    # END OF CLASS


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
        x_args = [x_tit, x_off, tit_size, center_x, lab_size, l_off_x, x_range, ndivx, max(tick_size, x_ticks), ]
        y_args = [y_tit, y_off, tit_size, center_y, lab_size, l_off_y, y_range, ndivy, max(tick_size, y_ticks), yax_col]
        z_args = [z_tit, z_off, tit_size, False, lab_size, None, z_range, None, max(tick_size, z_ticks)]
        for i, name in enumerate(['X', 'Y', 'Z']):
            format_axis(getattr(h, 'Get{}axis'.format(name))(), is_graph(h), *[x_args, y_args, z_args][i])
    except AttributeError or ReferenceError:
        pass
    set_time_axis(h, off=t_ax_off) if t_ax_off is not None else do_nothing()
    do(h.Sumw2, sumw2) if hasattr(h, 'Sumw2') else do_nothing()


def format_axis(axis, graph, title, tit_offset, tit_size, centre_title, lab_size, label_offset, limits, ndiv, tick_size, color=None):
    do(axis.SetTitle, title)
    do(axis.SetTitleOffset, tit_offset)
    do(axis.SetTitleSize, tit_size)
    axis.CenterTitle(centre_title)
    do(axis.SetLabelSize, lab_size)
    do(axis.SetLabelOffset, label_offset)
    if limits is not None:
        axis.SetLimits(*limits) if graph and 'xaxis' in axis.GetName() else axis.SetRangeUser(*limits)
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
# endregion
# ----------------------------------------


def set_2d_ranges(h, dx, dy):
    # find centers in x and y
    xmid, ymid = [(p.GetBinCenter(p.FindFirstBinAbove(0)) + p.GetBinCenter(p.FindLastBinAbove(0))) / 2 for p in [h.ProjectionX(), h.ProjectionY()]]
    format_histo(h, x_range=[xmid - dx, xmid + dx], y_range=[ymid - dy, ymid + dx])


def adapt_z_range(h, n_sigma=2):
    values = get_2d_hist_vec(h)
    m, s = mean_sigma(values[5:-5])
    z_range = [min(values).n, .8 * max(values).n] if s > m else [m - n_sigma * s, m + n_sigma * s]
    format_histo(h, z_range=z_range)


def find_range(values, lfac=.2, rfac=.2, thresh=.02):
    v = array(sorted(values))
    xmin, xmax = v[int(thresh * v.size)], v[int(v.size - thresh * v.size)]
    return increased_range([xmin, xmax], lfac, rfac)


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


def make_graph_args(x, y, ex=None, ey=None):
    if len(list(x)) != len(list(y)):
        log_warning('Arrays have different size!')
        return []
    lx = len(x)
    x = x if type(x[0]) in [Variable, AffineScalarFunc] else [make_ufloat(tup) for tup in zip(x, zeros(lx) if ex is None else ex)]
    y = y if type(y[0]) in [Variable, AffineScalarFunc] else [make_ufloat(tup) for tup in zip(y, zeros(lx) if ey is None else ey)]
    return [lx, array([v.n for v in x], 'd'), array([v.n for v in y], 'd'), array([v.s for v in x], 'd'), array([v.s for v in y], 'd')]


def set_titles(status=True):
    gStyle.SetOptTitle(status)


def load_resolution(default=800):
    try:
        from screeninfo import get_monitors
        return int(round_up_to(get_monitors()[0].height / 2, 100))
    except Exception as err:
        warning(err)
        return default


def get_graph_vecs(g):
    return get_graph_x(g), get_graph_y(g)


def get_graph_x(g):
    return array([make_ufloat([g.GetX()[i], g.GetEX()[i]]) for i in range(g.GetN())]) if 'Error' in g.ClassName() else array([make_ufloat(g.GetX()[i]) for i in range(g.GetN())])


def get_graph_y(g):
    return array([make_ufloat([g.GetY()[i], g.GetEY()[i]]) for i in range(g.GetN())]) if 'Error' in g.ClassName() else array([make_ufloat(g.GetY()[i]) for i in range(g.GetN())])


def get_hist_vec(p, err=True):
    return array([make_ufloat([p.GetBinContent(ibin), p.GetBinError(ibin)]) if err else p.GetBinContent(ibin) for ibin in range(1, p.GetNbinsX() + 1)])


def get_hist_args(p, err=True):
    return array([make_ufloat([p.GetBinCenter(ibin), p.GetBinWidth(ibin) / 2]) if err else p.GetBinCenter(ibin) for ibin in range(1, p.GetNbinsX() + 1)])


def get_hist_vecs(p, err=True):
    return get_hist_args(p, err), get_hist_vec(p, err)


def get_bin_entries(h):
    return array([h.GetBinEntries(ibin) for ibin in range(1, h.GetNbinsX() + 1)])


def get_h_values(h):
    return get_graph_y(h) if 'Graph' in h.ClassName() else get_hist_vec(h)


def get_h_args(h):
    return get_graph_x(h) if 'Graph' in h.ClassName() else get_hist_args(h)


def get_2d_hist_vec(h):
    xbins, ybins = range(1, h.GetNbinsX() + 1), range(1, h.GetNbinsY() + 1)
    return array([make_ufloat([h.GetBinContent(xbin, ybin), h.GetBinError(xbin, ybin)]) for xbin in xbins for ybin in ybins if h.GetBinContent(xbin, ybin)])


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
        m, s = mean_sigma(y)
        scale = val / (y[where(x == min(x))[0]] if to_low_flux else m)
    for i in range(x.size):
        gr.SetPoint(i, gr.GetX()[i], gr.GetY()[i] * scale)
        gr.SetPointError(i, gr.GetErrorX(i), gr.GetErrorY(i) * scale) if 'Error' in gr.ClassName() else do_nothing()
    return scale


def get_pull(h, name, bins, fit=True):
    set_root_output(False)
    h_out = TH1F('hp{}'.format(name[:3]), name, *bins)
    values = array([h.GetBinContent(ibin + 1) for ibin in range(h.GetNbinsX())], 'd')
    h_out.FillN(values.size, values, full(values.size, 1, 'd'))
    h_out.Fit('gaus', 'q') if fit else do_nothing()
    format_histo(h_out, x_range=increased_range([values.min(), values.max()], .1, .3))
    return h_out


def markers(i):
    return ((range(20, 24) + [29, 33, 34]) * 2)[i]


def fit_bucket(histo, show=True):
    set_root_output(False)
    h = histo
    format_histo(h, rebin=int(h.GetBinCenter(h.FindLastBinAbove(h.GetMaximum() * .02))) / 40)
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


def fill_hist(h, x, y=None, zz=None):
    x, y, zz = array(x).astype('d'), array(y).astype('d'), array(zz).astype('d')
    if h.ClassName == 'TProfile2D':
        for i in range(x.size):
            h.Fill(x[i], y[i], zz[i])
    elif 'TH1' in h.ClassName():
        h.FillN(x.size, x, ones(x.size))
    elif any(name in h.ClassName() for name in ['TH2', 'TProfile']):
        h.FillN(x.size, x, y, ones(x.size))
    else:
        h.FillN(x.size, x, y, zz, ones(x.size))


def set_palette(pal):
    gStyle.SetPalette(pal)


def is_graph(h):
    return 'Graph' in h.ClassName()


def update_canvas(c=None):
    c = choose(c, get_last_canvas())
    c.Modified()
    c.Update()


def get_color_gradient(n):
    stops = array([0., .5, 1], 'd')
    green = array([0. / 255., 200. / 255., 80. / 255.], 'd')
    blue = array([0. / 255., 0. / 255., 0. / 255.], 'd')
    red = array([180. / 255., 200. / 255., 0. / 255.], 'd')
    color_gradient = TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 255)
    color_table = [color_gradient + ij for ij in xrange(255)]
    return color_table[0::(len(color_table) + 1) / n]


if __name__ == '__main__':
    z = Draw()
