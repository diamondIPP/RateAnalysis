#!/usr/bin/env python
# --------------------------------------------------------
#       Class for saving the root plots
# created on September 25th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from os.path import expanduser, basename
from helpers.info_legend import InfoLegend
from helpers.draw import *


class SaveDraw(Draw):

    Save = True

    ServerMountDir = None
    MountExists = None

    def __init__(self, analysis=None, results_dir=None, sub_dir=''):
        super(SaveDraw, self).__init__(analysis.MainConfig.FileName)

        # INFO
        SaveDraw.Save = Draw.Config.get_value('SAVE', 'save', default=False)
        self.Legends = InfoLegend(analysis)

        # Results
        self.ResultsDir = join(Draw.Dir, 'Results', choose(results_dir, default=analysis.TCString))
        self.SubDir = str(sub_dir)

        # Server
        SaveDraw.ServerMountDir = expanduser(Draw.Config.get_value('SAVE', 'server mount directory', default=None))
        SaveDraw.server_is_mounted()
        self.ServerDir = SaveDraw.load_server_save_dir(analysis)

    # ----------------------------------------
    # region SET
    def set_sub_dir(self, name):
        self.SubDir = name

    def set_results_dir(self, name):
        self.ResultsDir = join(Draw.Dir, 'Results', name)
    # endregion SET
    # ----------------------------------------

    @staticmethod
    def server_is_mounted():
        if SaveDraw.MountExists is not None:
            return SaveDraw.MountExists
        SaveDraw.MountExists = dir_exists(join(SaveDraw.ServerMountDir, 'data'))
        if not SaveDraw.MountExists:
            warning('Diamond server is not mounted in {}'.format(SaveDraw.ServerMountDir))

    @staticmethod
    def load_server_save_dir(ana):
        if not SaveDraw.MountExists or SaveDraw.ServerMountDir is None or not hasattr(ana, 'DUT'):
            return
        run_string = f'RP-{ana.Ensemble.Name.lstrip("0").replace(".", "-")}' if hasattr(ana, 'RunPlan') else str(ana.Run.Number)
        return join(SaveDraw.ServerMountDir, 'content', 'diamonds', ana.DUT.Name, ana.TCString, run_string)

    # ----------------------------------------
    # region SAVE
    def save_full(self, h, filename, cname='c', **kwargs):
        self(h, **prep_kw(kwargs, show=False, save=False))
        self.save_plots(None, full_path=join(self.Dir, filename), show=False, cname=cname, **kwargs)

    def histo(self, histo, file_name=None, show=True, all_pads=False, prnt=True, save=True, info_leg=True, *args, **kwargs):
        c = super(SaveDraw, self).histo(histo, show, *args, **kwargs)
        histo.SetTitle('') if not Draw.Title else do_nothing()
        if info_leg:
            self.Legends.draw(c, all_pads, show and Draw.Legend)
        self.save_plots(file_name, prnt=prnt, show=show, save=save)
        return c

    def save_plots(self, savename, sub_dir=None, canvas=None, full_path=None, prnt=True, ftype=None, show=True, save=True, cname=None, **kwargs):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        kwargs = prep_kw(kwargs, save=save, prnt=prnt, show=show)
        if not kwargs['save'] or not SaveDraw.Save or (savename is None and full_path is None):
            return
        canvas = get_last_canvas() if canvas is None else canvas
        if cname is not None:
            canvas.SetName(cname)
        update_canvas(canvas)
        try:
            self.__save_canvas(canvas, sub_dir=sub_dir, file_name=savename, ftype=ftype, full_path=full_path, **kwargs)
            return Draw.add(canvas)
        except Exception as inst:
            warning('Error saving plots ...:\n  {}'.format(inst))

    def __save_canvas(self, canvas, file_name, res_dir=None, sub_dir=None, full_path=None, ftype=None, prnt=True, show=True, **kwargs):
        """should not be used in analysis methods..."""
        _ = kwargs
        file_path = join(choose(res_dir, self.ResultsDir), choose(sub_dir, self.SubDir), file_name) if full_path is None else full_path
        file_name = basename(file_path)
        ensure_dir(dirname(file_path))
        info(f'saving plot: {file_name}', prnt=prnt and self.Verbose)
        canvas.Update()
        Draw.set_show(show)  # needs to be in the same batch so that the pictures are created, takes forever...
        set_root_warnings(False)
        for f in choose(make_list(ftype), default=['pdf'], decider=ftype):
            canvas.SaveAs('{}.{}'.format(file_path, f.strip('.')))
        self.save_on_server(canvas, file_name, save=full_path is None)
        Draw.set_show(True)

    def save_on_server(self, canvas, file_name, save=True):
        if self.ServerDir is not None and save:
            canvas.SetName('c')
            fname = join(self.ServerDir, f'{basename(file_name)}.root')
            canvas.SaveAs(fname)
            if not Draw.Show:
                print(join('https://diamond.ethz.ch', 'psi2', fname[len(self.ServerMountDir) + 1:].replace('.root', '.html')))

    @staticmethod
    def save_last(canvas=None, ext='pdf'):
        filename = input(f'Enter the name of the {ext}-file: ')
        choose(canvas, get_last_canvas()).SaveAs(f'{filename.split(".")[0]}.{ext}')
    # endregion SAVE
    # ----------------------------------------


if __name__ == '__main__':

    from src.analysis import Analysis
    z = SaveDraw(Analysis(), join(Draw.Dir, 'config', 'main.ini'))
