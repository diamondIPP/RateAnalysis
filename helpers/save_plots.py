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

    def __call__(self, *args, **kwargs):
        return self.histo(*args, **kwargs)

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
        SaveDraw.MountExists = dir_exists(join(SaveDraw.ServerMountDir, 'Diamonds'))
        if not SaveDraw.MountExists:
            warning('Diamond server is not mounted in {}'.format(SaveDraw.ServerMountDir))

    @staticmethod
    def load_server_save_dir(ana):
        if not SaveDraw.MountExists or SaveDraw.ServerMountDir is None or not hasattr(ana, 'DUT'):
            return
        run_string = 'RunPlan{}'.format(ana.RunPlan.lstrip('0')) if hasattr(ana, 'RunPlan') else str(ana.Run.Number)
        return join(SaveDraw.ServerMountDir, 'Diamonds', ana.DUT.Name, 'BeamTests', make_tc_str(ana.TCString, long_=False), run_string)

    # ----------------------------------------
    # region SAVE
    def histo(self, histo, file_name=None, show=True, all_pads=False, prnt=True, save=True, info_leg=True, *args, **kwargs):
        c = super(SaveDraw, self).histo(histo, show, *args, **kwargs)
        if info_leg:
            self.Legends.draw(c, all_pads, show and Draw.Legend)
        self.save_plots(file_name, prnt=prnt, show=show, save=save)
        return c

    def distribution(self, values, binning=None, title='', file_name=None, show=True, prnt=True, save=True, w=1, h=1, *args, **kwargs):
        th = super(SaveDraw, self).distribution(values, binning, title, show=show, w=w, h=h, *args, **kwargs)
        self.save_plots(file_name, prnt=prnt, show=show, save=save)
        return th

    def save_plots(self, savename, sub_dir=None, canvas=None, prnt=True, ftype=None, show=True, save=True):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        if not save or not SaveDraw.Save or savename is None:
            return
        canvas = get_last_canvas() if canvas is None else canvas
        update_canvas(canvas)
        try:
            self.__save_canvas(canvas, sub_dir=sub_dir, file_name=savename, ftype=ftype, prnt=prnt, show=show)
            return Draw.add(canvas)
        except Exception as inst:
            warning('Error saving plots ...:\n  {}'.format(inst))

    def __save_canvas(self, canvas, file_name, res_dir=None, sub_dir=None, ftype=None, prnt=True, show=True):
        """should not be used in analysis methods..."""
        file_path = join(choose(res_dir, self.ResultsDir), choose(sub_dir, self.SubDir), file_name)
        ensure_dir(dirname(file_path))
        info('saving plot: {}'.format(file_name), prnt=prnt and SaveDraw.Verbose)
        canvas.Update()
        set_root_output(show)  # needs to be in the same batch so that the pictures are created, takes forever...
        set_root_warnings(False)
        for f in choose(make_list(ftype), default=['pdf'], decider=ftype):
            canvas.SaveAs('{}.{}'.format(file_path, f.strip('.')))
        self.save_on_server(canvas, file_name, ftype)
        set_root_output(True)

    def save_on_server(self, canvas, file_name, ftype=None):
        if self.ServerDir is not None:
            for ft in ['pdf', 'png'] if ftype is None else make_list(ftype):
                canvas.SaveAs('{}.{}'.format(join(self.ServerDir, basename(file_name)), ft))

    @staticmethod
    def server_pickle(old_path, value):
        if SaveDraw.server_is_mounted():
            picklepath = join(SaveDraw.ServerMountDir, 'Pickles', basename(dirname(old_path)), basename(old_path))
            do_pickle(picklepath, do_nothing, value)
    # endregion SAVE
    # ----------------------------------------


if __name__ == '__main__':

    from src.analysis import Analysis
    z = SaveDraw(Analysis(), join(Draw.Dir, 'config', 'main.ini'))