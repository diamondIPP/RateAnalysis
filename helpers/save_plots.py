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
    TCString = None

    ServerMountDir = None
    MountExists = None

    def __init__(self, analysis=None, results_dir=None, sub_dir=''):
        super(SaveDraw, self).__init__(analysis.MainConfig.FileName)

        # INFO
        SaveDraw.Save = Draw.Config.get_value('SAVE', 'save', default=False)
        SaveDraw.TCString = analysis.TCString
        self.Legend = InfoLegend(analysis)
        Draw.Legend = self.Legend.is_active() and Draw.Legend

        # Results
        self.ResultsDir = join(Draw.Dir, 'Results', choose(results_dir, default=SaveDraw.TCString))
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
        run_string = 'RunPlan{}'.format(ana.RunPlan.lstrip('0')) if hasattr(ana, 'RunPlan') else str(ana.RunNumber)
        return join(SaveDraw.ServerMountDir, 'Diamonds', ana.DUT.Name, 'BeamTests', make_tc_str(SaveDraw.TCString, long_=False), run_string)

    # ----------------------------------------
    # region SAVE
    def histo(self, histo, file_name=None, show=True, all_duts=False, all_pads=False, prnt=True, *args, **kwargs):
        c = super(SaveDraw, self).histo(histo, show, *args, **kwargs)
        self.Legend.draw(c, all_pads, all_duts)
        if file_name is not None and SaveDraw.Save:
            self.save_plots(file_name, prnt=prnt, show=show)
        return c

    def save_plots(self, savename, sub_dir=None, canvas=None, prnt=True, show=True, save=True):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        if not save:
            return
        canvas = get_last_canvas() if canvas is None else canvas
        update_canvas(canvas)
        try:
            self.__save_canvas(canvas, sub_dir=sub_dir, file_name=savename, prnt=prnt, show=show)
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
