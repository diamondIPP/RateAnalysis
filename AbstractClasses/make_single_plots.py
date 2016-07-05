#!/usr/bin/python

from PadAnalysis import PadAnalysis
from Elementary import Elementary

tc = '201510'
a = Elementary(tc)
a.print_banner('STARTING SINGLE PLOT GENERATION')
a.print_testcampaign()


def set_save_dirs(ana):
    ana.save_dir = '{info}_{rp}'.format(info=ana.make_info_string().strip('_'), rp=ana.run_number)
    ana.ana_save_dir = '{info}_{rp}'.format(info=ana.make_info_string().strip('_'), rp=ana.run_number)
    ana.set_save_directory('PlotsFelix')
    

# for run in [392, 398]:
for run in [171, 178, 392, 398]:
    z = PadAnalysis(run, 0)
    set_save_dirs(z)
    if run in [178, 398]:
        z.draw_bucket_waveforms(show=False, start=110000)
    z.draw_trigger_cell_vs_peakpos(show=False, corr=False)
    z.draw_trigger_cell_vs_peakpos(show=False, corr=True)
    z.draw_trigger_cell_vs_peakpos(show=False, corr=True, t_corr=True)
    z.draw_bucket_pedestal(show=False)
    z.draw_peak_timing(show=False)
    z.show_cut_contributions(flat=True, show=False)
    z.compare_consecutive_cuts(show=False, save_single=False)
    z.compare_consecutive_cuts(scale=True, show=False, save_single=False)
    z.draw_waveforms(20, show=False)
    z.Pulser.draw_waveforms(20, show=False)
    z.draw_waveforms(5000, show=False)
    z.draw_intlength_vs_triggercell(show=False)
    z.draw_peak_integrals(show=False, event=120040)
    z.draw_cut_means(show=False)
    z.draw_intdiff_vs_triggercell(show=False)
    z.draw_signal_map(draw_option='colz', show=False)
    z.draw_signal_vs_signale(show=False)
    z.draw_ped_sigma_selection(show=False)
