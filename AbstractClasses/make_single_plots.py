#!/usr/bin/python

from PadAnalysis import PadAnalysis
from Elementary import Elementary

tc = '201510'
a = Elementary(tc)
a.print_banner('STARTING SINGLE PLOT GENERATION')
a.print_testcampaign()


for run in [392, 398]:
    z = PadAnalysis(run, 0)
    z.save_dir = '{tc}_{dia}_{hv}_{rp}'.format(tc=z.TESTCAMPAIGN, dia=z.diamond_name, hv=z.bias, rp=z.run_number)
    z.draw_trigger_cell_vs_peakpos(show=False)
    z.draw_trigger_cell_vs_peakpos(show=False, corr=True)
    z.draw_bucket_pedestal(show=False)
    z.draw_peak_timing(show=False)
    z.show_cut_contributions(flat=True, show=False)
    z.compare_consecutive_cuts(show=False, save_single=False)
    z.compare_consecutive_cuts(scale=True, show=False, save_single=False)
    z.draw_waveforms(20, show=False)
    z.draw_waveforms(20, cut_string=z.Cut.generate_pulser_cut(), show=False)
    z.draw_waveforms(5000, show=False)
    z.draw_intlength_vs_triggercell(show=False)
    z.draw_peak_integrals(show=False, event=120040)
    z.draw_cut_means(show=False)
    z.draw_intdiff_vs_triggercell(show=False)
    z.draw_signal_map(draw_option='colz', show=False)

