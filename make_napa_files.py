import vc_parser

napa_sections = list(vc_parser.napa_region_section_filter['filter'])
# and these can just mix with the existing EMC data.
#A= vc_parser.get_EMC_CFF(sections=napa_sections, file_out_root='data/VC_CFF_timeseries_section')

# now, run some statistics:
for this_m in [6.0, 7.0]:
	magstr = str(this_m).replace('.','')
	weibull_recurr_output_dir = 'VC_Napa_m%s_recurr_figs' % magstr
	weibull_WT_output_dir = 'VC_Napa_m%s_WT_figs' % magstr
	
	
	#b = vc_parser.waiting_time_figs(section_ids=napa_sections, m0=this_m, output_dir=weibull_WT_output_dir)

	#b = vc_parser.recurrence_figs(section_ids=napa_sections, m0=this_m, output_dir=weibull_recurr_output_dir)
	
	output_dir='napa_EWT_m%s_time_figs' % magstr
	#b= vc_parser.EMC_EWT_figs(section_ids=napa_sections, m0=this_m, fits_data_file_CDF='%s/VC_CDF_Weibull_fits_dump.npy' % weibull_recurr_output_dir, WT_catalog_format='data/VC_CFF_timeseries_section_%d.npy', sim_file=vc_parser.default_sim_file, n_t0=10000, fnum=0, output_dir=output_dir, do_local_fit=False)
	#
	Z=vc_parser.simple_mpp_optimizer(sections='napa', section_names=None, start_year=10000., m0=this_m, b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=1000, keep_set=False, dump_file='dumps/napa_m%d_mpp_roc_optimize.npy' % this_m, n_cpus=None)
	
	try:
		ZZ=vc_parser.plot_best_opt_prams(scores_in=Z, plot_f_out='dumps/napa_roc_m%s.png' % str(this_m).replace('.', ''))
	except:
		ZZ=vc_parser.plot_best_opt_prams(scores_in='dumps/napa_m%d_mpp_roc_optimize.npy' % this_m, plot_f_out='dumps/napa_ros_m%s.png' % str(this_m).replace('.', ''))
	#

#b=vc_parser.waiting_time_figs(section_ids=napa_sections, m0=6.0, output_dir='VC_Napa_m6_WT_figs')

#b=vc_parser.recurrence_figs(section_ids=napa_sections, m0=6.0, output_dir='VC_Napa_m6_recurrence_figs')
