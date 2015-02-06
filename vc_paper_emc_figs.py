'''
# script to create figures and other things for VC El Mayor Cucapah and Napa papers and posters. 
# we want this to feel just like being inside vc_parser and/or vc_geodetic. we'll want to pull some 
# functions, etc. directly from there maybe, but mostly we'll be calling from those scripts in very specific, not
# portable ways. in the end, we want to push a button and have all the datasets, figs, etc. come out (mabye with a bit of layering).
'''
# to run headless, use this:
#import matplotlib
#matplotlib.use('Agg')
#import pylab as plt

from vc_parser import *
from vc_geodetic import *
import vc_geodetic
import vc_parser

sectionses = {'napa':napa_sections, 'emc':emc_sections}
#emc_event = {'lat': 32.128, 'lon':-115.303, 'mag':7.2, 'event_time':dtm.datetime(2010,4,4,3,40,41,tzinfo=pytz.timezone('US/Pacific'))}

def vc_map_with_fault_sections(map_flavor='napa', i_min=1000, i_max=4000, etas_gridsize=.1, f_out=None, plot_quake_dots=False, plot_section_labels=None, fault_colors=None, plot_sections=[], verbose=True, sim_file=default_sim_file, map_size=[8.,10.], map_res='i', map_padding = .7, n_cpus=None, fignum=0):
	# make an etas map from a vc catalog. we probably have to pare down the catalog, since simulated catalogs are super huge.
	# sectionses: a list of a list of sections[ [secs1], [secs2], etc.] to (over)plot separately.
	#
	# we'll want to control the color cycle:
	if fault_colors==None: fault_colors = mpl.rcParams['axes.color_cycle']
	if isinstance(fault_colors, str): fault_colors = [fault_colors]
	#colors_ =  mpl.rcParams['axes.color_cycle']		# use like: this_color = colors_[j%len(colors_)]
	# ... which is a list, like: ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	#
	# hard-code some section label position adjustments, as a factor of label_delta_{x,y}:
	section_label_adjustments = { 56:{'dx':3., 'dy':1.}, 57:{'dx':2., 'dy':0.}, 105:{'dx':-5., 'dy':2.}, 107:{'dx':-6., 'dy':4.5}, 110:{'dx':-8., 'dy':3.}, 112:{'dx':1., 'dy':-1.},
	115: {'dx':4., 'dy':-.5}, 116:{'dx':3., 'dy':-4.}, 149:{'dx':0., 'dy':4.}, 113:{'dx':-3., 'dy':1.}, 56:{'dx':2., 'dy':-1.}, 57:{'dx':2., 'dy':1.}, 84:{'dx':-3., 'dy':0.}, 73:{'dx':10., 'dy':1.}, 17:{'dx':0., 'dy':-1.}  }
		
	#	
	#
	# extras:
	if verbose: print "customize with map_flavor: ", map_flavor
	if map_flavor.lower()=='napa':
		plt.title('Virtual Quake ETAS map of Napa region\n\n')
		mainshock = {'event_time':dtm.datetime(2014, 8, 24, 3+7, 20, 44, tzinfo=pytz.timezone('UTC')), 'lat':38.22, 'lon':-122.313, 'mag':6.0}
		plot_sections = vc_parser.napa_sections
		#		
	if map_flavor.lower() == 'emc':
		plt.title('Virtual Quake ETAS map of El Mayor-Cucapah region\n\n')
		mainshock = {'event_time':dtm.datetime(2010, 4, 10, 3+7, 40, 41, tzinfo=pytz.timezone('UTC')), 'lat':32.128, 'lon':-115.303, 'mag':7.2}
		#
		# for now, hard code this. we'll do conditional bits later:
		plot_sections = vc_parser.emc_sections
		#
		#if plot_sections!=None and len(plot_sections)==0:
		#	plot_sections = [vc_parser.get_nearest_section_ids(n_sections=5, lat_0=mainshock['lat'], lon_0=mainshock['lon'], verbose=False)]
		#
	#
	if verbose: print "set up section_labels: "
	# section labels:
	if plot_section_labels==None and plot_sections!=None:
		plot_section_labels = [str(x) for x in plot_sections]
	while len(plot_section_labels)<len(plot_sections): plot_section_labels += [str(plot_sections[len(plot_section_labels)])]
	# now, we should have equal length sections and section_labels for plotting...
	#
	if verbose: print "section labels set up. get vc_ETAS_catalog()."
	# first, get a catalog (we can do this inline wiht the map script, but this provides some added flexibility):
	
	#catalog, sweep_index = vc_ETAS_catalog(section_ids=sectionses[map_flavor], sim_file=default_sim_file, start_year=10000., end_year=None, n_cpus=None, fignum=0, map_size=[10,8], etas_mc=3.0, sweeps_index=None)
	# now, make a map using part of this catalog:
	#
	#
	if verbose: print "... and now, make seismicity_map(), then do some local plotting."
	#
	#my_map = vc_parser.seismicity_map(section_ids=sectionses[map_flavor], start_date=10000., etas_catalog={'catalog':list(catalog[i_min:i_max]), 'sweep_index':sweep_index}, p_map=0., etas_gridsize=etas_gridsize)
	#
	ll_range = vc_parser.get_fault_model_extents(section_ids=plot_sections, sim_file=sim_file, n_cpus=n_cpus)
	print "ll_range: ", ll_range
	#
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	#
	plt.figure(fignum)
	plt.ion()
	plt.clf()
	my_map = vc_parser.vc_basemap( projection='cyl', llcrnrlon=ll_range['lon_min']-map_padding, llcrnrlat=ll_range['lat_min']-map_padding, urcrnrlon=ll_range['lon_max']+map_padding, urcrnrlat=ll_range['lat_max']+map_padding, lon_0=lon_0, lat_0=lat_0, resolution=map_res)
	
	#
	#
	ms_lon, ms_lat = my_map(mainshock['lon'], mainshock['lat'])
	#quakes_x, quakes_y = my_map(catalog['lon'][i_min:i_max], catalog['lat'][i_min:i_max])
	ln_ms1 = plt.plot(ms_lon, ms_lat, 'k*', ms=18, zorder=6, alpha=.9)
	ln_ms2 = plt.plot(ms_lon, ms_lat, 'r*', ms=15, zorder=7, alpha=.9)
	#
	# and later on, we might need some idea of the spatial extents of the map. for now, use the spread on quake-dots:
	#x_min, x_max = min(quakes_x), max(quakes_x)
	#y_min, y_max = min(quakes_y), max(quakes_y)
	x_min, y_min = my_map(ll_range['lon_min'], ll_range['lat_min'])
	x_max, y_max = my_map(ll_range['lon_max'], ll_range['lat_max'])
	label_delta_x = (x_max - x_min)*.0005
	label_delta_y = (y_max - - y_min) * .0005		# just taking a stab at this scaling for now...
	#
	if plot_quake_dots: ln_quakes = plt.plot(quakes_x, quakes_y, 'm.', zorder=5, alpha=.6)
	#
	if verbose: print "now, plot section traces:"
	#
	# now, plot any special sections:
	for j_sections, sections in enumerate(plot_sections):
		#
		if verbose: print "plotting sections entry: ", sections, plot_section_labels[j_sections]
		#
		# get fault traces for (list of) section_id(s):
		ft = vc_parser.get_block_traces(section_ids=sections, fignum=None, lat_lon=True, sim_file=default_sim_file)
		this_color = fault_colors[j_sections%len(fault_colors)]	# note: if we provide only one color, this will always produce that entry.
		#
		section_x_min = None
		section_x_max = None
		section_y_min = None
		section_y_max = None
		for j, rw in enumerate(ft):
			X,Y,Z = zip(*rw)
			#print 'min-min: ', min(X), min(Y), min(section_x_min, min(X)), min(section_y_min, min(Y))
			if section_x_min == None: section_x_min = min(X)
			if section_y_min == None: section_y_min = min(Y)	# note: max(None, 1) = 1; min(None,1) = None (or something like it).
			section_x_min = min(section_x_min, min(X))
			section_y_min = min(section_y_min, min(Y))
			
			section_x_max = max(section_x_max, max(X))
			section_y_max = max(section_y_max, max(Y))
			#
			X,Y = my_map(X,Y)
			plt.plot(X,Y, '-%s' % this_color, lw=1.5)
			#
			# label?:
		if j_sections<len(plot_section_labels):
			#x = numpy.mean(X)
			#y = numpy.mean(Y)
			#			
			x = numpy.mean([section_x_min, section_x_max])
			y = numpy.mean([section_y_min, section_y_max])
			deltas = section_label_adjustments.get(sections, {'dx':0., 'dy':0.})
			dx = deltas['dx']
			dy = deltas['dy']
			#
			plt.text(x + label_delta_x*(1.0 + dx), y + label_delta_y*(1.0 + dy), '%s' % plot_section_labels[j_sections], color=this_color)
	#
	if f_out!=None:
		plt.savefig(f_out)
	#
	return my_map

def vc_etas_RI_map(map_flavor='napa', i_min=1000, i_max=4000, etas_gridsize=.1, plot_sections=[], f_out=None, plot_quake_dots=False, fault_colors = None):
	# note, here "RI" means "Relative Intensity", NOT "recurrence interval"
	# make an etas map from a vc catalog. we probably have to pare down the catalog, since simulated catalogs are super huge.
	# sectionses: a list of a list of sections[ [secs1], [secs2], etc.] to (over)plot separately.
	#
	# we'll want to control the color cycle:
	#colors_ =  mpl.rcParams['axes.color_cycle']		# use like: this_color = colors_[j%len(colors_)]
	# ... which is a list, like: ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	if fault_colors==None: fault_colors = mpl.rcParams['axes.color_cycle']	
	if isinstance(fault_colors, str): fault_colors = [fault_colors]
	#	
	# first, get a catalog (we can do this inline wiht the map script, but this provides some added flexibility):
	section_ids = sectionses.get(map_flavor, None)
	if section_ids==None:
		with vc_parser.h5py.File(default_sim_file) as vc_data:
			print "getting all section_ids from %s" % default_sim_file
			section_ids = list(set(vc_data['block_info_table']['section_id']))
	
	#catalog, sweep_index = vc_ETAS_catalog(section_ids=sectionses[map_flavor], sim_file=default_sim_file, start_year=10000., end_year=None, n_cpus=None, fignum=0, map_size=[10,8], etas_mc=3.0, sweeps_index=None)
	catalog, sweep_index = vc_ETAS_catalog(section_ids=section_ids, sim_file=default_sim_file, start_year=10000., end_year=None, n_cpus=None, fignum=0, map_size=[10,8], etas_mc=3.0, sweeps_index=None)
	#
	# now, make a map using part of this catalog:
	#
	my_map = vc_parser.seismicity_map(section_ids=sectionses[map_flavor], start_date=10000., etas_catalog={'catalog':list(catalog[i_min:i_max]), 'sweep_index':sweep_index}, p_map=0., etas_gridsize=etas_gridsize)
	#
	# extras:
	if map_flavor.lower()=='napa':
		plt.title('Virtual Quake ETAS map of Napa region\n\n')
		mainshock = {'event_time':dtm.datetime(2014, 8, 24, 3+7, 20, 44, tzinfo=pytz.timezone('UTC')), 'lat':38.22, 'lon':-122.313, 'mag':6.0}
		#		
	if map_flavor.lower() == 'emc':
		plt.title('Virtual Quake ETAS map of El Mayor-Cucapah region\n\n')
		mainshock = {'event_time':dtm.datetime(2010, 4, 10, 3+7, 40, 41, tzinfo=pytz.timezone('UTC')), 'lat':32.128, 'lon':-115.303, 'mag':7.2}
		if plot_sections!=None and len(plot_sections)==0:
			plot_sections = [vc_parser.get_nearest_section_ids(n_sections=5, lat_0=mainshock['lat'], lon_0=mainshock['lon'], verbose=False)]
		#
	#if map_falvor.lower() == 'alcal':
	#	plt.title('Virtual Quake ETAS map of California\n\n')
	#	# for now, use the emc mainshock:
	#	mainshock = {'event_time':dtm.datetime(2010, 4, 10, 3+7, 40, 41, tzinfo=pytz.timezone('UTC')), 'lat':32.128, 'lon':-115.303, 'mag':7.2}
	#	# get all fault sections:
	#	with vc_parser.h5py.File(default_sim_file) as vc_data:
	#		section_ids = list(set(vc_data['block_info_table']['section_id']))
	#
	ms_lon, ms_lat = my_map.cm(mainshock['lon'], mainshock['lat'])
	quakes_x, quakes_y = my_map.cm(catalog['lon'][i_min:i_max], catalog['lat'][i_min:i_max])
	ln_ms1 = plt.plot(ms_lon, ms_lat, 'k*', ms=18, zorder=6, alpha=.9)
	ln_ms2 = plt.plot(ms_lon, ms_lat, 'r*', ms=15, zorder=7, alpha=.9)
	#
	if plot_quake_dots: ln_quakes = plt.plot(quakes_x, quakes_y, 'm.', zorder=5, alpha=.6)
	#
	# now, plot any special sections:
	for sections in plot_sections:
		#
		# get fault traces for (list of) section_id(s):
		ft = vc_parser.get_fault_traces(section_ids=sections, fignum=None, lat_lon=True, sim_file=default_sim_file)
		# break ft into sets of lines. there will be one dict-entry per fault (so we might get less than n_sections values).
		# we'll need to set up an iterator for colors as well. look around for a variable called _colors or colors_ (or something like
		# that) for instructions...
		#
		# each entry in ft will be a list of list of vectors [ v1, v2, ...] -->
		# [
		#  [
		#   [x00, y00], [x01,y01]
		#   ]
		#  [ [x10, y10], [x11, y11]
		#   ]
		#  [ ...etc...
		#   ]
		#  ]
		#
		for j, (fault_id, fault) in enumerate(ft.iteritems()):
			#this_color = 'r'
			this_color = fault_colors[j*len(fault_colors)]	# note: if we provide only one color, this will always produce that entry.
			X,Y=[],[]
			for vec in fault:
				for compon in vec:
					X+=[compon[0]]
					Y+=[compon[1]]
				#
			#
			X,Y = my_map.cm(X,Y)
			plt.plot(X,Y, '-%s' % this_color, lw=1.5)
			#
			# label?:
			if plot_section_labels:
				x = numpy.mean(X)
				y = numpy.mean(Y)
	#
	if f_out!=None:
		plt.savefig(f_out)
	#
	return my_map

def vc_surface_deformation(map_flavor = 'napa', map_resolution='i'):
	#	
	#blockwise_obj = blockwise_slip(sim_file=default_sim_file, faults=None, sections=None, pipe=None, f_pickle_out='dumps/blockwise_slip_0.pkl', plot_factor=1.0)
	# this takes a while, so a pre-calculated version might be preferable:
	#blockwise_obj = blockwise_slip(sim_file=default_sim_file, sections=napa_sections, f_pickle_out=None, plot_factor=1.0)
	#
	# pre-calculated:
	blockwise_obj='dumps/blockwise_slip_Nov11.pkl'
	disp_field = vc_geodetic.slip_field(blockwise_obj='dumps/blockwise_slip_Nov11.pkl', section_ids=sectionses[map_flavor], dx=1.5, dy=1.5, i_start=0, i_stop=-1, f_out='temp/slip_field.pkl')	# dx, dy are factors of fault-length, so \Delta_x = block_x_len * dx
	#
	zz=plot_disp_vector_field(okada_disps=disp_field, fignum=0, plot_factor=50., section_ids=sectionses[map_flavor], map_resolution=map_resolution)
	#
	# and now decorate it:
#
#def vc_etas_RI_map_emc(gridsize=.1):
#	# the final script for the EMC etas map...
#	#
#	#secs = vc_parser.get_nearest_section_ids(n_sections=5, lat=emc_event['lat'], lon=emc_event['lon'])
#	#
#	return vc_etas_RI_map(map_flavor='emc', i_min=1000, i_max=4000, etas_gridsize=.1, plot_sections=[secs])
#	#
#
def closest_emc_fault_section_figs(n_sections=5, sections=None, m0=7.0):
	'''
	# a few plots based on the nearest sections to the EMC event.
	# sections: overrides n_sections? maybe we'll truncate. anyway, you can provide sections or leave it None and fetch them.
	'''
	
	output_dir = 'figs_gji/emc_RI_WT_figs_n_%d' % n_sections

	if sections==None:
		sections = vc_parser.get_nearest_section_ids(lat_0=32.128, lon_0=-115.303, section_ids=vc_parser.emc_sections, dist_mode=0, n_sections=n_sections, verbose=False)
	print "sections: ", sections
	#
	#################
	#q=waiting_time_figs(section_ids=sections, file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], keep_figs=False, output_dir='emc_WT_figs_n_%d' % n_sections, mc_nits=100000, n_cpus=None)
	#
	# ... but we should already have these from gji_RI_probabilities(), and we don't use the individual figures anyway...
	#q1=conditional_RI_figs(section_ids=sections,file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=output_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None)
	
	q1=conditional_RI_figs(section_ids=[sections],file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=output_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None)
	plt.title('CDF $P(\Delta t_r)$ for 5 nearest fault segments (%s)' % (', '.join([str(x) for x in sections])))
	#ax = plt.gca()
	#ax.set_xlim(0., 175.)
	#plt.draw() 
	#
	plt.savefig('%s/EMC_conditional_ri_aggregate_%d_nearest.png' % (output_dir, n_sections))
	
	#
	EWTf = EMC_EWT_figs(section_ids=sections, m0=m0, fits_data_file_CDF='CDF_EMC_figs/VC_CDF_Weibull_fits_dump.npy', WT_catalog_format='data/VC_CFF_timeseries_section_%d.npy', sim_file=default_sim_file, n_t0=10000, fnum=0, output_dir=output_dir, do_local_fit=False)
	#
	#
	b=vc_parser.expected_waiting_time_t0(section_ids=sections, do_local_fit=False)
	plt.xlabel('Time since last $m>%.2f$ event, $t_0$ (years)' % m0)
	plt.ylabel('Expected interval $\Delta t$')
	plt.title('EMC expected waiting times for fault sections:\n' + ', '.join(map(str, sections)) + '\n')
	plt.savefig('%s/emc_ewt_n_%d_composite.png' % (output_dir, n_sections))
	
	return EWTf

#
def gji_RI_probabilities(section_ids=vc_parser.emc_sections, output_dir='figs_gji', production_gen_figs='/home/myoder/Dropbox/Research/VC/VC_EMC_gji_yoder/general_figs/'):
	#
	# ... and we should write a script to plot pre-fit data. it's probably better to just write a separate script, since matching the 
	# existing t_0 values with the parent (fit + fig generating script) is tricky (since they're floats). the better approach is to separate
	# the fitting and plotting parts... so do this in PyVQ.
	#
	gji_pri_figs_dir = '%s/pri' % output_dir
	#
	# get individual section figures:
	z=conditional_RI_figs(section_ids=section_ids,file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=gji_pri_figs_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None)
	#
	# now get the aggregate figure:
	z=conditional_RI_figs(section_ids=[section_ids],file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=gji_pri_figs_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None)
	plt.title('CDF $P(\Delta t_r)$ for EMC region, aggregated')
	ax = plt.gca()
	ax.set_xlim(0., 175.)
	plt.draw()
	#
	aggregate_file_name = '%s/EMC_conditional_ri_aggregate.png' % output_dir
	plt.savefig(aggregate_file_name)
	#
	# now, copy the primary figures to the production environment:
	for fname in ['RI_conditional_CDF_m70_section_%d' % s_id for s_id in (111, 149, 123, 124, 30)] + [aggregate_file_name]:
		# and eventualy, finish this...
		pass

def EMC_WT_dist_Appendix(wt_dir='figs_gji/pri', output_file = 'figs_gji/pri/appendix_wt_probs.tex', section_ids=vc_parser.emc_sections):
	'''
	# make at least the framework for an appendix of all the WT distribution figures.
	# file_format='VC_CDF_WT_m70_section_%d.png'
	#
	'''
	# formerly: wt_dir='VC_WT_probs_20150109'
	#
	src_template_file = '/home/myoder/Dropbox/Research/VC/EMC_VC_paper/VC_forecasting_yoder.tex'
	header_string = ''
	with open(src_template_file,'r') as f:
		for rw in f:
			if rw.startswith('\\begin{abstract}'): break
			header_string += rw
		
	#
	with open(output_file, 'w') as f:
		f.write('% waiting times appendix, generated by vc_paper_emc_figs.EMC_WT_dist_Appendix()\n\n\n')
		#
		f.write(header_string)
		#
		f.write('\n\n')
		f.write('\\section{Appendix A: Waiting-time probability distributions for fault segments in the El Mayor-Cucapah region}\n')
		#
		for i, section_id in enumerate(section_ids):
			f.write('\\begin{figure}\n')
			f.write('\\noindent \\textbf{a)}\\includegraphics[width=2.5in]{VC_CDF_WT_m70_section_%d.png}\n' % section_id)
			f.write('\\textbf{b)}\\includegraphics[width=2.5in]{../expected_waiting_time_figs/EWT_m0_70_section_%d}\n' % section_id)
			f.write('\\caption{\\textbf{a)} Waiting time distributions and \\textbf{b)} Expected waiting times for fault section %d.}\n' % section_id)
			f.write('\\label{fig:appendix_wt_%d}\n\\end{figure}\n' % section_id)
			f.write('\n\n')
			#
			if i%3==0: f.write('\\clearpage\n\n')
			#
		f.write('\n\n')
		f.write('\\end{document}')
		#
	#


#def EMC_exp_WT_Appendix(wt_dir='expected_waiting_time_figs', wt_prob_dir='VC_WT_probs_20150109', output_file = 'appendix_exp_wt.tex', section_ids=vc_parser.emc_sections):
def EMC_exp_WT_Appendix(wt_dir='\WTexpfigs', wt_prob_dir='\WTPfigs', output_file = 'vc_appendix_emc_PRI.tex', section_ids=vc_parser.emc_sections):
	'''
	# make at least the framework for an appendix of all the WT distribution figures.
	# file_format='VC_CDF_WT_m70_section_%d.png'
	#
	'''
	#
	manuscript_folder = '/home/myoder/Dropbox/Research/VC/VC_EMC_gji_yoder'
	#
	src_template_file = '%s/VC_forecasting_yoder.tex' % manuscript_folder
	header_string = ''
	with open(src_template_file,'r') as f:
		for rw in f:
			if rw.startswith('\\begin{abstract}'): break
			header_string += rw
		
	#
	with open(output_file, 'w') as f:
		f.write('% waiting times appendix B, generated by vc_paper_emc_figs.EMC_exp_WT_Appendix()\n\n\n')
		#
		f.write(header_string)
		#
		f.write('\n\n')
		f.write('\\textbf{Appendix A: Waiting-time distributions and expected waiting-times for fault segments in the El Mayor-Cucapah region}\n')
		#f.write('\\clearpage\n\n')
		#
		for i, section_id in enumerate(section_ids):
			#
			if i%3==0: f.write('\\clearpage\n\n')
			#
			f.write('\\begin{figure}\n')
			#f.write('\\noindent \\textbf{a)}\\includegraphics[width=2.5in]{VC_CDF_WT_m70_section_%d.png}\n' % section_id)
			f.write('\\noindent \\textbf{a)}\\includegraphics[width=2.5in]{%s/VC_CDF_WT_m70_section_%d.png}\n' % (wt_prob_dir, section_id))
			#f.write('\\textbf{b)}\\includegraphics[width=2.5in]{../expected_waiting_time_figs/EWT_m0_70_section_%d}\n' % section_id)
			f.write('\\textbf{b)}\\includegraphics[width=2.5in]{%s/EWT_m0_70_section_%d}\n' % (wt_dir, section_id))
			f.write('\\caption{\\textbf{a)} Waiting time distributions and \\textbf{b)} Expected waiting times for fault section %d.}\n' % section_id)
			f.write('\\label{fig:appendix_pri_%d}\n\\end{figure}\n' % section_id)
			f.write('\n\n')
			
			#
		f.write('\n\n')
		f.write('\\end{document}')
		#
	#
	#os.system('cp %s %s' % (output_file, '%s/%s' % (manuscript_folder, output_file)))
	os.system('cp %s %s/%s' % (output_file, manuscript_folder, output_file))


def EMC_EWT_figs(section_ids=None, m0=7.0, fits_data_file_CDF='CDF_EMC_figs/VC_CDF_Weibull_fits_dump.npy', WT_catalog_format='data/VC_CFF_timeseries_section_%d.npy', sim_file=default_sim_file, n_t0=10000, fnum=0, output_dir='expected_waiting_time_figs', do_local_fit=False):
	# make a set of "Expected Waiting Time" figures.
	# (wrapper for expected_waiting_time_t0() )
	#
	if os.path.isdir(output_dir)==False:
		os.mkdir(output_dir)
	#
	if isinstance(section_ids, str):
		if section_ids.upper()=='EMC':
			section_ids = [{'EMC':list(napa_sections)}] + list(napa_sections)
		if section_ids.lower() == 'napa':
			section_ids = [{'Napa':list(napa_sections)}] + list(napa_sections)

	#
	if section_ids==None:
		# get 'em all:
		with h5py.File(sim_file) as vcdata:
			section_ids=list(set(vcdata['block_info_table']['section_id']))
		#
	#
	# now, for each section_id, get EWT and draw a picture.
	for sec_id in section_ids:
		# first, is sec_id a number or a list? expected_waiting_time_t0() expects a list of section_ids.
		name_str = None
		#
		# handle the input a bit. we can provide just a list of section ids or section id lists or we can mix in some dicts:
		#  [1,2,3, [4,5,6], 7, [8,9,10], {'teens':[14,15,16,17,18,19]}, ... ]
		# 
		# a list of lists (or mi
		if isinstance(sec_id, dict):
			# we've been given a dictionary and some instructions. 
			# assume for now that it's a single pair.
			for key,val in sec_id.iteritems():
				# just keep the last one...
				name_str = str(key)
				sec_id=val
		if not hasattr(sec_id, '__len__'): sec_id=[sec_id]
		if name_str == None: 
			#name_str=str(*sec_id)
			name_str = '_'.join([str(x) for x in sec_id])
		#
		EWT = expected_waiting_time_t0(section_ids=sec_id, m0=m0, fits_data_file_CDF=fits_data_file_CDF, WT_catalog_format=WT_catalog_format, sim_file=sim_file, n_t0=n_t0, fnum=fnum, do_local_fit=do_local_fit)
		#
		plt.figure(fnum)
		plt.title('Expected Waiting times: sections %s' % name_str)
		plt.xlabel('Time since last $m>%.2f$ event, $t_0$' % m0)
		plt.ylabel('Expected interval $\Delta t$')
		#
		plt.savefig('%s/EWT_m0_%s_section_%s.png' % (output_dir, str(m0).replace('.',''), name_str))

def create_ROC_figs_GT_data(section_ids = vc_parser.emc_sections, nits=2500, fnum=0, num_roc_points=100, output_dir = 'dumps/gji_roc_gt_detail', m0=7.0):
	'''
	# create a whole slew of ROC data. this will include the optimized "best fit" (using whatever metric) and also the raw, full MC output.
	# (note we could do this with simple_mpp_optimizer(), but for now let's just run some modest size figs (say nits=2500) or so
	# for proof of concept.
	#
	#call like: simple_metric_optimizer(CFF=None, m0=7.0, b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=1000, keep_set=False, set_name='data_set', dump_file=None, f_gt_lt=operator.gt, f_score=operator.div, section_id=16)
	#
	# note also that we can provide section_id lists like section_id=[1,2,3], and simple_metric_optimizer() will use combine_section_CFFs()
	# to assemble an aggregate catalog.
	#
	'''
	#
	#
	R = random.Random()
	#
	for sec_id in section_ids:
		if isinstance(sec_id, float): sec_id=int(sec_id)
		if isinstance(sec_id, int): sec_id=[sec_id]
		#
		sec_str = '_'.join([str(x) for x in sec_id])
		#
		f_output_name = '%s/roc_sec_%s_nits_%d.npy' % (output_dir, sec_str, nits)
		f_output_name_rand = '%s/roc_sec_rand_%s_nits_%d.npy' % (output_dir, sec_str, nits)
		#
		CFF = vc_parser.combine_section_CFFs(sec_id)
		opt_datas, raw_datas = simple_metric_optimizer(CFF=CFF, m0=m0, b_min=-.25, b_max=.25, d_b=.01, nyquist_min=.2, nyquist_max=1.5, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name, f_gt_lt=operator.gt, f_score=operator.div, section_id=sec_id)
		#
		# random forecast: literally chooses alert segments at random. lines up really nicely on H=F.
		roc_random = vc_parser.get_random_forecast_set(CFF=CFF, section_id=sec_id, m0=m0, P_min=0., P_max=1.0, set_name=None, nits=nits, format='recarray', do_roc_plot=False, fignum=None)
		#
		# now, randomize this catalog and plot to show contrast:
		# this randomized the order of an existing catalog: 1) extract intervals, then shuffle them, 2) shuffle the catalog, 3) reassign event_year
		# as event_year[i-1]+interval[i] (basically). SO, the distribution of intervals is maintained.
		#CFF_rand = vc_parser.randomize_CFF(section_id=sec_id)
		#
		# this draws intervals from a uniform distribution and assigns t_i = t_{i-1} + <dt>*random().
		#CFF_rand = vc_parser.random_like_CFF(CFF_in=None, section_id=sec_id)
		
		#
		#opt_datas_rand, raw_datas_rand = simple_metric_optimizer(CFF=CFF_rand, m0=7.0, b_min=-.25, b_max=.25, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name_rand, f_gt_lt=operator.gt, f_score=operator.div, section_id=None)
		
		#return raw_datas
		# for now, put this here (get it done the first time). later, we'll move this off-line and use pre-calculated data, etc.
		plotted = plot_section_ROC_curve(roc_data=raw_datas, section_id=None, fignum=fnum, num_points=num_roc_points, label_str='section_id: %s' % sec_str, markers='.-')
		
		#plotted_rand = plot_section_ROC_curve(roc_data=raw_datas_rand, section_id=None, fignum=fnum, do_clf=False, num_points=num_roc_points, label_str='(random)', markers='.')
		#
		plt.figure(fnum)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		plt.title('Optimal ROC for Section %s' % sec_str)
		plt.savefig('%s/roc_opt_sec_%s_nits_%d.png' % (output_dir, sec_str, nits))
		#
		plt.figure(fnum+1)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		plt.title('Raw ROC for Section %s' % sec_str)
		plt.savefig('%s/roc_raw_sec_%s_nits_%d.png' % (output_dir, sec_str, nits))
#
def create_ROC_figs_LT_data(section_ids = vc_parser.emc_sections, nits=2500, fnum=0, num_roc_points=100, output_dir = 'dumps/gji_roc_lt_detail', m0=7.0):
	'''
	# for LT metric (acceleration): 
	#create a whole slew of ROC data. this will include the optimized "best fit" (using whatever metric) and also the raw, full MC output.
	# (note we could do this with simple_mpp_optimizer(), but for now let's just run some modest size figs (say nits=2500) or so
	# for proof of concept.
	#
	#call like: simple_metric_optimizer(CFF=None, m0=7.0, b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=1000, keep_set=False, set_name='data_set', dump_file=None, f_gt_lt=operator.gt, f_score=operator.div, section_id=16)
	#
	# note also that we can provide section_id lists like section_id=[1,2,3], and simple_metric_optimizer() will use combine_section_CFFs()
	# to assemble an aggregate catalog.
	#
	'''
	#
	R=random.Random()
	for sec_id in section_ids:
		if isinstance(sec_id, float): sec_id=int(sec_id)
		if isinstance(sec_id, int): sec_id=[sec_id]
		#
		sec_str = '_'.join([str(x) for x in sec_id])
		#
		f_output_name = '%s/roc_sec_lt_%s_nits_%d.npy' % (output_dir, sec_str, nits)
		f_output_name_rand = '%s/roc_sec_lt_rand_%s_nits_%d.npy' % (output_dir, sec_str, nits)
		#
		CFF = vc_parser.combine_section_CFFs(sec_id)
		#
		opt_datas, raw_datas = simple_metric_optimizer(CFF=CFF, section_id=sec_id, m0=m0, b_min=-.25, b_max=.25, d_b=.01, nyquist_min=.2, nyquist_max=1.2, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name, f_gt_lt=operator.lt, f_score=operator.div)
		#
		roc_random = vc_parser.get_random_forecast_set(CFF=CFF, section_id=sec_id, m0=m0, P_min=0., P_max=1.0, set_name=None, nits=nits, format='recarray', do_roc_plot=False, fignum=None)
		#
		# now, randomize this catalog and plot to show contrast:
		#CFF_rand = vc_parser.random_like_CFF(CFF_in=None, section_id=sec_id)
		#
		#opt_datas_rand, raw_datas_rand = simple_metric_optimizer(CFF=CFF_rand, m0=7.0, b_min=-.25, b_max=.25, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name_rand, f_gt_lt=operator.lt, f_score=operator.div, section_id=None)
		
		
		#return raw_datas
		# for now, put this here (get it done the first time). later, we'll move this off-line and use pre-calculated data, etc.
		plotted =      plot_section_ROC_curve(roc_data=raw_datas, section_id=None, fignum=fnum, num_points=num_roc_points, markers='.-')
		#plotted_rand = plot_section_ROC_curve(roc_data=raw_datas_rand, section_id=None, fignum=fnum, do_clf=False, num_points=num_roc_points, label_str='(random)',markers='.')
		plt.figure(fnum)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		plt.title('Optimal ROC (PSA) for Section %s' % sec_str)
		plt.legend(loc=0, numpoints=1)
		plt.savefig('%s/roc_psa_lt_opt_sec_%s_nits_%d.png' % (output_dir, sec_str, nits))
		#
		plt.figure(fnum+1)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		plt.title('Raw ROC (PSA) for Section %s' % sec_str)
		plt.legend(loc=0, numpoints=1)
		plt.savefig('%s/roc_psa_lt_raw_sec_%s_nits_%d.png' % (output_dir, sec_str, nits))	
###################################
#
# IAGS paper (short, letter bit for IAGS special publication):

#def iags_waiting_time_probabilities():
#	A=vc_parser.waiting_time_figs(section_ids=[24, 25], output_dir='figs_iags')
#	#
#	return A
def iags_RI_probabilities(section_ids=[123, 111, vc_parser.emc_sections], output_dir='figs_iags'):
	# use: def conditional_RI_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir='VC_CDF_WT_figs', mc_nits=100000, n_cpus=None, start_year=10000, end_year=None)
	#
	#for i, sec_id in enumerate(section_ids):
	
	z=conditional_RI_figs(section_ids=section_ids,file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=output_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None)

#
#def iags_expected_waiting_times(output_dir='figs_iags', section_ids = [123, 124, [123, 124], 125, [123, 124, 125]], m0=7.0, fnum=0):
def iags_expected_waiting_times(output_dir='figs_iags', section_ids = [123, 111, vc_parser.emc_sections], m0=7.0, fnum=0):
	#
	# make expected waiting time figures (yellow envelopes wiht black median line...)
	#
	plt.ion()
	#
	for sec_id in section_ids:
		# pass these as lists. if an integer (or string) is passed, wrap it as a list.
		if not hasattr(sec_id, 'append'): sec_id=[int(sec_id)]
		#
		#plt.figure(fnum)
		#plt.clf()
		B=vc_parser.expected_waiting_time_t0(section_ids=sec_id, m0=m0, do_local_fit=False, fnum=fnum)
		# now, annotate:
		plt.xlabel('Time $t_0$ since previous $m > %.2f$ earthquake' % m0)
		plt.ylabel('Expected interval $\\Delta t$ to next $m > %.2f$ earthquake' % m0)
		plt.title('Expected waiting times for section id(s): %s' % ', '.join([str(x) for x in sec_id]))
		#
		plt.savefig('%s/expected_waiting_time_m%s_sec_%s.png' % (output_dir, str(m0), '_'.join([str(x) for x in sec_id])))
		#

	
