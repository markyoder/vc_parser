'''
# script to create figures and other things for VC El Mayor Cucapah and Napa papers and posters. 
# we want this to feel just like being inside vc_parser and/or vc_geodetic. we'll want to pull some 
# functions, etc. directly from there maybe, but mostly we'll be calling from those scripts in very specific, not
# portable ways. in the end, we want to push a button and have all the datasets, figs, etc. come out (mabye with a bit of layering).
'''
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
	section_label_adjustments = {116:{'dx':0., 'dy':-2.}, 110:{'dx':-3., 'dy':3.}, 107:{'dx':-1, 'dy':3.}, 149:{'dx':0., 'dy':3.}, 113:{'dx':-3., 'dy':1.}}
		
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
	catalog, sweep_index = vc_ETAS_catalog(section_ids=sectionses[map_flavor], sim_file=default_sim_file, start_year=10000., end_year=None, n_cpus=None, fignum=0, map_size=[10,8], etas_mc=3.0, sweeps_index=None)
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

	if sections==None:
		sections = vc_parser.get_nearest_section_ids(lat_0=32.128, lon_0=-115.303, section_ids=vc_parser.emc_sections, dist_mode=0, n_sections=n_sections, verbose=False)
	print "sections: ", sections
	#
	#################
	q=waiting_time_figs(section_ids=sections, file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], keep_figs=False, output_dir='emc_WT_figs_n_%d' % n_sections, mc_nits=100000, n_cpus=None)
	#
	EWTf = EMC_EWT_figs(section_ids=sections, m0=m0, fits_data_file_CDF='CDF_EMC_figs/VC_CDF_Weibull_fits_dump.npy', WT_catalog_format='data/VC_CFF_timeseries_section_%d.npy', sim_file=default_sim_file, n_t0=10000, fnum=0, output_dir='expected_emc_waiting_time_figs_n_%d' % n_sections, do_local_fit=False)
	#
	#
	b=vc_parser.expected_waiting_time_t0(section_ids=sections, do_local_fit=False)
	plt.xlabel('Time since last $m>%.2f$ event, $t_0$ (years)' % m0)
	plt.ylabel('Expected interval $\Delta t$')
	plt.title('EMC expected waiting times for fault sections:\n' + ', '.join(map(str, sections)) + '\n')
	plt.savefig('expected_emc_waiting_time_figs_n_%d/emc_ewt_n%d_composite.png' % (n_sections, n_sections))
	
	return EWTf

def EMC_WT_dist_Appendix(wt_dir='VC_WT_probs_20150109', output_file = 'VC_WT_probs_20150109/appendix_wt_probs.tex', section_ids=vc_parser.emc_sections):
	'''
	# make at least the framework for an appendix of all the WT distribution figures.
	# file_format='VC_CDF_WT_m70_section_%d.png'
	#
	'''
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


def EMC_exp_WT_Appendix(wt_dir='expected_waiting_time_figs', output_file = 'expected_waiting_time_figs/appendix_exp_wt.tex', section_ids=vc_parser.emc_sections):
	'''
	# make at least the framework for an appendix of all the WT distribution figures.
	# file_format='VC_CDF_WT_m70_section_%d.png'
	#
	'''
	#
	src_template_file = '/home/myoder/Dropbox/Research/VC/EMC_VC_paper/VC_forecasting_yoder.tex'
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
		f.write('\\section{Appendix B: Expected waiting-times for fault segments in the El Mayor-Cucapah region}\n')
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
		
