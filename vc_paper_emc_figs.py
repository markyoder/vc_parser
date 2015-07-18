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
#
# some decorator stuff:
# set font sizes for figures:
fs_title  = 16
fs_label = 16.
fs_legend = 14.

def make_gji_figs():
	#
	# non-conditional CDFs:
	fs_title=16
	fs_legend=12
	fs_label=16

	# of course, choose output directories appropriately...
	
	# make the map(s) with fault sections:
	a=vc_map_with_fault_sections(map_flavor='emc', f_out='figs_gji/rev1/maps')
	# make the etas map (but only if necessary; this one takes a while...)
	#vc_etas_RI_map(map_flavor='napa', i_min=1000, i_max=10000, etas_gridsize=.1, plot_sections=[], f_out=None, plot_quake_dots=False, fault_colors = None):
	# the third map, i think, has to be mande manually using G.Earth.
	#
		#
	a=vc_parser.recurrence_figs(section_ids=vpp.vc_parser.emc_sections, output_dir='figs_gji/rev1/CDF_figs', fs_title=fs_title, fs_legend=fs_legend, fs_label=fs_label)
	z=vfp.EMC_EWT_figs(output_dir='figs_gji/rev1/EWTs', section_ids=vfp.emc_sections)
	waiting_time_aggregate(section_ids=vc_parser.emc_sections, do_local_fit=False, fs_title=fs_title, fs_legend=fs_legend, fs_label=fs_label, f_out='figs_gji/rev1/EWTs/EMC_expected_intervals.png')
	a=gji_RI_probabilities(section_ids=vc_parser.emc_sections, output_dir='figs_gji/rev1/pri', fs_legend=fs_legend, fs_label=fs_label)
	
	
	# it might be neccessary to re-make the appendix .tex:
	# EMC_WT_dist_Appendix(wt_dir='figs_gji/pri', output_file = 'figs_gji/pri/appendix_wt_probs.tex', section_ids=vc_parser.emc_sections)
	#
	# "closest 5" figures:
	a=closest_emc_fault_section_figs(output_dir='figs_gji/rev1/closest_sections/', fs_legend=fs_legend, fs_label=fs_label)
	#
	# surface deformation:
	# (note: 1) not fully save, etc. scripted, 2) we're pulling this fig from the paper since we don't really develop it.
	#vc_surface_deformation(map_flavor = 'napa', map_resolution='i')
	#
	#
	# create all the ROC data (note: use delta_b1 parameter to adjust the b_intiaite/b_maintain alert thresholds):
	#create_ROC_figs_LT_data(section_ids = vc_parser.emc_sections, nits=1000, fnum=0, num_roc_points=100, output_dir = 'dumps/gji_roc_lt_500', output_descriptions=[], fig_title_strs=[], m0=7.0, delta_b1=0.)
	# aggreagate:
	create_ROC_figs_LT_data(section_ids = [vc_parser.emc_sections], nits=1500, fnum=0, num_roc_points=100, output_dir = 'dumps/gji_roc_lt_500_b', output_descriptions=[], fig_title_strs=['EMCaggregate'], m0=7.0, delta_b1=0.)
	#
	# now, we need to create the ROC tabular data and the composite figures -- which should be straight forward, since the runner script (i think) compiles those data.
	#
	#best_fit_ROC_table({and get the calling parameters...})
	# best_fit_ROC_table(section_ids=vc_parser.emc_sections + ['EMC'], output_file='dumps/gji_roc_lt_500/roc_bestfits', m0=7.0, roc_file_format='dumps/gji_roc_lt_500/roc_sec_lt_%s_nits_1000_allprams.npy', longtable=True, lbl_str='tbl:roc_prams')
	z = best_fit_ROC_table(section_ids=vc_parser.emc_sections + ['EMC'], output_file='dumps/gji_roc_lt_500_b/roc_bestfits', m0=7.0, roc_file_format='dumps/gji_roc_lt_500_b/roc_sec_lt_%s_nits_1000_allprams.npy', longtable=True, lbl_str='tbl:roc_prams')
	#
	# section 16 ROC:
	z=roc_figure(section_id=16, title_str='Section 16 ROC', output_dir='figs_gji/rev1/ROCs_b/', fname=None, roc_data='dumps/gji_roc_lt_500_b/roc_sec_lt_16_nits_1000_allprams.npy')
	# full catalog ROC:
	a=roc_figure(roc_data='dumps/gji_roc_lt_500_b/roc_sec_lt_EMC_nits_1000_allprams.npy', title_str='EMC Aggregate', section_id=vfp.emc_sections, output_dir='figs_gji/rev1/ROCs_b/', fname='ROC_scatter_EMCsections.png')
	#
	rocs = plot_best_roc(n_rank=5, save_figs=True, b_0=0., nyquist_factor=.5, output_dir='figs_gji/rev1/ROCs', input_data_format='dumps/gji_roc_lt_500/roc_sec_lt_*_allprams*.npy')
	#
	z = gji_forecast_fig(fignum=0, section_id=16, f_out = 'figs_gji/rev1/forecast_section_16.png', time_range=(13400., 14850.), opt_data='dumps/gji_roc_lt_500_b/roc_sec_lt_%d_nits_1000_allprams.npy' )

def waiting_time_aggregate(section_ids=vc_parser.emc_sections, do_local_fit=False, fs_title=fs_title, fs_legend=fs_legend, fs_label=fs_label, f_out='figs_gji/rev1/EWTs/EMC_expected_intervals.png'):
	#
	b=vc_parser.expected_waiting_time_t0(section_ids=vc_parser.emc_sections, do_local_fit=False)

	plt.xlabel('Time since last $m>%.2f$ event, $t_0$ (years)' % 7.0, size=fs_label)
	plt.ylabel('Expected interval $\Delta t$', size=fs_label)
	plt.title('Expected waiting times for all EMC fault sections (aggregated)\n', size=fs_title)
	if not f_out==None:
		pth, fnm = os.path.split(f_out)
		if not os.path.isdir(pth): os.makedirs(pth)
		#
		plt.savefig(f_out)


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
	section_label_adjustments = {  17:{'dx':0., 'dy':-1.}, 56:{'dx':2., 'dy':-1.}, 57:{'dx':2., 'dy':1.}, 105:{'dx':-5., 'dy':2.}, 107:{'dx':-6., 'dy':4.5}, 110:{'dx':-8., 'dy':3.}, 111:{'dx':-1., 'dy':-1.}, 112:{'dx':1., 'dy':-1.}, 113:{'dx':-3., 'dy':1.}, 114:{'dx':-3., 'dy':3.}, 115: {'dx':4., 'dy':-.5}, 116:{'dx':3., 'dy':-4.}, 149:{'dx':0., 'dy':4.}, 73:{'dx':10., 'dy':1.}, 84:{'dx':-3., 'dy':0.}  }
		
	#	
	#
	# extras:
	if verbose: print "customize with map_flavor: ", map_flavor
	if map_flavor.lower()=='napa':
		plt.title('Virtual Quake ETAS map of Napa region\n\n', size=fs_title)
		mainshock = {'event_time':dtm.datetime(2014, 8, 24, 3+7, 20, 44, tzinfo=pytz.timezone('UTC')), 'lat':38.22, 'lon':-122.313, 'mag':6.0}
		plot_sections = vc_parser.napa_sections
		#		
	if map_flavor.lower() == 'emc':
		plt.title('Virtual Quake ETAS map of El Mayor-Cucapah region\n\n', size=fs_title)
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
	if verbose: print "... and now, make seismicity_map(), then do some local plotting... (or maybe not. i think we skip this since it takes forever)."
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
		f_path, f_name = os.path.split(f_out)
		if not os.path.isdir(f_path): os.makedirs(f_path)
		plt.savefig(f_out)
	#
	return my_map

def vc_etas_RI_map(map_flavor='napa', i_min=1000, i_max=4000, etas_gridsize=.1, plot_sections=[], f_out=None, plot_quake_dots=False, fault_colors = ['b'], plot_section_labels=False):
	# note, here "RI" means "Relative Intensity", NOT "recurrence interval"
	# make an etas map from a vc catalog. we probably have to pare down the catalog, since simulated catalogs are super huge.
	# sectionses: a list of a list of sections[ [secs1], [secs2], etc.] to (over)plot separately.
	#
	# we'll want to control the color cycle:
	#colors_ =  mpl.rcParams['axes.color_cycle']		# use like: this_color = colors_[j%len(colors_)]
	# ... which is a list, like: ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	if fault_colors==None: fault_colors = mpl.rcParams['axes.color_cycle']	
	if isinstance(fault_colors, str):
		# take a guess at the format:
		fault_colors = fault_colors.replace(' ', '').split(',')	#returns as list.
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
	# note: this takes a LONG time for the current ETAS model (which can be very much optimized, but is not yet). we should script up a bit to save the
	# [[xyz]] data so we can modify this map more easily.
	#
	if f_out!=None: fout_xyz = f_out[:-4]+'.xyz'
	my_map = vc_parser.seismicity_map(section_ids=sectionses[map_flavor], start_date=10000., etas_catalog={'catalog':list(catalog[i_min:i_max]), 'sweep_index':sweep_index}, p_map=0., etas_gridsize=etas_gridsize, etas_output_xyz=fout_xyz)
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
			this_color = fault_colors[j%len(fault_colors)]	# note: if we provide only one color, this will always produce that entry.
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
		f_path, f_name = os.path.split(f_out)
		if not os.path.isdir(f_path): os.makedirs(f_path)
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
def closest_emc_fault_section_figs(n_sections=5, sections=None, m0=7.0, output_dir='figs_gji/rev1/closest_sections/', fs_title=16, fs_label=16, fs_legend=12):
	'''
	# a few plots based on the nearest sections to the EMC event.
	# sections: overrides n_sections? maybe we'll truncate. anyway, you can provide sections or leave it None and fetch them.
	'''
	
	#output_dir = 'figs_gji/emc_RI_WT_figs_n_%d' % n_sections
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)

	if sections==None:
		sections = vc_parser.get_nearest_section_ids(lat_0=32.128, lon_0=-115.303, section_ids=vc_parser.emc_sections, dist_mode=0, n_sections=n_sections, verbose=False)
	print "sections: ", sections
	#
	#################
	#q=waiting_time_figs(section_ids=sections, file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], keep_figs=False, output_dir='emc_WT_figs_n_%d' % n_sections, mc_nits=100000, n_cpus=None)
	#
	# ... but we should already have these from gji_RI_probabilities(), and we don't use the individual figures anyway...
	fig_kwargs = {'fs_label':fs_label, 'fs_title':fs_title, 'fs_legend':fs_legend}
	#q1=conditional_RI_figs(section_ids=sections,file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=output_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None, **fig_kwargs)
	#
	q1=conditional_RI_figs(section_ids=[sections],file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=output_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None, **fig_kwargs)
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
	plt.xlabel('Time since last $m>%.2f$ event, $t_0$ (years)' % m0, size=fs_label)
	plt.ylabel('Expected interval $\Delta t$', size=fs_label)
	plt.title('EMC expected waiting times for fault sections:\n' + ', '.join(map(str, sections)) + '\n', size=fs_title)
	plt.savefig('%s/emc_ewt_n_%d_composite.png' % (output_dir, n_sections))
	
	return EWTf

#
def gji_RI_probabilities(section_ids=vc_parser.emc_sections, output_dir='figs_gji/rev1/pri/', production_gen_figs='/home/myoder/Dropbox/Research/VQ/VC_EMC_gji_yoder/general_figs/', fs_title=16, fs_label=16, fs_legend=12):
	#
	# ... and we should write a script to plot pre-fit data. it's probably better to just write a separate script, since matching the 
	# existing t_0 values with the parent (fit + fig generating script) is tricky (since they're floats). the better approach is to separate
	# the fitting and plotting parts... so do this in PyVQ.
	#
	if not os.path.isdir(output_dir): os.makedirs(output_dir)
	#gji_pri_figs_dir = '%s/pri' % output_dir
	gji_pri_figs_dir = output_dir
	#
	# get individual section figures:
	plt.figure(0, figsize=(9,8))
	fig_kwargs = {'fs_label':fs_title, 'fs_title':fs_title, 'fs_legend':fs_legend}
	#
	z=conditional_RI_figs(section_ids=section_ids,file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=gji_pri_figs_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None, **fig_kwargs)
	#
	# now get the aggregate figure:
	z=conditional_RI_figs(section_ids=[section_ids],file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir=gji_pri_figs_dir, mc_nits=100000, n_cpus=None, start_year=10000, end_year=None, **fig_kwargs)
	plt.title('CDF $P(\Delta t_r)$ for EMC region, aggregated', size=fs_title)
	#plt.xlabel('$m=7.0$ Recurrence interval $\\Delta t_r$ (years)', size=fs_label)
	#plt.ylabel('Probability $P(\\Delta t_r)$', size=fs_label)
	#
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

# this isn't quite right for what we'r trying to do. "pre-fit" should just be a line in the regurlar figure fitter that loads fit prams and skips fitting.
def PRI_faultwise_composite(fit_data='VC_CDF_WT_figs/VC_CDF_WT_fits_m70_section_-2.npy', fname_out='figs_gji/rev1/pri/EMC_waiting_times_faultwise_composite.png', section_ids=vc_parser.emc_sections, m0=7.0, fs_title=16, fs_label=16, fs_legend=12, sim_file=vc_parser.default_sim_file, fnum=0):
	# , t0s=[0., 180., 360., 720., 900.]
	# make the faultwise aggregated PRI distribution plot. we probably don't need to make this too slick, since i'm not sure this figure will ever be made again.
	# basically, get all fault-wise intervals on all sections; then sort them.
	# load the pre-fit data. we'd nominally have to figure out how to recalculate these if we were going to use this figure again, but actually, by itself i
	# don't think it's a great figure. it should probably be a plot of interval ratios (aka, dt/<dt> for each fault section). but even then, the interesting figures
	# are for each fault segment. this is only meant to serve as a transition from regional to faultwise analyses.
	
	# (script under development).
	#
	colors_ =  mpl.rcParams['axes.color_cycle']
	str_m0 = str(m0).replace('.','')[:2]
	if not fname_out==None:
		pth, fnm = os.path.split(fname_out)
		if not os.path.isdir(pth): os.makedirs(pth)
	#
	fitses = numpy.load(fit_data)	# for now, require this...
	#
	# now, get all intervals:
	delta_ts = []
	for sec_id in section_ids:
		X=vc_parser.get_CFF_on_section(sim_file=sim_file, section_id=sec_id)
		tm0 = [rw['event_year'] for rw in X if rw['event_magnitude']>=m0]
		#TdT += [[t, t-tm0[j]] for j,t in enumerate(tm0[1:])]
		delta_ts += [t-tm0[j] for j,t in enumerate(tm0[1:])]
	#
	delta_ts.sort()
	max_delta_ts = delta_ts[-1]*1.2
	X_fit = numpy.arange(0., max_delta_ts, max_delta_ts/10000.)
	#
	plt.figure(fnum, figsize=(12,10))
	plt.clf()
	#
	chi_t0 = fitses[0]['chi']
	beta_t0 = fitses[0]['beta']
	#
	N=float(len(delta_ts))
	for k,rw in enumerate(fitses):
		#
		clr = colors_[k%len(colors_)]
		t0=rw['t0']
		weib_chi=rw['chi']
		weib_beta = rw['beta']
		#sec_id = rw['section_id']
		#
		dts = [dt for dt in delta_ts if dt>=t0]
		n = float(len(dts))
		if n==0.: continue
		plt.plot(dts, [float(x)/n for x in xrange(len(dts))], '-', lw=2, label='data, $t_0=%.3f$' % t0, color=clr)
		#
		# use vc_parser.f_weibull(x=None, chi=1.0, beta=1.0, x0=None) for weibull fits.
		plt.plot(filter(lambda x:x>=t0,X_fit), [vc_parser.f_weibull(x=x, chi=weib_chi, beta=weib_beta, x0=t0) for x in X_fit if x>=t0], '-.', label='$t_0=%.2f$, $\\beta=%.3f$, $\\tau=%.3f$' % (t0, weib_beta, weib_chi), lw=2, color=clr)
		plt.plot(filter(lambda x:x>=t0,X_fit), [vc_parser.f_weibull(x=x, chi=chi_t0, beta=beta_t0, x0=t0) for x in X_fit if x>=t0], '-.', label='$t_0=%.2f$, $\\beta=%.3f$, $\\tau=%.3f$' % (t0, weib_beta, weib_chi), lw=2, color=clr)
		
	ax=plt.gca()
	ax.set_xlim(left=0, right=4000.)
	ax.set_ylim(0, 1.1)
	plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t$ (years)' % (m0), size=fs_label)
	plt.ylabel('Cumulative Probability $P(\\Delta t)$', size=fs_label)
	plt.title('Conditional Recurrence Probabilities for $m>%.2f$, faultwise mean' % (m0), size=fs_title)
	plt.legend(loc=0, numpoints=1, prop={'size':fs_legend})
	#
	if fname_out!=None:
		plt.savefig(fname_out)
	
	return delta_ts
	

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
def gji_forecast_fig(fignum=0, section_id=16, f_out = 'dumps/figs_gji/forecast_section_16.png', time_range=(13400., 14850.), opt_data='dumps/gji_roc_lt_500_b/roc_sec_lt_%d_nits_1000_allprams.npy',f_gt_lt=operator.lt ):
	#
	# "earthquake predictability" forecast time-series figure:
	#
	# plot_psa_metric_figure(CFF=None, section_id=None, m0=7.0, fignum=0, nyquist_factor=None, b_0=None, opt_data='dumps/gji_roc_lt_500/roc_sec_lt_%d_nits_1000_allprams.npy', lw=2., **kwargs)
	A=vc_parser.plot_psa_metric_figure(section_id=section_id, fignum=fignum, opt_data=opt_data,f_gt_lt=f_gt_lt)
	#
	# axes are added in order "rate", then "b", then "mags" clones 'rate'
	#f=plt.gcf()
	f = plt.figure(fignum)
	ax_b = f.axes[1]
	ax_rate = f.axes[0]
	#ax_mags = f.axes[2]	# but maybe this is handled internally, and no axis object is actualyl created? becasue this gives "index out of range".
	#
	ax_rate.set_xlim(time_range)		 # but note that these all share the same x values
	ax_b.set_ylim(-1.7, 1.25)
	#
	plt.draw()
	#
	# now, save the file.
	if f_out!=None:
		path_name = os.path.dirname(f_out)
		if not os.path.isdir(path_name): os.makedirs(path_name)
		print "saving to: ", f_out
		plt.savefig(f_out)
	#
	return f

#def EMC_exp_WT_Appendix(wt_dir='expected_waiting_time_figs', wt_prob_dir='VC_WT_probs_20150109', output_file = 'appendix_exp_wt.tex', section_ids=vc_parser.emc_sections):
def EMC_exp_WT_Appendix(wt_dir='\WTexpfigs', wt_prob_dir='\WTPfigs', output_file = 'vc_appendix_emc_PRI.tex', section_ids=vc_parser.emc_sections):
	'''
	# make at least the framework for an appendix of all the WT distribution figures.
	# file_format='VC_CDF_WT_m70_section_%d.png'
	#
	'''
	#
	manuscript_folder = '/home/myoder/Dropbox/Research/VQ/VC_EMC_gji_yoder'
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
	# just hard-code font sizes for now:
	# ... and we'll just define these globally (but we could override them here).
	#fs_title  = 16
	#fs_labels = 16.
	#fs_legend = 14.
	#
	if os.path.isdir(output_dir)==False:
		os.makedirs(output_dir)
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
		plt.title('Expected Waiting times: sections %s' % name_str, size=fs_title)
		plt.xlabel('Time since last $m>%.2f$ event, $t_0$' % m0, size=fs_label)
		plt.ylabel('Expected interval $\Delta t$', size=fs_label)
		#
		plt.savefig('%s/EWT_m0_%s_section_%s.png' % (output_dir, str(m0).replace('.',''), name_str))
#
def best_fit_roc_array(section_ids=vc_parser.emc_sections + ['EMC'], roc_file_format='dumps/gji_roc_lt_500/roc_sec_lt_%s_nits_1000_allprams.npy'):
	best_fit_array=[]
	cols = ('section_id', 'b', 'nyquist_factor', 'n_predicted', 'n_total', 'time_alert', 'time_total', 'H', 'F', 'score_lin')
	for sec_id in section_ids:
		rocs = numpy.load(roc_file_format % str(sec_id))
		rocs.sort(key=lambda rw:rw['H']-rw['F'])
		#
		# best prams are at the end. we might average over the last N or some other statistic. for now, just take the best fit.
		best_fits = rocs[-1]
		#
		# ... and just spell out the output variables:
		b = best_fits['b']
		alpha = best_fits['nyquist_factor']
		n_predicted = best_fits['n_predicted']
		n_total = best_fits['n_predicted'] + best_fits['n_missed']
		time_alert = best_fits['total_alert_time']
		time_total = best_fits['total_time']
		#
		if best_fits.has_key('H'):
			H = best_fits['H']
		else:
			H = n_predicted/n_total
		if best_fits.has_key('F'):
			F=best_fits['F']
		else:
			F = time_alert/time_total
		#
		best_fit_array += [[sec_id, b, alpha, n_predicted, n_total, time_alert, time_total, H, F, H-F]]
	#
	col_types = [type(x).__name__ for x in best_fit_array[0]]
	col_types[0]='S16'
	best_fit_array = numpy.core.records.fromarrays(zip(*best_fit_array), names=cols, formats=col_types)
	#
	return best_fit_array
#
def best_fit_ROC_table(section_ids=vc_parser.emc_sections + ['EMC'], output_file='dumps/gji_roc_lt_500/roc_bestfits', m0=7.0, roc_file_format='dumps/gji_roc_lt_500/roc_sec_lt_%s_nits_1000_allprams.npy', longtable=True, lbl_str='tbl:roc_prams'):
	'''
	# make a table of best fit ROC values; output in .tex format (we have this somewhere, but it's going to work a bit differently here).
	# opt for an input table/list of dicts or something, or a use filenames. use the raw roc data and, as usual, select the best-fit
	# parameter set by sorting.
	# note: m0 should be in the data file, but it's not (something to fix later -- probably add hierarchy to the data files (or 
	# we might be getting rid of them entirely if the C++ access is fast enough), something like {data:[present list or dict of dicts], meta_data:{}
	#pass
	#
	# we'll make two outputs: 1) simple CSV, 2) .tex format. do we have a csv->tex script already?
	# final copies of output files will (eventually) end up in: /home/myoder/Dropbox/Research/VQ/VC_EMC_gji_yoder/general_figs
	'''
	#
	fname_csv = '%s.csv' % (output_file)
	fname_tex = '%s.tex' % (output_file)
	#
	output_cols_csv = ['Sec. ID', 'b_0', 'alpha_n', 'N_predicted', 'N_total', 'Delta t_alert', 'Delta t_total', 'H', 'F', 'H-F']
	output_cols_tex = ['\\textbf{Sec. ID}', '$b_0$', '$\\alpha_n$', '$N_{predicted}$', '$N_{total}$', '$\\Delta t_{alert}$', '$\\Delta t_{total}$', '\\textbf{H}', '\\textbf{F}', '\\textbf{H-F}']
	#
	with open(fname_csv, 'w') as f:
		f.write('#ROC best fits for sections: %s\n' % ', '.join([str(x) for x in section_ids]))
		f.write('#!%s\n' % '\t'.join([str(x) for x in output_cols_csv]))
		f.write('#\n')
	with open(fname_tex, 'w') as f:
		f.write('%' + 'ROC best fits for sections: %s\n' % ', '.join([str(x) for x in section_ids]))
		if longtable:
			f.write('%' + '\n\\begin{longtable}{|%s|}\n' % '|'.join(['c' for x in output_cols_tex]))
			f.write('\caption{Best fit parameters and ROC scores for faults in Virtual EMC (our VQ simulation of the El Mayor-Cucapah region. The section labeled ``EMC'' shows the best fit parameters for the aggregated EMC catalog -- all fault segments together.}\n\label{%s} \\\ \n' % lbl_str)
			#
			# first, write column headers. note the \endfirsthead tag at the end.
			f.write('%s \\\ \n' % ' & '.join(['\multicolumn{1}{c}{%s}' % col for col in output_cols_tex]))
			f.write('\hline\hline\n')
			f.write('\endfirsthead\n')
			f.write('\n')
			#
			# now, write the contunuation header:
			f.write('\multicolumn{%d}{c}{{Table \\ref{%s} -- continued from previous page.}} \\\ \n' % (len(output_cols_tex), lbl_str))
			f.write('%s \\\ \n' % ' & '.join(['\multicolumn{1}{c}{%s}' % col for col in output_cols_tex]))
			f.write('\hline\n')
			f.write('\endhead\n\n')
			#
			# and a footer...
			f.write('\n \multicolumn{%d}{c}{{ -- Table \\ref{%s} continued on next page -- }} \\\ \n' % (len(output_cols_tex), lbl_str))
			f.write('\endfoot\n\n')
			#
			# and last footer...
			f.write('\hline\n')
			f.write('\endlastfoot\n\n')
			#
			end_string = '\end{longtable}\n\n'
		else:	
			f.write('%' + '\n\\begin{table}\n')
			#f.write('\caption{Best fit parameters and ROC scores for faults in Virtual EMC.}\label{tbl:roc_prams}\n\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}\hline\n')
			f.write('\caption{Best fit parameters and ROC scores for faults in Virtual EMC.}\n\label{%s}\n\\begin{tabular}{|%s|}\n' % (lbl_str, '|'.join(['c' for x in output_cols_tex])))
			f.write('%s \\\ \n' % ' & '.join(output_cols_tex))
			f.write('\hline\hline\n')
			#
			end_string = '\end{tabular}\n\end{table}\n\n'
		

	#
	for sec_id in section_ids:
		rocs = numpy.load(roc_file_format % str(sec_id))
		rocs.sort(key=lambda rw:rw['H']-rw['F'])
		#
		# best prams are at the end. we might average over the last N or some other statistic. for now, just take the best fit.
		best_fits = rocs[-1]
		#
		# ... and just spell out the output variables:
		b = best_fits['b']
		alpha = best_fits['nyquist_factor']
		n_predicted = best_fits['n_predicted']
		n_total = best_fits['n_predicted'] + best_fits['n_missed']
		time_alert = best_fits['total_alert_time']
		time_total = best_fits['total_time']
		#
		if best_fits.has_key('H'):
			H = best_fits['H']
		else:
			H = n_predicted/n_total
		if best_fits.has_key('F'):
			F=best_fits['F']
		else:
			F = time_alert/time_total
		#
		with open(fname_csv, 'a') as f:
			f.write('%s\n' % '\t'.join([str(x) for x in [sec_id, b, alpha, n_predicted, n_total, time_alert, time_total, H, F, H-F]]))
		with open(fname_tex, 'a') as f:
			f.write('%s \\\ \n' % ' & '.join(['\\textbf{%s}' % str(x) if x==sec_id else '$%.2f$' % x if isinstance(x,float) else '$%s$' % x for j,x in enumerate([sec_id, b, alpha, n_predicted, n_total, time_alert, time_total, H, F, H-F])]))
			f.write('\\hline\n')
	
	#
	with open(fname_tex, 'a') as f:
		f.write(end_string)
	
	#
	production_fig_path = '/home/myoder/Dropbox/Research/VQ/VC_EMC_gji_yoder/general_figs'
	dir_tex, fname_tex = os.path.split(fname_tex)
	dir_csv, fname_csv = os.path.split(fname_csv)
	#print "fnames...", dir_tex, " ** ", fname_tex
	os.system('cp %s %s' % (os.path.join(dir_tex, fname_tex), os.path.join(production_fig_path, fname_tex)))
	os.system('cp %s %s' % (os.path.join(dir_csv, fname_csv), os.path.join(production_fig_path, fname_csv)))
	
	
#
#def plot_best_roc(n_rank=5, save_figs=True, b_0=0., nyquist_factor=.5, output_dir='dumps/figs_gji/', font_size=16):
def plot_best_roc(n_rank=5, save_figs=True, b_0=0., nyquist_factor=.5, output_dir='dumps/figs_gji/', input_data_format='dumps/gji_roc_lt_500_b/roc_sec_lt_*_allprams*.npy'):
	# some plots of best ROC parameters:
	# ROC for best optimized forecasts:
	#   1)simple plot, only best fit, 
	#   2)more complex fit with top n=5 fits for each fault segment
	#  3) 3-d fig of score vs (b, alpha)
	#  extra figs:
	#  4) score vs b
	#  5) score vs alpha (these don't show much of a pattern except that most segments have b_optimal >0. nyquist is 
	#     pretty evenly distributed .2 < alpha <1.2 . it might be a good idea to explore alpha > 1.2.
	#
	# note: i think the right way to set the default font is something like:
	# mpl.rcParams['font.size']=...
	'''
	if font_size!=None and (isinstance(font_size, int) or isinstance(font_size, float)):
		font_size=float(font_size)
		font_size_original = mpl.rcParams['font.size']
		mpl.rcParams['font.size'] = font_size
		#print "reset mpl font.size: ", (font_size_original, font_size, mpl.rcParams['font.size'])
	'''
	#
	my_files = glob.iglob(input_data_format)
	colors_ =  mpl.rcParams['axes.color_cycle']
	#
	output_dir=output_dir.strip()
	if output_dir=='': output_dir = './'
	if output_dir!='' and output_dir[-1]!='/': output_dir+='/'
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	#
	plt.figure(0)
	plt.clf()
	plt.plot(range(2), range(2), 'r-',lw=2, zorder=1)
	#
	mean_H=0.
	mean_F=0.
	total_n_predicted=0.
	total_n = 0.
	total_t_alert = 0.
	total_t = 0.
	#
	# for reference, 'default' roc with b=0, alpha=.5 (or whatever we come up with later)
	roc_default = ROC_single_prams(section_ids=emc_sections, b_0=b_0, nyquist_factor=nyquist_factor, m0=7.0, fignum=None)
	#
	col_names = ['H', 'F','b', 'nyquist_factor']
	my_lists = {key:[] for key in col_names}
	j=0	# make a place-holder for j in case the for loop has no values to loop over, so j never gets assigned.
	for j, fl in enumerate(my_files):
		this_color = colors_[j%len(colors_)]
		roc = numpy.load(fl)
		#
		roc.sort(key=lambda x:x['H']-x['F'])
		#
		[my_lists[key].append(roc[-1][key]) for key in col_names]
		#
		f = [rw['F'] for i,rw in enumerate(roc[-n_rank-1:])]
		h = [rw['H'] for i,rw in enumerate(roc[-n_rank-1:])]
		#
		mean_H+=numpy.mean(h)
		mean_F+=numpy.mean(f)
		#
		fh=zip(f,h)
		fh.sort(key=lambda x: x[0])
		f,h = zip(*fh)
		#
		plt.plot(f,h, '-', color=this_color)
		plt.plot(f,h, 'o', color=this_color, alpha=.35)
		plt.plot([f[0], f[-1]], [h[0], h[-1]], 'o', color=this_color)
	#
	j=(0 or j)
	mean_H/=(float(j+1))
	mean_F/=(float(j+1))
	#
	plt.plot([mean_F], [mean_H], 'r*', ms=15, zorder=4, alpha=.9)
	plt.plot([mean_F], [mean_H], 'k*', ms=18, zorder=3, alpha=.9)
	#
	plt.plot([numpy.median(my_lists['F'])], [numpy.median(my_lists['H'])], 'y*', ms=15, zorder=3, label='median $(<\\bar{F}>, <\\bar{H}>)$')
	plt.plot([numpy.median(my_lists['F'])], [numpy.median(my_lists['H'])], 'b*', ms=18, zorder=2)
	#
	plt.plot([rw['F'] for rw in roc_default.itervalues()], [rw['H'] for rw in roc_default.itervalues()], 'gs', zorder=1, alpha=.8, label='default ($b_0=%.2f$, $\\alpha_n=%.2f$)' % (b_0, nyquist_factor))
	#
	plt.xlabel('False alarm rate $F$', size=fs_label)
	plt.ylabel('Hit rate $H$', size=fs_label)
	plt.legend(loc='lower right', numpoints=1, prop={'size':fs_legend})
	if save_figs: plt.savefig('%sROC_EMC_faultwise_n_%d.png' % (output_dir, n_rank))
	#
	plt.figure(1)
	plt.clf()
	plt.plot(my_lists['F'], my_lists['H'], 'bo', zorder=2, label='optimized')
	plt.plot([rw['F'] for rw in roc_default.itervalues()], [rw['H'] for rw in roc_default.itervalues()], 'gs', zorder=1, alpha=.8, label='default ($b_0=%.2f$, $\\alpha_n=%.2f$' % (b_0, nyquist_factor))
	plt.plot([numpy.mean(my_lists['F'])], [numpy.mean(my_lists['H'])], 'r*', ms=15, zorder=3, label='mean $(<F>, <H>)$')
	plt.plot([numpy.mean(my_lists['F'])], [numpy.mean(my_lists['H'])], 'k*', ms=18, zorder=2)
	
	plt.plot([numpy.median(my_lists['F'])], [numpy.median(my_lists['H'])], 'y*', ms=15, zorder=3, label='median $(<\\bar{F}>, <\\bar{H}>)$')
	plt.plot([numpy.median(my_lists['F'])], [numpy.median(my_lists['H'])], 'b*', ms=18, zorder=2)
	
	plt.plot(range(2), range(2), 'r-',lw=2, zorder=1)
	plt.xlabel('False alarm rate $F$', size=fs_label)
	plt.ylabel('Hit rate $H$', size=fs_label)
	plt.legend(loc='lower right', numpoints=1, prop={'size':fs_legend})
	if save_figs: plt.savefig('%sROC_EMC_faultwise_simple.png'% (output_dir))
	#
	f=plt.figure(2)
	plt.clf()
	ax=f.add_subplot(111, projection='3d')
	scores = [h-f for h,f in zip(my_lists['H'], my_lists['F'])]
	#ax.plot(my_lists['b'], my_lists['nyquist_factor'], scores, '.')
	for j,x in enumerate(my_lists['b']):
		this_color = colors_[j%len(colors_)]
		#ax.plot([my_lists['b'][j], my_lists['b'][j]], [my_lists['nyquist_factor'][j], my_lists['nyquist_factor'][j]], [0., scores[j]], 'o-')
		ax.plot([my_lists['b'][j], my_lists['b'][j]], [my_lists['nyquist_factor'][j], my_lists['nyquist_factor'][j]], [0., scores[j]], '-', color=this_color)
		ax.plot([my_lists['b'][j]], [my_lists['nyquist_factor'][j]], [scores[j]], 'o', color=this_color)
	
	ax.set_xlabel('Interval slope $b$', size=fs_label)
	ax.set_ylabel('nyqist factor $\\alpha_n$', size=fs_label)
	ax.set_zlabel('score $H-F$', size=fs_label)
	
	plt.figure(3)
	plt.clf()
	plt.plot(my_lists['nyquist_factor'], scores, 'o')
	plt.xlabel('nyquist factor', size=fs_label)
	plt.ylabel('score $H-F$', size=fs_label)
	
	plt.figure(4)
	plt.clf()
	plt.plot(my_lists['b'], scores, 'o')
	plt.xlabel('$b_0$', size=fs_label)
	plt.ylabel('score $H-F$', size=fs_label)
	#
	print "\n\noutput report:"
	print "optimized:"
	print "mean_H,mean_F: %f, %f (%f)" % (mean_H, mean_F, mean_H-mean_F) 
	print "medain_H, median_F: %f, %f (%f) " % (numpy.median(my_lists['H']), numpy.median(my_lists['F']), numpy.median(my_lists['H'])-numpy.median(my_lists['F']))
	#return my_lists
	print "stdev(H-F)]: %f" % (numpy.std([h-f for f,h in zip(my_lists['F'], my_lists['H'])]))
	roc_default_means = {'H':numpy.mean([rw['H'] for rw in roc_default.itervalues()]), 'F':numpy.mean([rw['F'] for rw in roc_default.itervalues()])}
	roc_default_meds = {'H':numpy.median([rw['H'] for rw in roc_default.itervalues()]), 'F':numpy.median([rw['F'] for rw in roc_default.itervalues()])}
	print "Defaults:"
	print "mean_H,mean_F: %f, %f (%f) " % (roc_default_means['H'], roc_default_means['F'], roc_default_means['H']-roc_default_means['F'])
	print "median_H,median_F: %f, %f (%f) " % (roc_default_meds['H'], roc_default_meds['F'], roc_default_meds['H']-roc_default_meds['F'])
	print "stdev[(H-F)]: %f" % (numpy.std([h-f for h,f in zip([rw['H'] for rw in roc_default.itervalues()],[rw['F'] for rw in roc_default.itervalues()])]))
	#
	# return original font size:
	#if font_size!=None and (isinstance(font_size, int) or isinstance(font_size, float)):
	#	#font_size_original = mpl.rcParams['font.size']
	#	mpl.rcParams['font.size'] = font_size_original
	#
	#return my_files
	return my_lists, roc_default

# plotting helper functions:
def roc_figure(roc_data=None, roc_random=None, CFF=None, section_id=None, fignum=0, do_clf=True, label_str=None, title_str=None, markers='.-', n_rand=1000, m0=7.0, bin_size=.1, output_dir='figs_gji/rev1/', fname=None):
	# (aka, figs 16, 17 -- from original manuscript).
	# Single section (catalog) ROC figure (includes all ROC trials, line over max ROC, error envelope)
	# construct an ROC curve for a section (or more generally, a collection of ROC optimizer data collected by some unknown parsing).
	# roc_data input is the 'full_pram' set from simple_metric_optimizer()[1] or the output file (which can be a string-filename
	# or a roc_data = numpy.load(roc_data); see for example numpy.load('dumps/roc_detail_data/roc_sec_16_allprams.npy')
	#
	# ... and later on, code up something for section_id so we can build a roc_data set on demand, if necessary.
	# ... and have a look through this, particularly the section_id parameter, to be sure it behaves properly.
	#
	# handle some input parameters:
	'''
	font_size = (font_size or 12.)
	if font_size!=None and (isinstance(font_size, int) or isinstance(font_size, float)):
		font_size=float(font_size)
		font_size_original = mpl.rcParams['font.size']
		mpl.rcParams['font.size'] = font_size
	'''
	#
	if not os.path.isdir(output_dir): os.makedirs(output_dir)
	if fname==None: fname = 'ROC_scatter_sec_%d.png' % (0 or section_id)	# in case it's None...
	fout = os.path.join(output_dir, fname)
	#
	if label_str==None:
		label_str = 'Forecast score'
	if roc_data==None:
		section_id = (16 or section_id)
		roc_data = 'dumps/gji_roc_lt_500_b/roc_sec_lt_%d_nits_1000_allprams.npy' % section_id
	#
	if isinstance(roc_data, str): roc_data = numpy.load(roc_data)
	#
	if isinstance(roc_data[0],dict):
		# list of dicts. let's turn it into a rec-array.
		col_names = [key for key,val in roc_data[0].iteritems() if not isinstance(val,str)]		#roc_data[0].keys()
		#lst_data = [[rw[key] for key in col_names] for rw in roc_data]
		roc_data = numpy.core.records.fromarrays(zip(*[[rw[key] for key in col_names] for rw in roc_data]), names=col_names, formats=[type(roc_data[0][key]).__name__ for key in col_names])
	#
	# and a random forecast:
	# (and we need to save these so that we can do quick reformats).
	if roc_random==None: roc_random = vc_parser.get_random_forecast_set(CFF=CFF, section_id=section_id, m0=m0, P_min=.05, P_max=1.0, set_name=None, nits=n_rand, format='recarray', do_roc_plot=False, fignum=None)
	#
	# bin up the data for aggregate curves:
	num_points = int(1.0/bin_size)
	roc_curve_data = [[] for x in xrange(num_points)]
	roc_curve_rand = [[] for x in xrange(num_points)]
	#
	for i,rw in enumerate(roc_data):
		bin_num = min(num_points-1, int(rw['F']/bin_size))
		roc_curve_data[bin_num]+=[i]
	#
	for i,rw in enumerate(roc_random):
		bin_num = min(num_points-1, int(rw['F']/bin_size))
		#print bin_num
		roc_curve_rand[bin_num]+=[i]
	# so now, we have groups of row indices for each bin.
	#
	roc_X = numpy.arange(0., 1.+bin_size, bin_size)
	#roc_Y = [None] + [len(x) for j,x in enumerate(roc_curve_data)]
	roc_data_max = [0.] + [max([roc_data['H'][i] for i in rw]) if len(rw)>0 else None for j,rw in enumerate(roc_curve_data)]
	#roc_min_random = [0.] + [min([roc_random['H'][i] for i in rw]) if len(rw)>0 else None for j,rw in enumerate(roc_curve_rand)]
	#roc_max_random = [0.] + [max([roc_random['H'][i] for i in rw]) if len(rw)>0 else None for j,rw in enumerate(roc_curve_rand)]
	
	roc_rand_std  = [0.] + [numpy.std([roc_random['H'][i] for i in rw]) if len(rw)>0 else None for j,rw in enumerate(roc_curve_rand)]
	#roc_rand_mean = [0.] + [numpy.mean([roc_random['H'][i] for i in rw]) if len(rw)>0 else None for j,rw in enumerate(roc_curve_rand)]
	roc_rand_mean = roc_X
	roc_min_random = [x-v for x,v in zip(roc_rand_mean, roc_rand_std)]
	roc_max_random = [x+v for x,v in zip(roc_rand_mean, roc_rand_std)]
	#
	max_lin = sorted([[h,f, h-f] for h,f in zip(roc_data['H'], roc_data['F'])], key=lambda x:x[2])[-1]
	#max_geom = sorted([[h,f, h/f] for h,f in zip(roc_data['H'], roc_data['F'])], key=lambda x:x[2])[-1]
	#
	plt.ion()
	plt.figure(fignum)
	plt.clf()
	plt.plot(range(2), range(2), 'r-', lw=2, zorder=4)
	plt.plot(roc_data['F'], roc_data['H'], 'b.', alpha=.8, zorder=3, label=label_str)
	#print roc_X, roc_data_max
	plt.plot(roc_X, [x+.02 if x!=None else None for x in roc_data_max], 'b-', lw=2.0, zorder=4)
	#
	plt.plot([max_lin[1]], [max_lin[0]], 'r*', ms=18, label='$max(H-F): H(%.2f)=%.2f$' % (max_lin[1], max_lin[0]), zorder=5)
	#plt.plot([max_geom[1]], [max_geom[0]], 'm*', ms=18, label='max(h/f)', zorder=5)
	#
	plt.plot(roc_random['F'], roc_random['H'], 'g.', alpha=.5, zorder=2, label='Random forecast')
	plt.plot(roc_X, roc_min_random, 'g-', alpha=.6, zorder=2,lw=1.5)
	plt.plot(roc_X, roc_max_random, 'g-', alpha=.6, zorder=2, lw=1.5)
	plt.fill_between(roc_X, roc_min_random, roc_max_random, color='g', alpha=.2)
	#
	plt.xlabel('False Alarm Rate $F$', size=fs_label)
	plt.ylabel('Hit Rate $H$', size=fs_label)
	plt.legend(loc=0, numpoints=1, prop={'size':fs_legend})
	if title_str!=None: plt.title(title_str, size=fs_title)
	#
	plt.savefig(fout)
	#
	# and pickle the roc_data (we'll figure out how to replot them later)... and for now, be sloppy with the name:
	with open('%s_data.pkl' % fout[:-4], 'wb') as f_pkl:
		cPickle.dump(roc_data, f_pkl)
	with open('%s_random.pkl' % fout[:-4], 'wb') as f_pkl:
		cPickle.dump(roc_random, f_pkl)	
	
	#
	#if font_size!=None and (isinstance(font_size, int) or isinstance(font_size, float)):
	#	#font_size_original = mpl.rcParams['font.size']
	#	mpl.rcParams['font.size'] = font_size_original
	#
	return roc_data, roc_random
		

#####
# data and preliminary figures:
#
def create_ROC_figs_GT_data(section_ids = vc_parser.emc_sections, nits=2500, fnum=0, num_roc_points=100, output_dir = 'dumps/gji_roc_gt_detail', output_descriptions=[], fig_title_strs=[], m0=7.0):
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
	if isinstance(output_descriptions, str): output_descriptions = [output_descriptions]
	if not hasattr(output_descriptions, '__len__'): output_descriptions=[]
	#
	if isinstance(fig_title_strs, str): fig_title_strs = [fig_title_strs]
	if not hasattr(fig_title_strs, '__len__'): fig_title_strs=[]
	#
	#
	R = random.Random()
	#
	#for sec_id in section_ids:
	for j, sec_id in enumerate(section_ids):
		if isinstance(sec_id, float): sec_id=int(sec_id)
		if isinstance(sec_id, int): sec_id=[sec_id]
		#
		if len(output_descriptions)>=(j+1) and output_descriptions[j]!=None:
			sec_str = output_descriptions[j]
		else:
			sec_str = '_'.join([str(x) for x in sec_id])
		#
		if len(fig_title_strs)>=(j+1) and fig_title_strs[j] != None:
			fig_title_str = fig_title_strs[j]
		else:
			fig_title_str = "Section %s" % sec_str		
		#
		sec_str = '_'.join([str(x) for x in sec_id])
		#
		f_output_name = '%s/roc_sec_%s_nits_%d.npy' % (output_dir, sec_str, nits)
		f_output_name_rand = '%s/roc_sec_rand_%s_nits_%d.npy' % (output_dir, sec_str, nits)
		#
		CFF = vc_parser.combine_section_CFFs(sec_id)
		opt_datas, raw_datas = simple_metric_optimizer(CFF=CFF, m0=m0, b_min=-.25, b_max=.25, d_b=.01, nyquist_min=.2, nyquist_max=1.5, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name, f_gt_lt=operator.gt, f_score=operator.sub, section_id=sec_id, opt_func=vc_parser.psa_forecast_1)
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
		plotted = plot_section_ROC_curve(roc_data=raw_datas, section_id=None, fignum=fnum, num_points=num_roc_points, label_str='PSQ forecast', markers='.-')
		
		#plotted_rand = plot_section_ROC_curve(roc_data=raw_datas_rand, section_id=None, fignum=fnum, do_clf=False, num_points=num_roc_points, label_str='(random)', markers='.')
		#
		plt.figure(fnum)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		#plt.title('Optimal ROC for Section %s' % sec_str)
		plt.title('Optimal ROC (PSA) for %s' % fig_title_str)
		plt.savefig('%s/roc_opt_sec_%s_nits_%d.png' % (output_dir, sec_str, nits))
		#
		plt.figure(fnum+1)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		#plt.title('Raw ROC for Section %s' % sec_str)
		plt.title('Raw ROC (PSA) for %s' % fig_title_str)
		plt.savefig('%s/roc_raw_sec_%s_nits_%d.png' % (output_dir, sec_str, nits))
#
def create_ROC_aggregate(section_ids=[vc_parser.emc_sections], nits=1000, fnum=0, num_roc_points=100, output_dir_lt='dumps/gji_roc_lt_500', output_dir_gt='dumps/gji_roc_gt_detail', m0=7.0):
	z_lt=create_ROC_figs_LT_data(section_ids=section_ids, nits=nits, fnum=fnum, num_roc_points=num_roc_points, output_dir=output_dir_lt, m0=m0)
	#
	# ... and why not... also do one for the GT metric:
	#z_gt=create_ROC_figs_GT_data(section_ids=section_ids, nits=nits, fnum=fnum, num_roc_points=num_roc_points, output_dir=output_dir_gt, m0=m0)
	#
	return None
#
def create_ROC_figs_LT_data(section_ids = vc_parser.emc_sections, nits=2500, fnum=0, num_roc_points=100, output_dir = 'dumps/gji_roc_lt_500', output_descriptions=[], fig_title_strs=[], m0=7.0, delta_b1=0., f_gt_lt=operator.lt):
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
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	if isinstance(output_descriptions, str): output_descriptions = [output_descriptions]
	if not hasattr(output_descriptions, '__len__'): output_descriptions=[]
	#
	if isinstance(fig_title_strs, str): fig_title_strs = [fig_title_strs]
	if not hasattr(fig_title_strs, '__len__'): fig_title_strs=[]
	#
	# a string to indicate whether we're doing LT or GT operations:
	str_lt_gt = {operator.lt:'lt', operator.gt:'gt'}.get(f_gt_lt, 'some_operator')

	#
	R=random.Random()
	for j, sec_id in enumerate(section_ids):
		if isinstance(sec_id, float): sec_id=int(sec_id)
		if isinstance(sec_id, int): sec_id=[sec_id]
		
		#
		if len(output_descriptions)>=(j+1) and output_descriptions[j]!=None:
			sec_str = output_descriptions[j]
		else:
			sec_str = '_'.join([str(x) for x in sec_id])
		#
		if len(fig_title_strs)>=(j+1) and fig_title_strs[j] != None:
			fig_title_str = fig_title_strs[j]
		else:
			fig_title_str = "Section %s" % sec_str
		#
		#print "sec_str: ", sec_str
		#print "fig_title: ", fig_title_str
		#continue
		#
		f_output_name = '%s/roc_sec_%s_%s_nits_%d.npy' % (output_dir, str_lt_gt, sec_str, nits)
		f_output_name_rand = '%s/roc_sec_%s_rand_%s_nits_%d.npy' % (output_dir, str_lt_gt, sec_str, nits)
		#
		CFF = vc_parser.combine_section_CFFs(sec_id)
		#
		opt_datas, raw_datas = simple_metric_optimizer(CFF=CFF, section_id=sec_id, m0=m0, b_min=-.25, b_max=.25, d_b=.01, delta_b1=delta_b1, nyquist_min=.2, nyquist_max=1.2, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name, f_gt_lt=f_gt_lt, f_score=operator.sub, opt_func=vc_parser.psa_forecast_1b)
		#
		roc_random = vc_parser.get_random_forecast_set(CFF=CFF, section_id=sec_id, m0=m0, P_min=0., P_max=1.0, set_name=None, nits=nits, format='recarray', do_roc_plot=False, fignum=None)
		#
		# now, randomize this catalog and plot to show contrast:
		#CFF_rand = vc_parser.random_like_CFF(CFF_in=None, section_id=sec_id)
		#
		#opt_datas_rand, raw_datas_rand = simple_metric_optimizer(CFF=CFF_rand, m0=7.0, b_min=-.25, b_max=.25, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=nits, keep_set=True, set_name=None, dump_file=f_output_name_rand, f_gt_lt=operator.lt, f_score=operator.div, section_id=None)
		
		
		#return raw_datas
		# for now, put this here (get it done the first time). later, we'll move this off-line and use pre-calculated data, etc.
		plotted =      plot_section_ROC_curve(roc_data=raw_datas, section_id=None, fignum=fnum, num_points=num_roc_points, markers='.-', label_str='PSA forecast')
		#plotted_rand = plot_section_ROC_curve(roc_data=raw_datas_rand, section_id=None, fignum=fnum, do_clf=False, num_points=num_roc_points, label_str='(random)',markers='.')
		plt.figure(fnum)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		#plt.title('Optimal ROC (PSA) for Section %s' % sec_str)
		plt.title('Optimal ROC (PSA) for %s' % fig_title_str)
		plt.legend(loc=0, numpoints=1)
		plt.savefig('%s/roc_psa_%s_opt_sec_%s_nits_%d.png' % (output_dir, str_lt_gt, sec_str, nits))
		#
		plt.figure(fnum+1)
		plt.plot(roc_random['F'], roc_random['H'], '.', alpha=.6, zorder=1, label='Random forecast')
		#plt.title('Raw ROC (PSA) for Section %s' % sec_str)
		plt.title('Raw ROC (PSA) for %s' % fig_title_str)
		plt.legend(loc=0, numpoints=1)
		plt.savefig('%s/roc_psa_%s_raw_sec_%s_nits_%d.png' % (output_dir, str_lt_gt, sec_str, nits))	
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

	
