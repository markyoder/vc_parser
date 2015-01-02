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

def vc_etas_RI_map(map_flavor='napa', i_min=1000, i_max=4000, etas_gridsize=.1, plot_sections=[], f_out=None, plot_quake_dots=False):
	# make an etas map from a vc catalog. we probably have to pare down the catalog, since simulated catalogs are super huge.
	# sectionses: a list of a list of sections[ [secs1], [secs2], etc.] to (over)plot separately.
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
		# get fault traces:
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
		for fault_id, fault in ft.iteritems():
			this_color = 'r'
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
		
