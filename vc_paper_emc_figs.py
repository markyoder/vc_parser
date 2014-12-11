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

def vc_etas_RI_map(map_flavor='napa', i_min=1000, i_max=4000, etas_gridsize=.1):
	# make an etas map from a vc catalog. we probably have to pare down the catalog, since simulated catalogs are super huge.
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
		#
	#
	ms_lon, ms_lat = my_map.cm(mainshock['lon'], mainshock['lat'])
	quakes_x, quakes_y = my_map.cm(catalog['lon'][i_min:i_max], catalog['lat'][i_min:i_max])
	ln_ms1 = plt.plot(ms_lon, ms_lat, 'k*', ms=18, zorder=6, alpha=.9)
	ln_ms2 = plt.plot(ms_lon, ms_lat, 'r*', ms=15, zorder=7, alpha=.9)
	#
	ln_quakes = plt.plot(quakes_x, quakes_y, 'm.', zorder=5, alpha=.6)

def vc_surface_deformation(map_flavor = 'napa'):
	#	
	#blockwise_obj = blockwise_slip(sim_file=default_sim_file, faults=None, sections=None, pipe=None, f_pickle_out='dumps/blockwise_slip_0.pkl', plot_factor=1.0)
	# this takes a while, so a pre-calculated version might be preferable:
	#blockwise_obj = blockwise_slip(sim_file=default_sim_file, sections=napa_sections, f_pickle_out=None, plot_factor=1.0)
	#
	# pre-calculated:
	blockwise_obj='dumps/blockwise_slip_Nov11.pkl'
	disp_field = vc_geodetic.slip_field(blockwise_obj='dumps/blockwise_slip_Nov11.pkl', section_ids=sectionses[map_flavor], dx=2.0, dy=2.0, i_start=0, i_stop=-1, f_out='temp/slip_field.pkl')	# dx, dy are factors of fault-length, so \Delta_x = block_x_len * dx
	#
	zz=plot_disp_vector_field(okada_disps=disp_field, fignum=0, plot_factor=50., section_ids=sectionses[map_flavor])
	#
	# and now decorate it:
		
	
		
