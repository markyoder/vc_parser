'''
# script to create figures and other things for VC El Mayor Cucapah and Napa papers and posters. 
# we want this to feel just like being inside vc_parser and/or vc_geodetic. we'll want to pull some 
# functions, etc. directly from there maybe, but mostly we'll be calling from those scripts in very specific, not
# portable ways. in the end, we want to push a button and have all the datasets, figs, etc. come out (mabye with a bit of layering).
'''
from vc_parser import *
from vc_geodetic import *

def vc_etas_RI_map(sections=emc_sections):
	# make an etas map from a vc catalog. we probably have to pare down the catalog, since simulated catalogs are super huge.
	#
	# first, get a catalog (we can do this inline wiht the map script, but this provides some added flexibility):
	catalog = vc_ETAS_catalog(section_ids=sections, sim_file=default_sim_file, start_year=10000., end_year=None, n_cpus=None, fignum=0, map_size=[10,8], etas_mc=3.0, etas_contour_intervals=24, p_map=0.0, sweeps_index=None)
	#
	# now, make a map using part of this catalog:
	
