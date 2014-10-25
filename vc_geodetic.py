'''''''''''''''''''''''''''''''''''''''''''''
# vc_geodidic.py
# code property of University of California Davis, Dept. of Physics
# Rundle research group
# additional affiliations:
#  NASA Jet Propulsion Laboratory
#  NASA QuakeSim and E-DECIDER missions
#
# author: Mark R. Yoder, PhD
# email:  mark.yoder@gmail.com
#         mryoder@ucdavis.edu
#
# summary:
# produce various forms of (simulated) geodetic data, figures, etc. from the Virtual California (VC)
# earthquake simulator (or more generally, the Virtual Quake (VQ) earthquake simulator framework). The 
# basic strategy will be to infer surface displacement, analogous to GPS station displacements, from
# VC event-slip data.
#
'''''''''''''''''''''''''''''''''''''''''''''

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.font_manager as mfont
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import mpl_toolkits.basemap as bmp
#from mpl_toolkits.basemap import Basemap

import matplotlib as mpl
from matplotlib import cm
import itertools
plt.ion()
#
import math
import h5py
import numpy
import scipy
import scipy.optimize as spo
import operator
import glob
import random
#
#import cStringIO
import sys
import json
import cPickle
import time
import os
import datetime as dtm
import pytz
#
import vc_parser
#
#import imp
#import inspect		# use this module to determine a function's parameters and similar tasks.
import multiprocessing as mpp
#
#import ANSStools
#import BASScast
#
napa_region_section_filter = {'filter':set([45, 50, 172, 171, 170, 169, 44, 168, 167, 139, 40, 142, 41, 46])}

emc_section_filter = {'filter': (16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 69, 70, 73, 83, 84, 92, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 123, 124, 125, 126, 149)}
allcal_full_mks = '../ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'
default_sim_file=allcal_full_mks
#
emc_sections = emc_section_filter['filter']
napa_sections = napa_region_section_filter['filter']
colors_ =  mpl.rcParams['axes.color_cycle']
#
def blocks_test(n_cycles=1):
	#
	print "doing vanilla version:"
	total_time=0.0
	for i in xrange(n_cycles):
		t0=time.time()
		B=get_blocks_dict()
		t1=time.time()-t0
		print "delta_t: %f" % t1
		total_time+=t1
	print "mean time: %f" % (total_time/float(n_cycles))
	#
	total_time=0.
	print "doing Pool.apply_async() version"
	for i in xrange(n_cycles):
		t0=time.time()
		B=get_blocks_dict_pool()
		t1=time.time()-t0
		print "delta_t: %f" % t1
		total_time+=t1
	print "mean time: %f" % (total_time/float(n_cycles))
	#
	#
	total_time=0.
	print "doing Process() version"
	for i in xrange(n_cycles):
		t0=time.time()
		B=get_blocks_dict_process()
		t1=time.time()-t0
		print "delta_t: %f" % t1
		total_time+=t1
	print "mean time: %f" % (total_time/float(n_cycles))
	#
	print "finished."
#
#
def get_blocks_dict_process(sim_file=default_sim_file, faults=None, sections=None, n_cpus=None):
	# 
	if n_cpus==None: n_cpus = mpp.cpu_count()
	#
	
	# pool with 1 cpu is slower by about a factor of 1.5
	pipe_func, pipe_proc = mpp.Pipe()
	process = mpp.Process(target = get_blocks_dict, kwargs={'sim_file':sim_file, 'faults':faults, 'sections':sections, 'pipe':pipe_func})
	
	process.start()
	#process.join()
	#
	return pipe_proc.recv()
#
def get_blocks_dict_pool(sim_file=default_sim_file, faults=None, sections=None, n_cpus=None):
	# 
	if n_cpus==None: n_cpus = mpp.cpu_count()
	#
	
	# pool with 1 cpu is slower by about a factor of 1.5
	pool = mpp.Pool(n_cpus)
	X = []
	j=0
	
	X += [pool.apply_async(func=get_blocks_dict, kwds={'sim_file':sim_file, 'faults':faults, 'sections':sections})]
	pool.close()
	#
	#return pool.get()
	if len(X)==1: return X[0].get()
	if len(X)>1:
		# reduce() syndax?
		r_dict = X[0].get()
		for x in X[1:]:
			r_dict.update(x.get())
		#
		return r_dict
		
def get_blocks_dict(sim_file=default_sim_file, faults=None, sections=None, pipe=None):
	# return a dictionary of block dict-like objects (they may be recArrays.
	# this step can likely be skipped, since block_info_table is typically written
	# in block_id order, and so already indexed. this, however, is not necessarilly
	# enforced.
	#
	
	#block_dict = {}		#--> {block_id:{values...}}
	#
	with h5py.File(sim_file, 'r') as vc_data:
		block_info = vc_data['block_info_table']
		'''
		for rw in block_info:
			#block_dict[rw['block_id']] = rw		# if we want this to actually be a dict, it might be necessary to 
												# loop over vc_data['info_table'].dtype.names
			block_dict[rw['block_id']] = {name:rw[name] for name in block_info.dtype.names}
		'''
		block_dict = {rw['block_id']:{name:rw[name] for name in block_info.dtype.names} for rw in block_info if (sections==None or rw['section_id'] in sections) or (faults==None or rw['fault_id'] in faults)}
		#
	#
	if pipe!=None:
		pipe.send(block_dict)
	else:
		return block_dict
#
def block_slip_for_mpp(rws, block_info, events_data):
	return_positions = {key:[key, 0., 0., 0., 0.] for key in block_info.keys()}
	for i_rw, rw in enumerate(rws):
		#
		# we'll need to get the direction of motion... which nominally will need to be encoded into block_info.
		# for now, assume everything moves in the same direction, phi=0., theta=pi/2. (see above)
		block_id=rw['block_id']
		event_number = rw['event_number']
		#event_time = vc_data['event_table'][event_number]['event_year']
		event_time = events_data[event_number]
		#
		slip = rw['slip']
		#
		theta = block_info[block_id]['slip_theta']
		phi   = block_info[block_id]['slip_phi']
		#
		#x0, y0, z0 = block_info[key]['positions'][-1][1:4]
		x0, y0, z0 = return_positions[key][-1][2:5]
		#
		return_positions[block_id] += [[block_id, event_time, slip*math.cos(theta)*math.cos(phi) + x0, slip*math.cos(theta)*math.sin(phi) + y0, slip*math.sin(theta) + z0]]
		#
		if i_rw%10**5==0:
			print 'rw: %d (dt=%f)' % (i_rw, time.time()-t0)
	return return_positions

#
def blockwise_slip_mpp(sim_file=default_sim_file, faults=None, sections=None, pipe=None, f_pickle_out='dumps/blocwise_slip.pkl', n_cpus=None):
	# this might be working, but it pulls a lot of memory (>4GB), so needs to be tested on another machine.
	if n_cpus == None: n_cpus = mpp.cpu_count()
	#
	t0=time.time()
	print "getting blocks_dict...", t0
	block_info = get_blocks_dict(sim_file=sim_file, faults=faults, sections=sections)
	print "block_dict fetched: %f/%f" % (time.time(), time.time()-t0)
	#
	t0=time.time()
	print "block info fetched. assign mean values to blocks :: %f" % t0
	#
	# add mean position to block_info:
	#for key in block_info.keys():
	for key, rw in block_info.iteritems():
		# items() should be faster than the keys() approach...
		#rw=block_info[key]
		mean_x = numpy.mean([rw['m_x_pt%d' % j] for j in [1,2,3,4]])
		mean_y = numpy.mean([rw['m_y_pt%d' % j] for j in [1,2,3,4]])
		mean_z = numpy.mean([rw['m_z_pt%d' % j] for j in [1,2,3,4]])
		#
		block_info[key].update({'mean_x':mean_x, 'mean_y':mean_y, 'mean_z':mean_z})
		block_info[key]['slip_phi'] = 0.
		block_info[key]['slip_theta'] = math.pi/2.
		#
		# and set a field for a slip-sequence.
		block_info[key]['positions'] = [[0.0, mean_x, mean_y, mean_z]]
	#
	print 'blocks finished. index events', time.time()-t0
	
	t0=time.time()
	#
	# so now, we can spin through the event_sweep and generate positional vectors. how do we get directon?
	# use neighboring elements? this will be a bit tricky since elements may be stacked on top of one another.
	# note, however, fault/section/block do contain rake/strike, etc. info...
	#
	# now, spin through event_sweep_table and add to the slip-sequence.
	# note: we'll use the event_sweep_table and the event_table. nominally (see discussion of blocks_dict)
	# we should use a dict (or similar structure) to explicitly index event_table. however, event_table SHOULD
	# be pseudo indexed simply by order -- event_number is sequential and unique; the table is ordered by
	# event_number.
	#
	t1=time.time()
	print "block info ready. now sweep events :: %f/%f" % (t1, t1-t0)
	t0=time.time()
	#
	with h5py.File(sim_file, 'r') as vc_data:
		t0=time.time()
		events_data = {key:val for key,val in zip(*[vc_data['event_table']['event_number'], vc_data['event_table']['event_year']])}
		#events_data = vc_data['event_table']['event_year']
		# vc_data['event_table'][event_number]['event_year']
		print "events are indexed: %f" % (time.time()-t0)
		t0-time.time()
		#
		N=len(vc_data['event_sweep_table'])
		dN = N/n_cpus
		pool = mpp.Pool(n_cpus)
		i=0
		pool_results = []
		while dN*(i+1)<N:
			pool_results += [pool.apply_async(block_slip_for_mpp, args=(vc_data['event_sweep_table'][i*dN:(i+1)*dN], block_info, events_data))]
		pool.close()
		#
		print "now, compile from mpp results: ", time.time()-t0
		t0=time.time()
		for res in pool_results:
			R = res.get()
			for key in R.keys():
				if block_info[key].has_key('positions'): 
					block_info[key]['positions'].update(R[key])
				else:
					block_info[key]['positions'] = R[key]
		#
		'''
		for i_rw, rw in enumerate(vc_data['event_sweep_table']):
			#
			# we'll need to get the direction of motion... which nominally will need to be encoded into block_info.
			# for now, assume everything moves in the same direction, phi=0., theta=pi/2. (see above)
			block_id=rw['block_id']
			event_number = rw['event_number']
			#event_time = vc_data['event_table'][event_number]['event_year']
			event_time = events_data[event_number]
			
			slip = rw['slip']
			#
			theta = block_info[block_id]['slip_theta']
			phi   = block_info[block_id]['slip_phi']
			#
			#x0=block_info[key]['positions'][-1][1]
			#y0=block_info[key]['positions'][-1][2]
			#z0=block_info[key]['positions'][-1][3]
			x0, y0, z0 = block_info[key]['positions'][-1][2:5]
			#
			block_info[block_id]['positions'] += [[block_id, event_time, slip*math.cos(theta)*math.cos(phi) + x0, slip*math.cos(theta)*math.sin(phi) + y0, slip*math.sin(theta) + z0]]
			#
			if i_rw%10**5==0:
				print 'rw: %d (dt=%f)' % (i_rw, time.time()-t0)
		'''
		#
	print "finished spinning all the sweep events: ", time.time()-t0
	#
	t0=time.time()
	print "converting positions to recarrays."
	#
	# and for convenience, convert ['positions'] to a recarray:
	for key in block_info.iterkeys():
		# outputs = numpy.core.records.fromarrays(zip(*outputs), names=output_names, formats = [type(x).__name__ for x in outputs[0]])
		block_info[key]['positions'] = numpy.core.records.fromarrays(zip(*block_info[key]['positions']), names=['block_id', 'event_year', 'x', 'y', 'z'], formats = [type(x).__name__ for x in block_info[key]['positions'][0]] )
	#
	print "finished: %f" % (time.time()-t0)
	#
	print "try to pickle:"
	try:
		with open(f_pickle_out, 'w') as f_pickle:
			cPickle.dump(block_info, f_pickle)
		#
	except:
		print "failed to pickle..."
	
	return block_info
#
def blockwise_slip(sim_file=default_sim_file, faults=None, sections=None, pipe=None, f_pickle_out='dumps/blockwise_slip_0.pkl', plot_factor=1.0):
	t0=time.time()
	print "getting blocks_dict...", t0
	block_info = get_blocks_dict(sim_file=sim_file, faults=faults, sections=sections)
	#with h5py.File(sim_file, 'r') as vc_data:
	#	block_info=vc_data['block_info'].copy()
	
	print "block_dict fetched: %f/%f" % (time.time(), time.time()-t0)
	#
	t0=time.time()
	print "block info fetched. assign mean values to blocks :: %f" % t0
	#
	# add mean position to block_info:
	#for key in block_info.keys():
	#for key, rw in block_info.items():
	#for key in block_info.iterkeys():
	for key, rw in block_info.iteritems():
		#rw = block_info[key]
		# items() should be faster than the keys() approach...
		#rw=block_info[key]
		#
		#mean_x = numpy.mean([rw['m_x_pt%d' % j] for j in [1,2,3,4]])
		#mean_y = numpy.mean([rw['m_y_pt%d' % j] for j in [1,2,3,4]])
		#mean_z = numpy.mean([rw['m_z_pt%d' % j] for j in [1,2,3,4]])
		# in a single list comprehension:
		mean_x, mean_y, mean_z = [numpy.mean([rw['m_%s_pt%d' % (xyz, j)] for j in [1,2,3,4]]) for xyz in ('x', 'y', 'z')]
		#
		block_info[key].update({'mean_x':mean_x, 'mean_y':mean_y, 'mean_z':mean_z})
		#
		# slip angles: we'll need to get actual geodetic information about the faults from here. the rake/dip angles appear
		# to speak in terms of stress accumulation and are transformed "along the fault trace", or something like that. aka,
		# a vertical strike/slip fault has phi=0. or pi, theta = pi/2 (where, of course, the angle=0 position can be defined
		# arbitratily.
		#
		# for now, use holding values. going forward, we'll need to pull from elsewhere or do some fitting.
		# note: there are "m_trace_flag_pt{x} == 1" valuse for _pt1 and _pt4, but none for _pt2, _pt3, specifically
		# N=13482 N_1,4 = 2815 (about 1/4.7)... and i think the best thing to do is to pick this up in a second loop through the dict.
		# note: the "depth_id" field indicates the depth of stacked elements. aka, at x,y, there area  bunch of z-varying elements
		# with das_id={0,1,2,3..}, all with the same das_id. each (x,y) step along the segment is a unique das_id.
		#
		# that said, each block consists of 4 (x,y,z) vertices. we can determine the strike and dip from this plane.
		# it looks like a typical block is like {1:TL, 2:BL, 3:BR, 4:TR}. we can also solve this more generally via planar geometry,
		# or we could use a little trick of mean \delta{x,y,z} to solve for directions which requires only that pt. 1 be
		# considered the "origin", but this can create weird artifacts when elements are tilted. for now, just assume geometry.
		#
		# strike is mean surface angle (spherical phi) between (1,4) and (2,3).
		# dip is mean vertical angle (spherical theta) between (1,2) and (4,3) (or maybe opposite of that?)
		# so, calc dx, dy, dz for both pairs and then spherical coords accordingly.
		dx = rw['m_x_pt4'] - rw['m_x_pt1'] + rw['m_x_pt3'] - rw['m_x_pt2']
		dy = rw['m_y_pt4'] - rw['m_y_pt1'] + rw['m_y_pt3'] - rw['m_y_pt2']
		#dz = rw['m_z_pt4'] - rw['m_z_pt1'] + rw['m_z_pt3'] - rw['m_z_pt2']
		fault_phi = math.atan(dy/dx)
		#
		dx = -rw['m_x_pt4'] - rw['m_x_pt1'] + rw['m_x_pt3'] + rw['m_x_pt2']
		dy = -rw['m_y_pt4'] - rw['m_y_pt1'] + rw['m_y_pt3'] + rw['m_y_pt2']
		dz = -rw['m_z_pt4'] - rw['m_z_pt1'] + rw['m_z_pt3'] + rw['m_z_pt2']
		fault_theta = math.atan(math.sqrt(dx*dx + dy*dy)/dz)	# this is (should be) the fault dip... right?
		#														# but fault dip is a listed parameter (with rake).
		#
		#fault_phi   = - .75*math.pi
		#fault_theta = - .5*math.pi
		#
		# ... actually, this is all wrong (the right approach to calculations, but not the right angles).
		# now, these angles may need to be adjusted for direction and origin.
		block_info[key]['slip_theta'] = block_info[key]['dip_rad']  + fault_theta		# tpyically pi/2 for vertical strike/slip faults.
		block_info[key]['slip_phi']   = block_info[key]['rake_rad'] + fault_phi			# typically pi for strike/slip along the fault.
		# now we need unit vectors for slip on this block (aka, combining strike, dip, rake):
		#
		# and set a field for a slip-sequence.
		block_info[key]['positions'] = [[0.0, mean_x, mean_y, mean_z, 0.]]	# [[time, x,y,z, slip]]
	#
	print 'blocks finished. index events', time.time()-t0
	t0=time.time()
	#
	# so now, we can spin through the event_sweep and generate positional vectors. how do we get directon?
	# use neighboring elements? this will be a bit tricky since elements may be stacked on top of one another.
	# note, however, fault/section/block do contain rake/strike, etc. info...
	#
	# now, spin through event_sweep_table and add to the slip-sequence.
	# note: we'll use the event_sweep_table and the event_table. nominally (see discussion of blocks_dict)
	# we should use a dict (or similar structure) to explicitly index event_table. however, event_table SHOULD
	# be pseudo indexed simply by order -- event_number is sequential and unique; the table is ordered by
	# event_number.
	#
	t1=time.time()
	print "block info ready. now sweep events :: %f/%f" % (t1, t1-t0)
	t0=time.time()
	#
	with h5py.File(sim_file, 'r') as vc_data:
		t0=time.time()
		events_data = {key:val for key,val in zip(*[vc_data['event_table']['event_number'], vc_data['event_table']['event_year']])}
		#events_data = vc_data['event_table']['event_year']
		# vc_data['event_table'][event_number]['event_year']
		print "events are indexed: %f" % (time.time()-t0)
		t0-time.time()
		#
		for i_rw, rw in enumerate(vc_data['event_sweep_table']):
			#
			# we'll need to get the direction of motion... which nominally will need to be encoded into block_info.
			# for now, assume everything moves in the same direction, phi=0., theta=pi/2. (see above)
			block_id=rw['block_id']
			event_number = rw['event_number']
			#event_time = vc_data['event_table'][event_number]['event_year']
			event_time = events_data[event_number]
			
			slip = rw['slip']*plot_factor
			#
			theta = block_info[block_id]['slip_theta']
			phi   = block_info[block_id]['slip_phi']
			#
			#x0=block_info[key]['positions'][-1][1]
			#y0=block_info[key]['positions'][-1][2]
			#z0=block_info[key]['positions'][-1][3]
			x0, y0, z0 = block_info[block_id]['positions'][-1][1:4]
			#
			#block_info[block_id]['positions'] += [[block_id, event_time, slip*math.cos(theta)*math.cos(phi) + x0, slip*math.cos(theta)*math.sin(phi) + y0, slip*math.sin(theta) + z0]]
			block_info[block_id]['positions'] += [[event_time, slip*math.cos(theta)*math.cos(phi) + x0, slip*math.cos(theta)*math.sin(phi) + y0, slip*math.sin(theta) + z0, slip]]
			#
			if i_rw%10**5==0:
				print 'rw: %d (dt=%f)' % (i_rw, time.time()-t0)
		
		#
	print "finished spinning all the sweep events: ", time.time()-t0
	#
	t0=time.time()
	print "converting positions to recarrays."
	#
	# and for convenience, convert ['positions'] to a recarray:
	for key in block_info.iterkeys():
		# outputs = numpy.core.records.fromarrays(zip(*outputs), names=output_names, formats = [type(x).__name__ for x in outputs[0]])
		#block_info[key]['positions'] = numpy.core.records.fromarrays(zip(*block_info[key]['positions']), names=['block_id', 'event_year', 'x', 'y', 'z'], formats = [type(x).__name__ for x in block_info[key]['positions'][0]] )
		pos_col_names = ['event_year', 'x', 'y', 'z', 'slip']
		block_info[key]['positions'] = numpy.core.records.fromarrays(zip(*block_info[key]['positions']), names=pos_col_names, formats = [type(x).__name__ for x in block_info[key]['positions'][0]] )
	#
	print "finished: %f" % (time.time()-t0)
	#
	print "try to pickle:"
	try:
		with open(f_pickle_out, 'w') as f_pickle:
			cPickle.dump(block_info, f_pickle)
		#
	except:
		print "failed to pickle..."
	
	return block_info
#
def plot_blockwise_slip(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_start=0, i_stop=None, do_return=False, fnum=0, sim_file=default_sim_file, map_size=[8.,10.]):
	# eventually, add section and faultwise filters...
	#
	#
	n_cpus=None
	plt.ion()
	plt.figure(fnum)
	plt.clf()
	#
	# blockwise_obj is a dict (or dict-like) object, with keys: BWS[section_id]
	if isinstance(blockwise_obj, str):
		blockwise_obj = numpy.load(blockwise_obj)
	#return blockwise_obj
	if sections==None:
		sections = blockwise_obj.keys()
	#
	# draw a map:
	ll_range = vc_parser.get_fault_model_extents(section_ids=sections, sim_file=sim_file, n_cpus=n_cpus)
	#
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	#
	print "make map..."
	plt.figure(fnum, figsize=map_size)
	#bm = vc_basemap(llcrnrlon=ll_range['lon_min'], llcrnrlat=ll_range['lat_min'], urcrnrlon=ll_range['lon_max'], urcrnrlat=ll_range['lat_max'], lon_0=lon_0, lat_0=lat_0, resolution='i', projection='cyl')
	#bm = vc_parser.vc_basemap( projection='cyl', llcrnrlon=ll_range['lon_min'], llcrnrlat=ll_range['lat_min'], urcrnrlon=ll_range['lon_max'], urcrnrlat=ll_range['lat_max'], lon_0=lon_0, lat_0=lat_0, resolution='l')
	#
	print "map drawn. now, plot fault positions..."
	#plot_initial_section_positions(blockwise_obj=blockwise_obj, sections=sections, faults=faults, i_range=[i_start, i_start+1], fignum=fnum)
	#
	for j, key in enumerate(sections):
		posis = blockwise_obj[key]['positions']
		this_color = colors_[j%len(colors_)]
		#plt.plot(bm(posis['x'][0:1], posis['y'][0:1]), '.', alpha=.6, color=this_color, zorder=15)
		#plt.plot(bm(posis['x'][i_start:(len(posis) or i_stop)], posis['y'][i_start:(len(posis) or i_stop)]), '-', alpha=.3, color=this_color, zorder=15)
		plt.plot(posis['x'][0:1], posis['y'][0:1], '.', alpha=.6, color=this_color, zorder=15)
		plt.plot(posis['x'][i_start:(len(posis) or i_stop)], posis['y'][i_start:(len(posis) or i_stop)], '-', alpha=.3, color=this_color, zorder=15)
		#
	#
	if do_return: return blockwise_obj
#
def plot_initial_section_positions(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_range=[0,1], fignum=3):
	# vceventually, add section and faultwise filters...
	blockwise_obj=('dumps/blockwise_slip.pkl' or blockwise_obj)
	#
	plt.figure(fignum)
	plt.clf()
	#
	# blockwise_obj is a dict (or dict-like) object, with keys: BWS[section_id]
	if isinstance(blockwise_obj, str):
		blockwise_obj = numpy.load(blockwise_obj)
	#for i, rw in blockwise_obj.iteritems():
	for key in blockwise_obj.iterkeys():
		#rw = blockwise_obj[key]
		#x,y = rw['positions'][i_range[0]:i_range[1]]['x'], rw['positions'][i_range[0]:i_range[1]]['y']]
		#
		#plt.plot([x], [y], '.')
		plt.plot(blockwise_obj[key]['positions'][i_range[0]:i_range[1]]['x'], blockwise_obj[key]['positions'][i_range[0]:i_range[1]]['y'], '.-')
	#
	return blockwise_obj
#

