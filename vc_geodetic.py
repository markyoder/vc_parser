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
def block_slip_for_mpp(rws, block_info, events_data, plot_factor=1.0):
	return_positions = {key:[key, 0., 0., 0., 0., 0.] for key in block_info.keys()}
	for i_rw, rw in enumerate(rws):
		#
		# we'll need to get the direction of motion... which nominally will need to be encoded into block_info.
		# for now, assume everything moves in the same direction, phi=0., theta=pi/2. (see above)
		block_id=rw['block_id']
		event_number = rw['event_number']
		#event_time = vc_data['event_table'][event_number]['event_year']
		event_time = events_data[event_number]
		#
		slip = rw['slip']*plot_factor
		#
		theta = block_info[block_id]['slip_theta']
		phi   = block_info[block_id]['slip_phi']
		#
		#x0, y0, z0 = block_info[key]['positions'][-1][1:4]
		x0, y0, z0 = return_positions[key][-1][2:5]
		slip_vector = [slip*x for x in block_info['block_id']['slip_vector']]	# this should be a numpy array, but [] will work with any iteratable obj.
		#
		return_positions[block_id] += [[block_id, event_time, x0 + slip_vector[0], y0 + slip_vector[1], z0 + slip_vector[2], slip]]
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
		#block_info[key]['slip_phi'] = 0.
		#block_info[key]['slip_theta'] = math.pi/2.
		#
		dx, dy, dz = [(rw['m_%s_pt4' % xyz] - rw['m_%s_pt1' % xyz] + rw['m_%s_pt3' % xyz] - rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		fault_phi = math.atan(dy/dx)		# strike angle...
		vector_len = math.sqrt(dx*dx + dy*dy + dz*dz)		# ... though note that dz will always be zero, given the layout of the blocks.
		slip_strike_vector = math.cos(block_info[key]['rake_rad'])*numpy.array([x/vector_len for x in (dx, dy, dz)])	# (aka, normal-vector compon. * strike_compon)
		#
		#dx = -rw['m_x_pt4'] - rw['m_x_pt1'] + rw['m_x_pt3'] + rw['m_x_pt2']
		#dy = -rw['m_y_pt4'] - rw['m_y_pt1'] + rw['m_y_pt3'] + rw['m_y_pt2']
		#dz = -rw['m_z_pt4'] - rw['m_z_pt1'] + rw['m_z_pt3'] + rw['m_z_pt2']
		dx, dy, dz = [(-rw['m_%s_pt4' % xyz] - rw['m_%s_pt1' % xyz] + rw['m_%s_pt3' % xyz] + rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		fault_theta = math.atan(math.sqrt(dx*dx + dy*dy)/dz)		# and this should be equal to the dip angle (does not change when we rotate throught strike)
		vector_len = math.sqrt(dx*dx + dy*dy + dz*dz)
		slip_thrust_vector = math.sin(block_info[key]['rake_rad'])*numpy.array([x/vector_len for x in (dx, dy, dz)])	# len=1 vector in thrust direction * rake compon.
		#
		block_info[key]['slip_theta'] = fault_theta		# temporarily for diagnostics. dip angles should be equivalent... right?
		block_info[key]['slip_phi']   = fault_phi
		# ... and slip vector:
		block_info[key]['slip_vector'] = slip_strike_vector + slip_thrust_vector
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
				block_info[key]['positions'].sort(order='event_time')
				
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
		# each block consists of 4 (x,y,z) vertices. we can determine the strike and dip from this plane.
		# it looks like a typical block is like {1:TL, 2:BL, 3:BR, 4:TR}. we can also solve this more generally via planar geometry,
		# or we could use a little trick of mean \delta{x,y,z} to solve for directions which requires only that pt. 1 be
		# considered the "origin", but this can create weird artifacts when elements are tilted. for now, just assume geometry.
		#
		# strike is mean surface angle (spherical phi) between (1,4) and (2,3).
		# dip is mean vertical angle (spherical theta) between (1,2) and (4,3) (or maybe opposite of that?)
		# so, calc dx, dy, dz for both pairs and then spherical coords accordingly.
		#dx = rw['m_x_pt4'] - rw['m_x_pt1'] + rw['m_x_pt3'] - rw['m_x_pt2']
		#dy = rw['m_y_pt4'] - rw['m_y_pt1'] + rw['m_y_pt3'] - rw['m_y_pt2']
		#dz = rw['m_z_pt4'] - rw['m_z_pt1'] + rw['m_z_pt3'] - rw['m_z_pt2']
		dx, dy, dz = [(rw['m_%s_pt4' % xyz] - rw['m_%s_pt1' % xyz] + rw['m_%s_pt3' % xyz] - rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		fault_phi = math.atan(dy/dx)		# strike angle...
		# strike normal vector. we could use this + a dot-product to get slip, or we can use the angle (above).
		vector_len = math.sqrt(dx*dx + dy*dy + dz*dz)		# ... though note that dz will always be zero, given the layout of the blocks.
		slip_strike_vector = math.cos(block_info[key]['rake_rad'])*numpy.array([x/vector_len for x in (dx, dy, dz)])	# (aka, normal-vector compon. * strike_compon)
		#
		#dx = -rw['m_x_pt4'] - rw['m_x_pt1'] + rw['m_x_pt3'] + rw['m_x_pt2']
		#dy = -rw['m_y_pt4'] - rw['m_y_pt1'] + rw['m_y_pt3'] + rw['m_y_pt2']
		#dz = -rw['m_z_pt4'] - rw['m_z_pt1'] + rw['m_z_pt3'] + rw['m_z_pt2']
		
		#dx, dy, dz = [(rw['m_%s_pt4' % xyz] + rw['m_%s_pt1' % xyz] - rw['m_%s_pt3' % xyz] - rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		#fault_theta = math.atan(math.sqrt(dx*dx + dy*dy)/dz)		# and this should be equal to the dip angle (does not change when we rotate throught strike)
		# thrust normal vector. we could use this + a dot-product to get slip, or we can use the angle (above).
		#vector_len = math.sqrt(dx*dx + dy*dy + dz*dz)
		#slip_thrust_vector = math.sin(block_info[key]['rake_rad'])*numpy.array([x/vector_len for x in (dx, dy, dz)])	# len=1 vector in thrust direction * rake component.
		
		# ... except, and crap, that in the model, the dip angle does not appear to be
		# accurately portrayed by the segment blocks. we'll have to take it from the 
		# dip_angle field value.
		#		
		# note that we can also calculate this from linear transformations:
		# 1) assume slip in the y^ direction. 2) rotate theta_rake about x^,  then 3)
		# theta_dip about y^, then theta_strike about z^ (rake, dip in block_info).
		# of course, these will need to be phase adjusted, but this is the idea.
		slip_thrust_vector = rotate_y([1., 0., 0.], block_info[key]['rake_rad'])
		slip_thrust_vector = rotate_x(slip_thrust_vector, block_info[key]['dip_rad'])
		slip_thrust_vector = rotate_z(slip_thrust_vector, fault_phi)
		dx = 1.0*math.cos(block_info[key]['dip_rad'])*math.cos(fault_phi)
		dy = 1.0*math.cos(block_info[key]['dip_rad'])*math.sin(fault_phi)
		dz = 1.0*math.sin(block_info[key]['dip_rad'])
		slip_thrust_vector = numpy.array([dx, dy, dz]) - numpy.array(slip_strike_vector)
		#
		#fault_phi   = - .75*math.pi
		#fault_theta = - .5*math.pi
		#
		# ... actually, this is all wrong (the right approach to calculations, but not the right angles).
		# now, these angles may need to be adjusted for direction and origin.
		# these should be the actual, spherical coords, direction of slip.

		#block_info[key]['slip_theta'] = fault_theta		# temporarily for diagnostics. dip angles should be equivalent... right?
		block_info[key]['slip_phi']   = fault_phi
		# ... and slip vector:
		block_info[key]['slip_vector'] = slip_strike_vector + slip_thrust_vector
		
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
			#
			slip = rw['slip']*plot_factor
			#
			#theta = block_info[block_id]['slip_theta']
			theta = block_info[block_id]['dip_rad'] + math.pi/2.
			phi   = block_info[block_id]['slip_phi']
			#
			#x0=block_info[key]['positions'][-1][1]
			#y0=block_info[key]['positions'][-1][2]
			#z0=block_info[key]['positions'][-1][3]
			x0, y0, z0 = block_info[block_id]['positions'][-1][1:4]
			#
			#block_info[block_id]['positions'] += [[block_id, event_time, slip*math.cos(theta)*math.cos(phi) + x0, slip*math.cos(theta)*math.sin(phi) + y0, slip*math.sin(theta) + z0]]
			# ... no, this is not quite it either. we have to do the rotation transform
			# properly.
			#block_info[block_id]['positions'] += [[event_time, slip*math.cos(block_info['rake_rad'])*math.cos(phi) + x0, slip*math.cos(block_info['rake_rad'])*math.sin(phi) + y0, slip*math.sin(theta) + z0, slip]]
			slip_vector = [slip*x for x in block_info[block_id]['slip_vector']]	# this should be a numpy array, but [] will work with any iteratable obj.
			block_info[block_id]['positions'] += [[event_time, x0 + slip_vector[0], y0 + slip_vector[1], z0 + slip_vector[2], slip]]
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
def plot_blockwise_slip(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_start=0, i_stop=None, do_return=False, fnum=0, sim_file=default_sim_file, map_size=[8.,10.], map_res='i', map_padding = .7):
	# eventually, add section and faultwise filters...
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
		block_ids = blockwise_obj.keys()
	else:
		with h5py.File(sim_file, 'r') as vch5:
			blocks = vch5['block_info_table']
			block_ids = [rw['block_id'] for rw in blocks if rw['section_id'] in sections] 
	#
	# draw a map:
	ll_range = vc_parser.get_fault_model_extents(section_ids=sections, sim_file=sim_file, n_cpus=n_cpus)
	#
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	#
	print "make map..."
	plt.figure(fnum, figsize=map_size)
	plt.clf()
	#bm = vc_basemap(llcrnrlon=ll_range['lon_min'], llcrnrlat=ll_range['lat_min'], urcrnrlon=ll_range['lon_max'], urcrnrlat=ll_range['lat_max'], lon_0=lon_0, lat_0=lat_0, resolution='i', projection='cyl')
	#bm = lambda *args: args	# use this if there is no map...
	#
	bm = vc_parser.vc_basemap( projection='cyl', llcrnrlon=ll_range['lon_min']-map_padding, llcrnrlat=ll_range['lat_min']-map_padding, urcrnrlon=ll_range['lon_max']+map_padding, urcrnrlat=ll_range['lat_max']+map_padding, lon_0=lon_0, lat_0=lat_0, resolution=map_res)
	#
	print "map drawn. now, plot fault positions..."
	#plot_initial_section_positions(blockwise_obj=blockwise_obj, sections=sections, faults=faults, i_range=[i_start, i_start+1], fignum=fnum)
	#
	for j, key in enumerate(block_ids):
		#posis = blockwise_obj[key]['positions']
		posis = blockwise_obj[key]['positions'][(0 or i_start)]
		posis = numpy.append(posis, blockwise_obj[key]['positions'][(len(blockwise_obj[key]['positions']) or i_stop)-1])
		this_color = colors_[j%len(colors_)]
		# 
		#plt.plot(bm(posis['x'][0], posis['y'][0]), '.', alpha=.6, color=this_color, zorder=15)
		#plt.plot(bm(posis['x'], posis['y']), '-', alpha=.3, color=this_color, zorder=15)
		#
		# xy_to_lat_lon(x, y, sim_file=allcal_full_mks, lat0=None, lon0=None, chi=111.1, return_format='dict')
		X,Y = posis['x'][0], posis['y'][0]
		Y,X = vc_parser.xy_to_lat_lon(X, Y, sim_file=default_sim_file, return_format='list')
		X,Y = bm(X,Y)
		plt.plot([X], [Y], '.', alpha=.6, color=this_color, zorder=15)
		#
		#plt.plot(bm(posis['x'][0].tolist(), posis['y'][0].tolist()), '.', alpha=.6, color=this_color, zorder=15)
		X,Y = posis['x'], posis['y']
		Y,X = zip(*[vc_parser.xy_to_lat_lon(X[i],Y[i], sim_file=default_sim_file, return_format='list') for i in xrange(len(X))])
		#plt.plot(bm(posis['x'], posis['y']), '-', alpha=.3, color=this_color, zorder=15)
		plt.plot(X,Y, '-', alpha=.3, color=this_color, zorder=15)
		#
		#
		if j%10**5==0: print "plotted %d sections..." % j
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
#
def plot_blockwise_slip_3d(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_start=0, i_stop=None, do_return=False, fnum=0, sim_file=default_sim_file, fig_size=(10., 8.)):
	# eventually, add section and faultwise filters...
	#
	n_cpus=None
	plt.ion()
	f=plt.figure(fnum, figsize=fig_size)
	plt.clf()
	ax3d = f.add_subplot(111, projection='3d')
	#
	# blockwise_obj is a dict (or dict-like) object, with keys: BWS[section_id]
	if isinstance(blockwise_obj, str):
		blockwise_obj = numpy.load(blockwise_obj)
	#return blockwise_obj
	if sections==None:
		block_ids = blockwise_obj.keys()
	else:
		with h5py.File(sim_file, 'r') as vch5:
			blocks = vch5['block_info_table']
			block_ids = [rw['block_id'] for rw in blocks if rw['section_id'] in sections] 
	#
	ax3d.set_xlabel('$x$')
	ax3d.set_ylabel('$y$')
	ax3d.set_zlabel('$z$')
	for j, key in enumerate(block_ids):
		#posis = blockwise_obj[key]['positions']
		posis = blockwise_obj[key]['positions'][(0 or i_start)]
		posis = numpy.append(posis, blockwise_obj[key]['positions'][(len(blockwise_obj[key]['positions']) or i_stop)-1])
		this_color = colors_[j%len(colors_)]
		#
		ax3d.plot([posis['x'][0]], [posis['y'][0]], [posis['z'][0]], marker='.', alpha=.3, color=this_color)
		ax3d.plot(posis['x'], posis['y'], posis['z'], linestyle='-', alpha=.3, color=this_color)
		#
		#if j%10**5==0: print "plotted %d sections..." % j
	#
	if do_return: return blockwise_obj

#
def rotation_matrix_general(axis=None, theta=None):
	# generalize rotation matrix from:
	# http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
	#
	"""
	#Return the rotation matrix associated with counterclockwise rotation about
	#the given axis by theta radians.
    """
	axis = numpy.asarray(axis)
	theta = numpy.asarray(theta)
	axis = axis/math.sqrt(numpy.dot(axis, axis))
	a = math.cos(theta/2.)
	b, c, d = -axis*math.sin(theta/2.)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	
	return numpy.array([[aa+bb-cc-dd, 2.*(bc+ad), 2.*(bd-ac)],
			[2.*(bc-ad), aa+cc-bb-dd, 2.*(cd+ab)],
			[2.*(bd+ac), 2.*(cd-ab), aa+dd-bb-cc]])
#
def rotate_vector_general(vector=None, axis=None, theta=None):
	return numpy.dot(rotation_matrix_general(axis, theta), vector)
#
def rotate_x(vector=None, theta=0.):
	# probably ought to code these directly for speed. also, note the distinction between
	# vector and coordinate rotations.
	rotation_matrix_x = numpy.array([[1.0, 0., 0.], [0., math.cos(theta), -math.sin(theta)], [0., math.sin(theta), math.cos(theta)]])
	#
	return numpy.dot(rotation_matrix_x, vector)
	#return rotate_vector_general(vector=vector, axis=[1., 0., 0.], theta=theta)
#
def rotate_y(vector=None, theta=0.):
	# rotate about y axis (again, should hard-code these for speed. is there a numpy implementation?)
	#
	#M_y_rotation = numpy.array([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-math.sin(theta), 0., math.cos(theta)]]
	rotation_matrix_y = numpy.array([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-math.sin(theta), 0., math.cos(theta)]])
	#
	return numpy.dot(rotation_matrix_y, vector)
	#return rotate_vector_general(vector=vector, axis=[0., 1., 0.], theta=theta)
#
def rotate_z(vector=None, theta=0.):
	# ... and we should probably pre-define these, for speed...
	rotation_matrix_z = numpy.array([[math.cos(theta), -math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 1., 0.]])
	#
	return numpy.dot(rotation_matrix_z, vector)
	#return rotate_vector_general(vector=vector, axis=[0., 0., 1.], theta=theta)
