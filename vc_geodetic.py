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
def blockwise_slip_mpp(sim_file=default_sim_file, faults=None, sections=None, pipe=None, n_cpus=None, chunk_size=None):
	if n_cpus==None: n_cpus = mpp.cpu_count()
	#
	pool = mpp.Pool(n_cpus)
	#
	with h5py.File(sim_file, 'r') as vc_data:
		
		if chunk_size==None:
			chunk_size=len(vc_data['event_sweep_table'])/n_cpus
		#
		print "map to pool..."
		results = pool.map_async(get_block_slip, vc_data['event_sweep_table'], chunksize=chunk_size) # ... but i think this does not pickle
																									# because of the hdf5 file reference?
		pool.join()
		#slip_data = []
		#for rw in vc_data['event_sweep_table']:
		#	slip_data += [get_block_slip(rw)]
	#
	slip_data = results.get()
	slip_data.sort(key = lambda x: (x[0], x[1]))	# note: this should be smarter. maybe first transform to rec-array?
	#
	return slip_data
#
def blockwise_slip(sim_file=default_sim_file, faults=None, sections=None, pipe=None):
	print "getting blocks_dict..."
	block_info = get_blocks_dict(sim_file=sim_file, faults=faults, sections=sections)
	#
	print "block info fetched. assign mean values to blocks."
	#
	# add mean position to block_info:
	for key in block_info.keys():
		rw=block_info[key]
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
	print "block info ready. now sweep events."
	#
	with h5py.File(sim_file, 'r') as vc_data:
		for rw in vc_data['event_sweep_table']:
			#
			# we'll need to get the direction of motion... which nominally will need to be encoded into block_info.
			# for now, assume everything moves in the same direction, phi=0., theta=pi/2. (see above)
			block_id=rw['block_id']
			event_number = rw['event_number']
			event_time = vc_data['event_table'][event_number]['event_year']
			slip = rw['slip']
			#
			theta = block_info[block_id]['slip_theta']
			phi   = block_info[block_id]['slip_phi']
			#
			x0=block_info[key]['positions'][-1][1]
			y0=block_info[key]['positions'][-1][2]
			z0=block_info[key]['positions'][-1][3]
			#
			block_info[block_id]['positions'] += [[block_id, event_time, slip*math.cos(theta)*math.cos(phi) + x0, slip*math.cos(theta)*math.sin(phi) + y0, slip*math.sin(theta) + z0]]
		
		#
	return block_info
#
def get_block_slip(block_id_rw, sim_file=default_sim_file,):
	# in particular, for use with mpp...
	#
	with h5py.File(sim_file, 'r') as vc_data:
		block_id=block_id_rw['block_id']
		event_number = block_id_rw['event_number']
		event_time = vc_data['event_table'][event_number]['event_year']
		slip = block_id_rw['slip']
		#
		theta = 0.0
		phi   = math.pi/2.0
		#
		#x0=block_info[key]['positions'][-1][1]
		#y0=block_info[key]['positions'][-1][2]
		#z0=block_info[key]['positions'][-1][3]
		#
	#
	return [block_id, event_time, slip*math.cos(theta)*math.cos(phi), slip*math.cos(theta)*math.sin(phi), slip*math.sin(theta)]
		

