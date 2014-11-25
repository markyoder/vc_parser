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
from operator import itemgetter
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
import quakelib
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
year2sec = 3600.*24.*365.24

#
class meta_dict(dict):
	# a simple wrapper class to include some meta-data with a dictionary.
	# aka, return a dictionary containing indices, coordinates, and values for a displacement field.
	# also include the meta-data from which the data were generated. the object will, for most purposes,
	# behave just like a dictionary, and with minimal allowances, a dict. object can be substituted.
	#
	def __init__(self, dict_in=None, **kwargs):
		if dict_in==None: dict_in={}
		super(meta_dict, self).__init__(dict_in)
		#
		# ultimately, meta-data can be defined dynamically (it does not need to be pre-defined here).
		for key, val in kwargs.iteritems():
			setattr(self, key, val)
		
	@property
	def __doc__():
		return "a simple wrapper class to include some meta-data with a dictionary. aka, return a dictionary containing indices, coordinates, and values for a displacement field. also include the meta-data from which the data were generated. the object will, for most purposes, 	# behave just like a dictionary, and with minimal allowances, a dict. object can be substituted."
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
		event_time = events_data[event_number]['event_year']
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
		events_data = {key:{'event_year':val_year, 'event_magnitude':val_mag} for key,val_year,val_mag in zip(*[vc_data['event_table']['event_number'], vc_data['event_table']['event_year'], vc_data['event_table']['event_magnitude']])}
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
def plot_block_slip_vector(block_id=0, sim_file=default_sim_file):
	#
	# for development and testing: plot a single block and it's supposed slip direction(s).
	#
	with h5py.File(sim_file, 'r') as vc_data:
		block=vc_data['block_info_table'][block_id]
		#
		block_vertices = [[block['m_%s_pt%d' % (xyz, 1+corner_index)] for xyz in ('x', 'y', 'z')] for corner_index in xrange(4)]
		# get strike direction vector:
		dx_ph, dy_ph, dz_ph = [(block['m_%s_pt4' % xyz] - block['m_%s_pt1' % xyz] + block['m_%s_pt3' % xyz] - block['m_%s_pt2' % xyz])*.5 for xyz in ('x', 'y', 'z')]
		norm_len_strike = numpy.linalg.norm([dx_ph, dy_ph, dz_ph])
		dx_ph, dy_ph, dz_ph = [x/norm_len_strike for x in [dx_ph, dy_ph, dz_ph]]
		vec_len_strike = 1000.
		#
		# thrust direction vector:
		dx_th, dy_th, dz_th = [(-block['m_%s_pt4' % xyz] - block['m_%s_pt1' % xyz] + block['m_%s_pt3' % xyz] + block['m_%s_pt2' % xyz])*.5 for xyz in ('x', 'y', 'z')]
		norm_len_thrust = numpy.linalg.norm([dx_th, dy_th, dz_th])
		dx_th, dy_th, dz_th = [x/norm_len_thrust for x in [dx_th, dy_th, dz_th]]
		#vec_len_thrust = .5*norm_len_thrust
		vec_len_thrust = 1000.
		#vec_len = math.sqrt(dx*dx + dy*dy)
		print "vec_len: ", norm_len_strike, norm_len_thrust
		print "verts: ", dx_th, dy_th, dz_th
		dx_strike, dy_strike, dz_strike = [x*norm_len_strike for x in [dx_ph, dy_ph, dz_ph]]
		mean_x, mean_y, mean_z = [numpy.mean([block['m_%s_pt%d' % (xyz, j)] for j in [1,2,3,4]]) for xyz in ('x', 'y', 'z')]
		theta_rake = block['rake_rad'] + math.pi 	# correct for "backslip" direction.
		theta_dip = block['dip_rad']  -math.pi/2.
		#
	#
	# calculate some geodetic bits:
	fault_phi = math.atan(dy_ph/dx_ph)
	fault_theta_0 = math.acos(dz_th/math.sqrt(dx_th*dx_th + dy_th*dy_th + dz_th*dz_th))
	#
	strike_vector = numpy.array([dx_ph, dy_ph, dz_ph])*math.cos(theta_rake)
	thrust_vector = numpy.array([dx_th, dy_th, dz_th])*math.sin(theta_rake)
	#
	slip_vector_geom = strike_vector+thrust_vector
	
	print "prams: "
	print "dip: %f/%f/%f", (theta_dip, fault_theta_0, 1.0*math.pi-fault_theta_0)
	print "rake: ",theta_rake
	print "strike(ish): ", fault_phi
	#
	# now, try creating a slip vector from scratch using the various known angles...
	#v_len = vec_len
	#slip_vector_0 = [v_len, 0, 0.]
	#slip_vector = [v_len, 0, 0.]
	#slip_vector = rotate_y(slip_vector, theta_rake)
	#slip_vector_1 = slip_vector
	#slip_vector = rotate_x(slip_vector, block['dip_rad'])
	#slip_vector = rotate_z(slip_vector, fault_phi)
	#
	#thrust_phase=0.
	#dx_slip = slip_vector[0]
	#dy_slip = slip_vector[1]
	#dz_slip = slip_vector[2]
	#
	#slip_thrust_vector = numpy.array([dx, dy, dz]) - numpy.array(slip_vector)	
	#
	block_vertices += [block_vertices[0]]
	blocks_xyz = zip(*block_vertices)
	figs = []
	for i in xrange(4):
		plt.ion()
		figs+=[plt.figure(i)]
		#ax=plt.gca()
		figs[-1].clf()
		if i<3:
			plt.plot(blocks_xyz[i%3], blocks_xyz[(i+1)%3], 'o-')
			plt.xlabel('%s' % ['x','y','z'][i%3])
			plt.ylabel('%s' % ['x','y','z'][(i+1)%3])
	#
	mean_center = numpy.array([mean_x, mean_y, mean_z])
	#
	# and strike/thrust components:
	#dx_strike, dy_strike, dz_strike = [math.cos(-theta_rake)
	#
	ax3d = figs[-1].add_subplot(111, projection='3d')
	ax3d.set_xlabel('x')
	ax3d.set_ylabel('y')
	ax3d.set_zlabel('z')
	ax3d.plot(blocks_xyz[0], blocks_xyz[1], blocks_xyz[2], 'o-')
	mean_vector = numpy.array([mean_x, mean_y, mean_z])
	ax3d.plot([mean_x], [mean_y], [mean_z], 'o')
	#ax3d.plot([mean_x, mean_x+dx_strike], [mean_y, mean_y+dy_strike], [mean_z, mean_z+dz_strike], '-', lw=4)
	#ax3d.plot([mean_x, mean_x+dx_slip], [mean_y, mean_y+dy_slip], [mean_z, mean_z+dz_slip], 'd-', lw=2)
	
	#ax3d.plot([mean_x, mean_x+slip_vector_0[0]], [mean_y, mean_y+slip_vector_0[1]], [mean_z, mean_z+slip_vector_0[2]], '--', lw=2)
	#ax3d.plot([mean_x, mean_x+slip_vector_1[0]], [mean_y, mean_y+slip_vector_1[1]], [mean_z, mean_z+slip_vector_1[2]], '--', lw=2)
	
	X,Y,Z = zip(*[mean_vector, 1000.*numpy.array(strike_vector)+ mean_vector])
	ax3d.plot(X,Y,Z, '-', lw=2)
	X,Y,Z = zip(*[mean_vector, 1000.*numpy.array(thrust_vector)+ mean_vector])
	ax3d.plot(X,Y,Z, '-', lw=2)
	X,Y,Z = zip(*[mean_vector, 1000.*numpy.array(slip_vector_geom)+ mean_vector])
	ax3d.plot(X,Y,Z, '*-', lw=2)
	
	
	#
	#figs[0].plot(blocks_xyz[0], blocks_xyz[1], 'o-')
	#figs[1].plot(blocks_xyz[1], blocks_xyz[2], 'o-')
	#figs[2].plot(blocks_xyz[2], blocks_xyz[0], 'o-')
	
#
def blockwise_slip(sim_file=default_sim_file, faults=None, sections=None, pipe=None, f_pickle_out='dumps/blockwise_slip_0.pkl', plot_factor=1.0):
	# note: plot_factor is being (effectively) moved to the plotting routines, so should not typically be used here.
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
	for key, rw in block_info.iteritems():
		# this is probably a good candidate for MPP since we send over only one row (or a block of rows) and return a small dict.
		# (note there is an MPP script that calls this function itself, which is an imperfect approach).
		#
		# some of these angle and magnitude calculations probably need to be cleaned up. many of them can be found in quakelib.Simelement(),
		# and the rest of them probably should be moved to that class (see .dip(); we can add strike(), etc.)
		#
		#mean_x = numpy.mean([rw['m_x_pt%d' % j] for j in [1,2,3,4]])
		#mean_y = numpy.mean([rw['m_y_pt%d' % j] for j in [1,2,3,4]])
		#mean_z = numpy.mean([rw['m_z_pt%d' % j] for j in [1,2,3,4]])
		# in a single list comprehension:
		#
		mean_x, mean_y, mean_z = [numpy.mean([rw['m_%s_pt%d' % (xyz, j)] for j in [1,2,3,4]]) for xyz in ('x', 'y', 'z')]
		#
		# slip angles: get geodetic information about the faults from here.
		#
		# each block consists of 4 (x,y,z) vertices. we can determine the strike and dip from this plane.
		# it looks like a typical block is like {1:TL, 2:BL, 3:BR, 4:TR}. because rake angle is (or should be) defined
		# absolutely with respect to right-lateral or left-lateral faulting (so rake=0 is (say) "right" and rake=pi is "left"),
		# it might be better to define the fault geometry axes absolutely (aka the strike vector is always positive towards east and
		# the orthogonal thrust vector is always up.
		#
		# strike is mean surface angle (spherical phi) between (1,4) and (2,3).
		# dip is mean vertical angle (spherical theta) between (1,2) and (4,3) (or maybe opposite of that?)
		# ... but, it looks like the fault segments to not graphically represent dip, so we'll take that from block_info_table.
		# so, calc dx, dy, dz for both pairs and then spherical coords accordingly.
		#dx = rw['m_x_pt4'] - rw['m_x_pt1'] + rw['m_x_pt3'] - rw['m_x_pt2']
		#dy = rw['m_y_pt4'] - rw['m_y_pt1'] + rw['m_y_pt3'] - rw['m_y_pt2']
		#dz = rw['m_z_pt4'] - rw['m_z_pt1'] + rw['m_z_pt3'] - rw['m_z_pt2']
		#dx, dy, dz 
		strike_vector = [(rw['m_%s_pt4' % xyz] - rw['m_%s_pt1' % xyz] + rw['m_%s_pt3' % xyz] - rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		fault_phi = math.atan(strike_vector[1]/strike_vector[0])		# strike angle (plus phase) -- specifically "strike" from x axis (east), so strike + pi/2.
		strike_len = numpy.linalg.norm(strike_vector)
		strike_vector=numpy.array(strike_vector)/strike_len	# direction of strike
		#
		# get thrust vector:
		#dx_th, dy_th, dz_th 
		thrust_vector = [(-rw['m_%s_pt4' % xyz] - rw['m_%s_pt1' % xyz] + rw['m_%s_pt3' % xyz] + rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		dx_th, dy_th, dz_th = thrust_vector
		# spherical transformation:
		fault_theta = math.acos(dz_th/math.sqrt(dx_th*dx_th + dy_th*dy_th + dz_th*dz_th))
		thrust_len = numpy.linalg.norm(thrust_vector)
		thrust_vector = numpy.array(thrust_vector)/thrust_len
		#
		# actual slip in the strike and thrust directions:
		strike_slip_vector = math.cos(block_info[key]['rake_rad'])*strike_vector
		thrust_slip_vector = math.sin(block_info[key]['rake_rad'])*thrust_vector
		#total_slip_vector = strike_slip_vector + thrust_slip_vector
		#
		#block_info[key].update({'mean_x':mean_x, 'mean_y':mean_y, 'mean_z':mean_z})
		block_info[key]['slip_phi']   = fault_phi
		# ... and slip vector:
		#block_info[key]['slip_vector'] = slip_strike_vector + slip_thrust_vector
		block_info[key]['slip_vector'] = strike_slip_vector + thrust_slip_vector		# total slip vector
		#
		# and set a field for a slip-sequence.
		block_info[key]['positions'] = [[0.0, mean_x, mean_y, mean_z, 0., 'initial']]	# [[time, x,y,z, slip, slip_type]]
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
		#events_data = {key:val for key,val in zip(*[vc_data['event_table']['event_number'], vc_data['event_table']['event_year']])}
		#events_data = {key:val for key,val in zip(*[vc_data['event_table']['event_number'], {'event_year':vc_data['event_table']['event_year'], 'event_magnitude':vc_data['event_table']['event_magnitude']}])}		
		events_data = {key:{'event_year':val_year, 'event_magnitude':val_mag} for key,val_year,val_mag in zip(*[vc_data['event_table']['event_number'], vc_data['event_table']['event_year'], vc_data['event_table']['event_magnitude']])}
		#
		print "events are indexed: %f" % (time.time()-t0)
		t0-time.time()
		#
		for i_rw, rw in enumerate(vc_data['event_sweep_table']):
			#
			# we'll need to get the direction of motion... which nominally will need to be encoded into block_info.
			block_id=rw['block_id']
			event_number = rw['event_number']
			#
			# (see above definition for events_data: events_data ~ {'event_number':{event_year:yr, 'event_magnitude':mag}}
			event_time = events_data[event_number]['event_year']
			#
			slip = rw['slip']*plot_factor
			#
			phi   = block_info[block_id]['slip_phi']		# more or less strike angle (maybe exactly strike angle); angle on surface about z^.
			#
			# now, we want to add two slip vectors. 1) the aseismic_slip_vector = mean_slip_rate * \delta_t, 2) the seismic_slip_vector = earthquake_slip.
			#
			seismic_slip_vector = [slip*x for x in block_info[block_id]['slip_vector']]	# this should be a numpy array, but [] will work with any iteratable obj.
			#																	# (note that ['slip_vector'] is more correctly a slip_unit_vector.
			# calc aseismic slip:
			#delta_t = event_time - block_info[block_id]['positions'][-1][0]	# (we've not yet converted to a rec_array). native value is in m/s; convert to m/year
			aseismic_slip = block_info[block_id]['slip_velocity']*year2sec*(event_time - block_info[block_id]['positions'][-1][0])
			aseismic_slip_vector = [aseismic_slip*x for x in block_info[block_id]['slip_vector']]
			#
			# and just to be thorough, let's actually calculate the rupture duration to temporally separate the aseismic and seismic data pionts:
			rupture_delta_t = (10.0**(events_data[event_number]['event_magnitude']/2. - 2.3))/year2sec		# (see Yoder et al. 2014 ETAS). the exact value is not really 			
																					# important; it just needs to be small and preferably (??) not 0.
			# aseismic:
			x0, y0, z0 = block_info[block_id]['positions'][-1][1:4]
			block_info[block_id]['positions'] += [[event_time, x0 + aseismic_slip_vector[0], y0 + aseismic_slip_vector[1], z0 + aseismic_slip_vector[2], aseismic_slip, 'aseismic']]
			#
			# seismic (this is the new position after the event):
			x0, y0, z0 = block_info[block_id]['positions'][-1][1:4]
			block_info[block_id]['positions'] += [[event_time+rupture_delta_t, x0 + seismic_slip_vector[0], y0 + seismic_slip_vector[1], z0 + seismic_slip_vector[2], slip, 'seismic']]
			
			#block_info[block_id]['positions'] += [[event_time, x0 + slip_vector[0], y0 + slip_vector[1], z0 + slip_vector[2], slip]]
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
		pos_col_names = ['event_year', 'x', 'y', 'z', 'slip', 'slip_type']
		type_len = max([len(x[-1]) for x in block_info[key]['positions']])	# max length of any value in slip_type col. (needed for recarray typing) 
		block_info[key]['positions'] = numpy.core.records.fromarrays(zip(*block_info[key]['positions']), names=pos_col_names, formats = [type(x).__name__ for x in block_info[key]['positions'][0][0:-1]] + ['S%d' % (type_len)] )
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
def plot_blockwise_slip(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_start=0, i_stop=None, do_return=False, fnum=0, sim_file=default_sim_file, map_size=[8.,10.], map_res='i', map_padding = .7, plot_factor=100.0):
	# original "blockwise" map plotter. plots a 2D map; each block element is a separate plot (and so goes the color rotation).
	# originally, this plotted the full slip time-series -- which is fine by itself, but when we add the map, it fully explodes, so we
	# plot only the total slip vector (x[0], x[-1]).
	#
	# eventually, add section and faultwise filters...
	#
	n_cpus=None
	plt.ion()
	plt.figure(fnum)
	plt.clf()
	#
	if plot_factor==None: plot_factor = 1.0
	plot_factor = float(plot_factor)
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
	print "ll_range: ", ll_range
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
		#delta_X, delta_Y, delta_Z = [posis[xyz][-1]-posis[xyz][0] for xyz in ['x', 'y', 'z']]
		X,Y = posis['x'], posis['y']
		Y,X = zip(*[vc_parser.xy_to_lat_lon(X[i],Y[i], sim_file=default_sim_file, return_format='list') for i in xrange(len(X))])
		delta_X = (X[-1] - X[0]) * plot_factor
		delta_Y = (Y[-1] - Y[0]) * plot_factor
		#delta_z = posis['z'][-1] - posis['z'][0]
		#plt.plot(bm(posis['x'], posis['y']), '-', alpha=.3, color=this_color, zorder=15)
		#plt.plot(bm(posis['x'], posis['y']), '-', alpha=.3, color=delta_z, cmap=cm.coolwarm, zorder=15)
		#plt.plot(X,Y, '-', alpha=.3, color=this_color, zorder=15)
		#
		plt.arrow(X[0], Y[0], delta_X, delta_Y, head_width=.005, head_length=.01, fc=this_color, ec=this_color, zorder=15, alpha=.5)
		#
		#
		#if j%100==0: print "plotted %d sections..." % j
	#
	if do_return: return blockwise_obj
#
#
def plot_blockwise_slip_color_thrust(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_start=0, i_stop=None, do_return=False, fnum=0, sim_file=default_sim_file, map_size=[8.,10.], map_res='i', map_padding = .7, plot_factor=1.0):
	# eventually, add section and faultwise filters...
	#
	n_cpus=None
	
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
	print "ll_range: ", ll_range
	#
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	#
	print "make map..."
	#set up the map:
	#
	plt.ion()
	plt.figure(fnum, figsize=map_size)
	plt.clf()
	bm = vc_parser.vc_basemap( projection='cyl', llcrnrlon=ll_range['lon_min']-map_padding, llcrnrlat=ll_range['lat_min']-map_padding, urcrnrlon=ll_range['lon_max']+map_padding, urcrnrlat=ll_range['lat_max']+map_padding, lon_0=lon_0, lat_0=lat_0, resolution=map_res)
	#
	print "map drawn. now, plot fault positions..."
	#plot_initial_section_positions(blockwise_obj=blockwise_obj, sections=sections, faults=faults, i_range=[i_start, i_start+1], fignum=fnum)
	#
	# we want the delta_z range, which means that ultimately we need two loops through the data. first, let's try looping once to collect all the plotting arrays,
	# then a second loop to plot them. this may be memory problematic, particularly if we want to plot the full sequence as opposed to just the first and last points.
	#
	slip_vectors = []
	#min_delta_z = None
	#max_delta_z = None
	#
	for j, key in enumerate(block_ids):
		#posis = blockwise_obj[key]['positions']
		posis = blockwise_obj[key]['positions'][(0 or i_start)]
		posis = numpy.append(posis, blockwise_obj[key]['positions'][(len(blockwise_obj[key]['positions']) or i_stop)-1])
		#
		k_start = (0 or i_start)
		x0,y0,z0    = blockwise_obj[key]['positions'][k_start]['x'], blockwise_obj[key]['positions'][k_start]['y'], blockwise_obj[key]['positions'][k_start]['z']
		y0, x0 = vc_parser.xy_to_lat_lon(x0,y0, sim_file=default_sim_file, return_format='list')
		x0, y0 = bm(x0,y0)
		k_stop = (len(blockwise_obj[key]['positions']) or i_stop)-1
		x1,y1,z1    = blockwise_obj[key]['positions'][k_stop]['x'], blockwise_obj[key]['positions'][k_stop]['y'], blockwise_obj[key]['positions'][k_stop]['z']
		y1, x1 = vc_parser.xy_to_lat_lon(x1,y1, sim_file=default_sim_file, return_format='list')
		x1, y1 = bm(x1, y1)
		#
		#dx, dy dz = x1-x0, y1-y0, z1-z0
		slip_vectors += [[x0,y0,z0, (x1-x0), (y1-y0), (z1-z0)]]
		#min_delta_z = min(plot_arrays[-1][-1], (min_delta_z or plot_arrays[-1][-1]))
		#max_delta_z = max(plot_arrays[-1][-1], (max_delta_z or plot_arrays[-1][-1]))
		#
	#
	slip_vectors = numpy.core.records.fromarrays(zip(*slip_vectors), names=('x','y','z','dx', 'dy', 'dz'), formats=[type(x).__name__ for x in slip_vectors[0]])
	#
	# now, plot.
	#	
	min_delta_z = min(slip_vectors['dz'])
	max_delta_z = max(slip_vectors['dz'])
	#
	# ... and eventually, split this up into two functions? sounds like a great idea, but since the coordinates are tied to the map projection, it might get
	# a little messy or feel a little backwards. we may end up, basically, with some plotting functions that don't stand alone very well. in other words, 
	# we'll pass a bunch of parameters through this funciton (plot_factor, etc.), and mabye pass a "plotter" function as a parameter?
	#
	# set up colors:
	cmap_name = 'jet'		#'coolwarm'
	my_colormap = plt.get_cmap(cmap_name)
	cNorm = mcolor.Normalize(vmin=min_delta_z, vmax=max_delta_z)
	scalar_map = cm.ScalarMappable(norm=cNorm, cmap=my_colormap)
	#print "delta_z range: ", min_delta_z, max_delta_z
	#print "scalar_map: ", scalar_map.get_clim()
	#
	if plot_factor==None: plot_factor = 1.0
	plot_factor = float(plot_factor)
	#
	for j, rw in enumerate(slip_vectors):
		plt.plot(rw['x'], rw['y'], color=scalar_map.to_rgba(rw['dz']), alpha=.5)
		plt.arrow(rw['x'], rw['y'], rw['dx']*plot_factor, rw['dy']*plot_factor, color=scalar_map.to_rgba(rw['dz']), alpha=.5)
	#
	#if do_return: return blockwise_obj
	if do_return: return slip_vectors
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
def plot_blockwise_slip_3d(blockwise_obj='dumps/blockwise_slip.pkl', sections=None, faults=None, i_start=0, i_stop=None, do_return=False, fnum=0, sim_file=default_sim_file, fig_size=(10., 8.), plot_factor=1.0):
	# eventually, add section and faultwise filters...
	#
	n_cpus=None
	plt.ion()
	f=plt.figure(fnum, figsize=fig_size)
	plt.clf()
	ax3d = f.add_subplot(111, projection='3d')
	#
	if plot_factor==None: plot_factor = 1.0
	plot_factor = float(plot_factor)
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
		#ax3d.plot(posis['x'], posis['y'], posis['z'], linestyle='-', alpha=.3, color=this_color)
		#
		delta_X = (posis['x'][-1] - posis['x'][0]) * plot_factor
		delta_Y = (posis['y'][-1] - posis['y'][0]) * plot_factor
		delta_z = (posis['z'][-1] - posis['z'][0]) * plot_factor
		
		X = [posis['x'][0]+i*delta_X for i in xrange(2)]
		Y = [posis['y'][0]+i*delta_Y for i in xrange(2)]
		Z = [posis['z'][0]+i*delta_z for i in xrange(2)]
			
		#plt.plot(bm(posis['x'], posis['y']), '-', alpha=.3, color=this_color, zorder=15)
		#plt.plot(X,Y, '-', alpha=.3, color=this_color, zorder=15)
		#ax3d.arrow(posis['x'][0], posis['y'][0], posis['z'][0], delta_X, delta_Y, delta_z, head_width=.005, head_length=.01, fc=this_color, ec=this_color, zorder=15, alpha=.5)
		ax3d.plot(X,Y,Z, color=this_color, zorder=15, alpha=.5)
		#
		#if j%10**5==0: print "plotted %d sections..." % j
	#
	if do_return: return blockwise_obj
#
def slip_field_time_series(blockwise_obj=None, dx=None, dy=None, i_start=0, i_stop=-1, di=1, section_ids=None, plot_factor=10., f_out=None):
	# make a massive slip_field_time_series object. aka, call slip_field for some set of time-slices. this has potential to be HUGE
	# and really requires some sort of custom treatment depending on the specific question at hand, so we'll probably hold off on writing it for a while.
	#
	pass
#
def slip_field(blockwise_obj=None, dx=None, dy=None, i_start=0, i_stop=-1, section_ids=None, plot_factor=10., f_out=None):
	#+
	# use:
	# Vec<3> calc_displacement_vector(const Vec<3> location, const double c, const double dip, const double L, const double W, const double US, const double UD, const double UT, const double lambda, const double mu) throw(std::invalid_argument)
	#
	# args: blockwise_obj: a standard(ish) blockwise_obj (see blockwise_slip, etc.). trap it, however, so we can
	# also take a list of dicts, list, or recarray. and extract the [[dx, dy, dz] ] lists there-from.
	# can we/will we need to approximate from grid elements to full earthquake events. 
	#
	# ... basically, blockwise_obj should be a dictionary like {block_id:{...block info..., 'positions':[]}}
	#
	#pass
	# for our first iteration, (as per default values), assume we're plotting to total displacement of each block.
	# copy from QuakeLibEvent::calc_displacement_vector() (which we'll eventually need to get working... maybe sooner than later,
	# depending on speed) and QuakeLibOkada; then use quakelib.Okada.calc_displacement_vector() to do the Okada calculations.
	#
	# right now, positions are in x-y coordinates (meters?). in the next VC version, they'll be in lat/lon, so we'll need to sort that out.
	#
	i_start = (0 or i_start)
	i_stop = (0 or i_stop)
	#
	# a default blockwise_obj:
	if blockwise_obj==None:
		#blockwise_obj=numpy.load('dumps/blockwise_290ct.pkl')
		blockwise_obj = 'dumps/blockwise_290ct.pkl'
	#
	# blockwise_obj is a dict (or dict-like) object, with keys: BWS[section_id]
	if isinstance(blockwise_obj, str):
		blockwise_obj = numpy.load(blockwise_obj)
	#
	if section_ids==None:
		block_ids = blockwise_obj.keys()
	else:
		#with h5py.File(sim_file, 'r') as vch5:
		#	blocks = vch5['block_info_table']
		#	block_ids = [rw['block_id'] for rw in blocks if rw['section_id'] in section_ids]
		block_ids = [rw['block_id'] for rw in blockwise_obj.itervalues() if rw['section_id'] in section_ids]
	#
	block_L = numpy.linalg.norm(numpy.array([blockwise_obj[0]['m_x_pt4'], blockwise_obj[0]['m_y_pt4'], blockwise_obj[0]['m_z_pt4']]) - numpy.array([blockwise_obj[0]['m_x_pt1'], blockwise_obj[0]['m_y_pt1'], blockwise_obj[0]['m_z_pt1']]))
	block_W = numpy.linalg.norm(numpy.array([blockwise_obj[0]['m_x_pt4'], blockwise_obj[0]['m_y_pt4'], blockwise_obj[0]['m_z_pt4']]) - numpy.array([blockwise_obj[0]['m_x_pt3'], blockwise_obj[0]['m_y_pt3'], blockwise_obj[0]['m_z_pt3']]))
	if dx == None:
		#dx = block_L
		dx=1.0
	if dy==None:
		#dy=block_L
		dy=1.0
	dx = dx*block_L
	dy = dy*block_L
	#
	# for now, assume we're getting a proper "blockwise_obj"...
	try:
		# do we have a dict or recarray (aka, named columns)? we need to test more thoroughly, but this is a start.
		blockwise_obj[0]['positions'].dtype
		# and if this doesn't fail then...
		#x_index = 'x'
		#y_index = 'y'
		#z_index = 'z'
		#slip_index = 'slip'
		xyz_indices = ['x', 'y', 'z']
	except:
		# then guess...
		#x_index = 1
		#y_index = 2
		#z_index = 3
		#slip_index = 4
		xyz_indices = [1,2,3]
	#
	# try to use this derived dict. object instead (but note, you can use a regular dictionary object in place of meta_dict() ):
	# dx=None, dy=None, i_start=0, i_stop=-1, section_ids=None, plot_factor=10., f_out=None
	# (noting that the field_plot_prams dict. could alternatively be assigned as individual attributes: ...,dx=dx, dy=dy, ...
	# index this dictionary like {(i_x, i_y):{'dx':dx, 'dy':dy, 'dz':dz}, ...}, and calculate the index based on... dunno,
	# based on magnitude/L_r, but in this case, block-size, slip magnitude? use plot_factor...
	#
	try:
		disp_field = meta_dict({}, field_plot_prams = {'dx':dx, 'dy':dy, 'i_start':i_start, 'i_stop':i_stop, 'plot_factor':plot_factor, 'file_out':f_out})
	except:
		disp_field = {}
	#
	Okada_obj = quakelib.Okada()
	#
	#for block_id, block_data in blockwise_obj.iteritems():
	for j, block_id in enumerate(block_ids):
		block_data = blockwise_obj[block_id]
		#
		#dip = block_data['dip_rad']
		#strike = block_data['slip_phi']	# or some phase transformation of this...
		#print "starting block_id: ", block_id
		#
		pos_0 = numpy.array([block_data['positions'][i_start][i_x] for i_x in xyz_indices])
		pos_0_okada = [pos_0[0], pos_0[1], 0.]	# for use with Okada bits.
		pos_final = numpy.array([block_data['positions'][i_stop][i_x] for i_x in xyz_indices])
		total_slip_vector = pos_final-pos_0			# note: 1) these simulated positions are compiled in blockwise_slip(), and 2) we cast as numpy.array() here
													# becasue we'll probably want to dot-product this later (otherwise we'd use a light-n-fast list comprehension)
													# we can get relative distance components for Okada by dotting with this vector.
		slip_mag = numpy.linalg.norm(total_slip_vector)
		slip_unit_vector = block_data['slip_vector']
		#
		# calculate the field sites (and put them on a lattice):
		x0 = float(int((pos_0[0] - block_L*plot_factor)/dx))*dx
		y0 = float(int((pos_0[1] - block_L*plot_factor)/dy))*dy
		#
		x_max = pos_0[0] + block_L*plot_factor	# note: these max vals might be off-lattice, but since we use these vals as a limit (x<x_max)
		y_max = pos_0[1] + block_L*plot_factor  #     it won't matter, with the exception that our total sequence lengths will be +/- 1 element.
		#
		c = abs(block_data['m_z_pt2'])
		#
		this_element = quakelib.SimElement()
		this_element.set_rake(block_data['rake_rad'])
		#
		# set vertices:
		[this_element.set_vert(i, block_data['m_x_pt%d' % k], block_data['m_y_pt%d' % k], block_data['m_z_pt%d' % k]) for i,k in enumerate([1,2,4])]
		#
		#	print "vert: ", this_element.vert(i)
		this_element.set_lame_mu(block_data['lame_mu'])
		this_element.set_lame_lambda(block_data['lame_lambda'])
		#
		# for now, set okada vals here; get quakelib working later...
		# (calling the Okada calculations from SimElement() is not working, so calculate the Okada prams here)... for now.
		US = slip_mag*math.cos(block_data['rake_rad'])
		UD = slip_mag*math.sin(block_data['rake_rad'])
		UT= 0.0
		#
		L = numpy.linalg.norm(numpy.array([block_data['m_%s_pt4' % s] for s in ['x', 'y', 'z']]) - numpy.array([block_data['m_%s_pt1' % s] for s in ['x', 'y', 'z']]))
		W = numpy.linalg.norm(numpy.array([block_data['m_%s_pt2' % s] for s in ['x', 'y', 'z']]) - numpy.array([block_data['m_%s_pt1' % s] for s in ['x', 'y', 'z']]))
		c = abs(this_element.max_depth())
		#dip = this_element.dip()
		#
		this_y = y0
		while this_y<y_max:
			this_x = x0
			y_index = int(this_y/dy)
			#
			while this_x<x_max:
				x_index = int(this_x/dx)
				#disp_vector = numpy.array(this_x, this_y) - pos_0_okada		# displacement vector (between a site and the epicenter).
				disp_vector = get_okada_location(position_vec=[this_x, this_y, 0.], origin_vec=pos_0_okada, theta_CCW_from_x=block_data['slip_phi'])
				#
				# get okada displacement for this location:
				#positon_vector = quakelib.Vec3(*rotate_z(disp_vector, -block_data['slip_phi']))
				okada_disp = Okada_obj.calc_displacement_vector(quakelib.Vec3(*disp_vector), c, block_data['dip_rad'], block_L, block_W, slip_mag*math.cos(block_data['rake_rad']), slip_mag*math.sin(block_data['rake_rad']), 0.0, block_data['lame_lambda'], block_data['lame_mu'])
				okada_disp = rotate_z(okada_disp, block_data['slip_phi'])
				#
				# now, update disp_field values:
				if not disp_field.has_key((x_index, y_index)): disp_field[(x_index, y_index)]={'xyz':[this_x, this_y, 0.], 'd_xyz':numpy.array([0., 0., 0.])}
				disp_field[(x_index, y_index)]['d_xyz'] += numpy.array(okada_disp)
				#
				this_x += dx
			this_y += dy
		#
	if f_out!=None and isinstance(f_out, str):
		print "cPickle.dump() to: %s" % f_out
		with open(f_out, 'w') as fout:
			cPickle.dump(disp_field, fout)
		print "pickled..."
	#
	return disp_field
	#
def plot_disp_vector_field(okada_disps='dumps/okada_slips_allcal_2.pkl', fignum=0, plot_factor=1.0, z_colors=True, sim_file=default_sim_file, n_cpus=None, do_map=True):
	if isinstance(okada_disps, str):
		with open(okada_disps, 'r') as f:
			# not sure, but this might also load with numpy.load()
			okada_disps = cPickle.load(f)
		#
	#
	# the native vector-vield object is an indexed dictionary like:
	# {(x_i, y_i), {'d_xyz':[], 'xyz':[]} }
	# where x_i, y_i are lat,lon indices (bacially, truncated lat/lon integers to facitate aggregation.
	# we want to break this down into a recarray (see plot_disp_field() below)
	#
	try:
		if isinstance(okada_disps.keys()[0], tuple):
			okada_disps = [[x for x in rw['xyz']] + [dx for dx in rw['d_xyz']] for rw in okada_disps.itervalues()]
			okada_disps = numpy.rec.array(okada_disps, names=['x', 'y', 'z', 'dx', 'dy', 'dz'], formats=[type(x).__name__ for x in okada_disps[0]])
	except:
		pass
	#
	# data are like[ [x,y,z, dx, dy, dz, dxyz, dxy], ...] (where the last two are the xyz length and the xy length and,
	# because they are derived, may not always be present.
	#
	vec_colors=None
	if z_colors:
		# set up a color map based on vertical displacement:
		z_vals = okada_disps['dz']
		#log_z_vals = numpy.log10(z_vals)
		my_colormap = plt.get_cmap('jet')
		cNorm = mcolor.Normalize(vmin=min(z_vals), vmax=max(z_vals))
		scalar_map = cm.ScalarMappable(norm=cNorm, cmap=my_colormap)
		vec_colors=scalar_map.to_rgba(z_vals)
	#
	# start with a simple plot, no rules. later we'll look at color-coding the z component, etc.
	plt.figure(fignum)
	plt.ion()
	plt.clf()
	#print "plot factor: %f" % plot_factor
	#print okada_disps['dx'][0:5], plot_factor*okada_disps['dx']
	
	plt.quiver(okada_disps['x'], okada_disps['y'], [plot_factor*dx for dx in okada_disps['dx']], [plot_factor*dy for dy in okada_disps['dy']], pivot='tail', color=vec_colors, alpha=.9)
	#
	if not do_map: return okada_disps
	#
	# now, make a map:
	ll_range = vc_parser.get_fault_model_extents(section_ids=None, sim_file=sim_file, n_cpus=n_cpus)
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	map_size=[8,10]
	plt.figure(fignum+1, figsize=map_size)
	plt.clf()
	bm = vc_parser.vc_basemap( projection='cyl', llcrnrlon=ll_range['lon_min'], llcrnrlat=ll_range['lat_min'], urcrnrlon=ll_range['lon_max'], urcrnrlat=ll_range['lat_max'], lon_0=lon_0, lat_0=lat_0, resolution='i')
	
	#
	#xy_to_lat_lon(x, y, sim_file=allcal_full_mks, lat0=None, lon0=None, chi=111.1, return_format='dict')
	lat, lon = zip(*[vc_parser.xy_to_lat_lon(okada_disps['x'][i], okada_disps['y'][i], return_format='list') for i in xrange(len(okada_disps))])
	map_X, map_Y = bm(lon, lat)
	#
	# simple quiver type:
	#
	#plt.quiver(map_X, map_Y, [plot_factor*dx for dx in okada_disps['dx']], [plot_factor*dy for dy in okada_disps['dy']], pivot='tail', color=vec_colors, zorder=11, alpha=.6)
	#
	# on-map, proportioal mini-quivers:
	
	lat_tip, lon_tip = zip(*[vc_parser.xy_to_lat_lon(okada_disps['x'][i]+plot_factor*okada_disps['dx'][i], okada_disps['y'][i]+plot_factor*okada_disps['dy'][i], return_format='list') for i in xrange(len(okada_disps))])
	#
	map_X_tip, map_Y_tip = bm(lon_tip, lat_tip)
	#
	disp_norms = [numpy.linalg.norm(x) for x in zip(*[okada_disps['dx'], okada_disps['dy'], okada_disps['dz']])]
	median_displacement = numpy.median(disp_norms)
	plt.plot([x for i,x in enumerate(map_X) if disp_norms[i]>median_displacement] , [x for i,x in enumerate(map_Y) if disp_norms[i]>median_displacement], 'b.', zorder=11, alpha=.25)
	for i in xrange(len(map_X)):
		if disp_norms[i]<=median_displacement: continue
		plt.plot([map_X[i], map_X_tip[i]], [map_Y[i], map_Y_tip[i]], 'b-', zorder=11, alpha=.6)
	
	
	#
	return okada_disps

def plot_disp_field(okada_disps='dumps/okada_slips_allcal_2.pkl', plot_column = 'dxyz', fignum=0):
	# lots of details to be worked out here, not the least of which being the best return format. for slip_field(). for now, let's work with what we've got:
	# dict of dicts: {(x_index, y_index):{'d_xyz':array([dx, dy, dz]), 'xyz':array([x,y,z])}}
	# plot_column: 'dxyz', 'dx', 'dy', 'dz' are direct column plots. we'll add some code to permit other combinations like 'dxy'
	#
	# can we just load this? yeah, this works... and set default.
	if isinstance(okada_disps, str):
		with open(okada_disps, 'r') as f:
			okada_disps = cPickle.load(f)
	#
	# first, compile this into something we can plot:
	field_data = []
	for rw in okada_disps.itervalues():
		field_data += [[x for x in rw['xyz']] + [x for x in rw['d_xyz']]]
		field_data[-1] += [numpy.linalg.norm(field_data[-1][-3:])]
		field_data[-1] += [numpy.linalg.norm(field_data[-1][-2:])]
		
	field_data = numpy.core.records.fromarrays(zip(*field_data), names=['x', 'y', 'z', 'dx', 'dy', 'dz', 'dxyz', 'dxy'], formats=[type(x).__name__ for x in field_data[0]])
	field_data_prime = [x for x in field_data if x[plot_column]>0.]
	field_data_prime = numpy.core.records.fromarrays(zip(*field_data_prime), names=['x', 'y', 'z', 'dx', 'dy', 'dz', 'dxyz', 'dxy'], formats=[type(x).__name__ for x in field_data_prime[0]])
	#
	#return field_data
	#if plot_column in ('dxy', 'dyx'):
	#	z_vals = [math.sqrt(x['dx']**2. + x['dy']**2.) for x in field_data_prime]
	#else:
	z_vals = field_data_prime[plot_column]		# note: "plotting z", not vertical displacement. in this case, "z"
												# is the "norm" or absolute length of the displacement vector.
	#
	# create triang for contouring:
	#triang = tri.Triangulation(field_data['x'], field_data['y'])
	my_colormap = plt.get_cmap('jet')
	cNorm = mcolor.Normalize(vmin=min(z_vals), vmax=max(z_vals))
	scalar_map = cm.ScalarMappable(norm=cNorm, cmap=my_colormap)
	#
	# now, can we set dynamic alpha to imply relief?
	# min/max vertical displacement.
	vert_disp_min, vert_disp_max = [float(f(field_data_prime['z'])) for f in [min, max]]
	alpha_base = .5	# and we'll "shade" between alpha_base and 1.0.
	delta_vert = vert_disp_max - vert_disp_min
	alphas = [.5 + .5*(rw['z']-vert_disp_min)/delta_vert for rw in field_data_prime]
	#
	plt.figure(fignum)
	plt.clf()
	#plt.gca().set_aspect('equal')
	#plt.tricontourf(triang, field_data['dxyz'])
	#plt.colorbar()
	#
	plt.scatter(field_data_prime['x'], field_data_prime['y'], marker='.', color=scalar_map.to_rgba(z_vals),alpha=.7)
	# this doesn't quite work. we'd probably have to ploe each point indipendently, which might get memory intensive.
	#plt.scatter(field_data_prime['x'], field_data_prime['y'], marker='.', color=scalar_map.to_rgba(z_vals),alpha=alphas)
	'''
	try:
		plt.colorbar()
	except:
		print "colorbar() failed..."
	#
	'''
	#
	# and log-transformed:
	log_z_vals = [math.log10(x) for x in z_vals]
	
	my_colormap = plt.get_cmap('jet')
	cNorm = mcolor.Normalize(vmin=min(log_z_vals), vmax=max(log_z_vals))
	scalar_map = cm.ScalarMappable(norm=cNorm, cmap=my_colormap)
	plt.figure(fignum+1)
	plt.clf()
	plt.scatter(field_data_prime['x'], field_data_prime['y'], marker='.', color=scalar_map.to_rgba(log_z_vals),alpha=.7)
	return field_data	
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
	#rotation_matrix_x = numpy.array([[1.0, 0., 0.], [0., math.cos(theta), -math.sin(theta)], [0., math.sin(theta), math.cos(theta)]])
	if len(vector)==2:
		# forget about len()==1 for now...
		vector = numpy.array([0.] + vector)
	#
	#return numpy.dot(rotation_matrix_x, vector)
	#return rotate_vector_general(vector=vector, axis=[1., 0., 0.], theta=theta)
	return numpy.dot(numpy.array([[1.0, 0., 0.], [0., math.cos(theta), -math.sin(theta)], [0., math.sin(theta), math.cos(theta)]]), vector)
#
def rotate_y(vector=None, theta=0.):
	# rotate about y axis (again, should hard-code these for speed. is there a numpy implementation?)
	#theta = float(theta)%(math.pi*2.0)
	#
	#rotation_matrix_y = numpy.array([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-math.sin(theta), 0., math.cos(theta)]])
	if len(vector)==2:
		vector = numpy.array([vector[0], 0., vector[1]])
	#
	return numpy.dot(numpy.array([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-math.sin(theta), 0., math.cos(theta)]]), vector)
	#return rotate_vector_general(vector=vector, axis=[0., 1., 0.], theta=theta)
#
def rotate_z(vector=None, theta=0.):
	# ... and we should probably pre-define these, for speed...
	#rotation_matrix_z = numpy.array([[math.cos(theta), -math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]])
	if len(vector)==2: vector = numpy.append(vector, [0.])
	#
	return numpy.dot(numpy.array([[math.cos(theta), -math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]]), vector)
	#return rotate_vector_general(vector=vector, axis=[0., 0., 1.], theta=theta)
#
def draw_block_data(block_dict, fignum=4):
	# take-away message:
	# to calculate Okada:
	# 1) rotate the displacement vector by -theta_strike (aka, -block_data['slip_phi'])
	# 2) calc Okada, etc.
	# 3) rotate return vectors by block_data['slip_phi'] (aka, the pseudo-strike -- noting that this "strike" angle is measured from x^ (direction
	#    equals East), rather than y^ (North).
	#
	plt.figure(fignum)
	plt.clf()
	plt.ion()
	#
	B=block_dict	# shorthand...
	Xs = [[B['m_%s_pt%d' % (x,i)] for x in ['x','y','z']] for i in [1,2,3,4]]
	#
	zXs = zip(*Xs)
	plt.plot(zXs[0], zXs[1], '.')
	#
	vW = numpy.array(Xs[1])-numpy.array(Xs[0]) 
	vL = numpy.array(Xs[3]) - numpy.array(Xs[0])		# later versions of this will have only 3 points; this will be Xs[2] - Xs[0]
	W = numpy.linalg.norm(vW)
	L = numpy.linalg.norm(vL)
	#
	# plot the slip vector:
	plt.arrow(Xs[0][0], Xs[0][1], vL[0], vL[1], ls='dashed', color='r', head_width=25., head_length=50.)
	#
	# calculate the x and y components of the slip vector:
	theta_strike = B['slip_phi']
	x_prime = .5*L*math.cos(theta_strike)
	y_prime = .5*L*math.sin(theta_strike)
	#
	# plot along slip vector to show we got the angle (phase) right:
	plt.arrow(Xs[0][0], Xs[0][1], .5*vL[0], .5*vL[1], ls='dashed', color='m', head_width=25., head_length=50.)
	plt.plot([Xs[0][0]+x_prime], [Xs[0][1]+y_prime], 'o')
	#
	# slip-vector transformation:
	# rotate the slip vector into the X = [1,0,0] position (nominally we should be rotating the coord-system, but this is effectively the same thing...)
	x_pp, y_pp, z_pp = rotate_z([x_prime, y_prime], -theta_strike)
	plt.arrow(Xs[0][0], Xs[0][1], .5*vL[0], .5*vL[1], ls='dashed', color='m', head_width=25., head_length=50.)
	plt.arrow(Xs[0][0], Xs[0][1], x_pp, y_pp, ls='dashed', color='m', head_width=25., head_length=50.)
	#
	# random location:
	# pick a random point near the origin (block position):
	x_r = min(zXs[0]) + (max(zXs[0])-min(zXs[0]))*random.random()
	y_r = min(zXs[1]) + (max(zXs[1])-min(zXs[1]))*random.random()
	dx_r = x_r - Xs[0][0]
	dy_r = y_r - Xs[0][1]
	#
	# plot random point with "dotted" arrow:
	plt.plot(x_r, y_r, 'rs')
	plt.arrow(Xs[0][0], Xs[0][1], dx_r, dy_r, color='r', ls='dotted', head_width=25., head_length=50.)
	#
	# this rotates the vector to its position relative to the strike vector:
	dx_r_prime, dy_r_prime, dz_r_prime = rotate_z([dx_r, dy_r, 0.], -theta_strike)
	plt.arrow(Xs[0][0], Xs[0][1], dx_r_prime, dy_r_prime, color='c', ls='dotted', head_width=25., head_length=50.)
	print "rotated angle: %f (%f deg)" % (math.atan(dy_r_prime/dx_r_prime), 180.*math.atan(dy_r_prime/dx_r_prime)/math.pi)
	print "(%f, %f)" % (dx_r_prime, dy_r_prime)
	print "initial angle: %f (%f deg)" % (math.atan(dy_r/dx_r), 180.*math.atan(dy_r/dx_r)/math.pi)
	print "slip_phi: %f (%f deg)" % (theta_strike, 180.*theta_strike/math.pi)
	print "theta_initial - theta_rotation {should ==}: %f - %f =?= %f (%f)" % (math.atan(dy_r/dx_r), theta_strike, math.atan(dy_r/dx_r)-theta_strike, math.atan(dy_r_prime/dx_r_prime))
	
	print vW, vL, W, L
	
	return Xs
#
def get_okada_location(position_vec=[], origin_vec=[0., 0., 0.], theta_CCW_from_x=0.):
	# calculate relative position for Okada, etc. type calculations. Okada calculations assume strike lies along x axis, dip, rake, etc. will
	# be accounted for in the Okada calculation.
	#
	# position - origin:
	#print "okada_position: ", position_vec
	#print "okada_origin: ", origin_vec
	relative_position = [position_vec[i] - origin_vec[i] for i in xrange(len(position_vec))]	# or we could use numpy.array()
	#
	# rotate relative position vector:
	#rp_prime = rotate_z(relative_position, -theta_strike)
	#
	return rotate_z(relative_position, -theta_CCW_from_x)
#
def test_okada(blockwise_obj=None, block_ids=[], dx=None, dy=None, i_start=0, i_stop=-1, plot_factor=10., fnum=7):
	#
	if blockwise_obj==None: blockwise_obj=numpy.load('dumps/blockwise_290ct.pkl')
	if block_ids==None or block_ids==[]:
		#
		block_ids = list({rw['block_id'] for rw in blockwise_obj.itervalues()})
	#
	# some scripts to test our okada rotations, etc.
	#
	block_ids.sort()
	#
	i_start = (0 or i_start)
	i_stop = (0 or i_stop)
	# blockwise_obj is a dict (or dict-like) object, with keys: BWS[section_id]
	if isinstance(blockwise_obj, str):
		blockwise_obj = numpy.load(blockwise_obj)
	#
	block_L = numpy.linalg.norm(numpy.array([blockwise_obj[0]['m_x_pt4'], blockwise_obj[0]['m_y_pt4'], blockwise_obj[0]['m_z_pt4']]) - numpy.array([blockwise_obj[0]['m_x_pt1'], blockwise_obj[0]['m_y_pt1'], blockwise_obj[0]['m_z_pt1']]))
	block_W = numpy.linalg.norm(numpy.array([blockwise_obj[0]['m_x_pt4'], blockwise_obj[0]['m_y_pt4'], blockwise_obj[0]['m_z_pt4']]) - numpy.array([blockwise_obj[0]['m_x_pt3'], blockwise_obj[0]['m_y_pt3'], blockwise_obj[0]['m_z_pt3']]))
	#
	if dx == None:
		#dx = block_L
		dx=1.0
	if dy==None:
		#dy=block_L
		dy=1.0
	# if fractional values, assume we mean fraction of block size:
	#if dx>0. and dx<1.: dx = dx*block_L
	#if dy>0. and dy<1.: dy = dy*block_L
	dx = dx*block_L
	dy = dy*block_L
	#
	#
	print "block size: %f x %f" % (block_L, block_W)
	print "dx = %f, dy = %f" % (dx, dy)
	#
	# for now, assume we're getting a proper "blockwise_obj"...
	try:
		# do we have a dict or recarray (aka, named columns)? we need to test more thoroughly, but this is a start.
		blockwise_obj[0]['positions'].dtype
		# and if this doesn't fail then...
		#x_index = 'x'
		#y_index = 'y'
		#z_index = 'z'
		#slip_index = 'slip'
		xyz_indices = ['x', 'y', 'z']
	except:
		# then guess...
		#x_index = 1
		#y_index = 2
		#z_index = 3
		#slip_index = 4
		xyz_indices = [1,2,3]
	#
	# try to use this derived dict. object instead (but note, you can use a regular dictionary object in place of meta_dict() ):
	# dx=None, dy=None, i_start=0, i_stop=-1, sections=None, plot_factor=10., f_out=None
	# (noting that the field_plot_prams dict. could alternatively be assigned as individual attributes: ...,dx=dx, dy=dy, ...
	# index this dictionary like {(i_x, i_y):{'dx':dx, 'dy':dy, 'dz':dz}, ...}, and calculate the index based on... dunno,
	# based on magnitude/L_r, but in this case, block-size, slip magnitude? use plot_factor...
	#
	try:
		disp_field = meta_dict({}, field_plot_prams = {'dx':dx, 'dy':dy, 'i_start':i_start, 'i_stop':i_stop, 'plot_factor':plot_factor, 'file_out':f_out})
	except:
		disp_field = {}
	#
	Okada_obj = quakelib.Okada()
	plt.figure(fnum)
	plt.clf()
	ax = plt.gca()
	arrow_factor=100.
	#
	#for block_id, block_data in blockwise_obj.iteritems():
	for j, block_id in enumerate(block_ids):
		block_data = blockwise_obj[block_id]
		#
		# checking the shorter list...
		#if not (sections==None or block_data['section_id'] in sections): continue
		#
		#dip = block_data['dip_rad']
		#strike = block_data['slip_phi']	# or some phase transformation of this...
		#print "starting block_id: ", block_id
		#
		# block positions to calculate total slip, etc.
		pos_0 = numpy.array([block_data['positions'][i_start][i_x] for i_x in xyz_indices])
		pos_0_okada = [pos_0[0], pos_0[1], 0.]	# for use with Okada calcs.
		pos_final = numpy.array([block_data['positions'][i_stop][i_x] for i_x in xyz_indices])
		total_slip_vector = pos_final-pos_0			# note: 1) these simulated positions are compiled in blockwise_slip(), and 2) we cast as numpy.array() here
													# becasue we'll probably want to dot-product this later (otherwise we'd use a light-n-fast list comprehension)
													# we can get relative distance components for Okada by dotting with this vector.
		slip_mag = numpy.linalg.norm(total_slip_vector)	# note: this will be the "unit_slip" in SimElement::calc_displacement_vector
		slip_unit_vector = block_data['slip_vector']	# slip unit-vector. this is different than the slip_unit on the quakelib side, which is probably
														# "slip for this unit(block)".
		#
		# now, plot the source vector:
		this_color = colors_[j%len(colors_)]
		ax.plot(pos_0[0], pos_0[1], 'o', color=this_color)
		ax.arrow(pos_0[0], pos_0[1], 100*total_slip_vector[0], 100*total_slip_vector[1], head_width=.15, head_length=.25, fc=this_color, ec=this_color )
		#
		# calculate the field sites (and put them on a lattice):
		#
		x0 = float(int((pos_0[0] - block_L*plot_factor)/dx))*dx
		y0 = float(int((pos_0[1] - block_L*plot_factor)/dy))*dy
		#
		x_max = pos_0[0] + block_L*plot_factor
		y_max = pos_0[1] + block_L*plot_factor
		#
		c = abs(block_data['m_z_pt2'])
		#
		this_element = quakelib.SimElement()
		this_element.set_rake(block_data['rake_rad'])
		#this_element.set_dip(block_data['dip_rad'])
		#nothing = [this_element.set_vert(i, quakelib.Vec3(block_data['m_x_pt%d' % k], block_data['m_y_pt%d' % k], block_data['m_z_pt%d' % k])) for i,k in enumerate([1,2,4])]
		nothing = [this_element.set_vert(i, block_data['m_x_pt%d' % k], block_data['m_y_pt%d' % k], block_data['m_z_pt%d' % k]) for i,k in enumerate([1,2,4])]
		#
		#	print "vert: ", this_element.vert(i)
		this_element.set_lame_mu(block_data['lame_mu'])
		this_element.set_lame_lambda(block_data['lame_lambda'])
		#
		# for now, set okada vals here; get quakelib working later...
		US = slip_mag*math.cos(block_data['rake_rad'])
		UD = slip_mag*math.sin(block_data['rake_rad'])
		UT= 0.0
		#
		L = numpy.linalg.norm(numpy.array([block_data['m_%s_pt4' % s] for s in ['x', 'y', 'z']]) - numpy.array([block_data['m_%s_pt1' % s] for s in ['x', 'y', 'z']]))
		W = numpy.linalg.norm(numpy.array([block_data['m_%s_pt2' % s] for s in ['x', 'y', 'z']]) - numpy.array([block_data['m_%s_pt1' % s] for s in ['x', 'y', 'z']]))
		c = abs(this_element.max_depth())
		#dip = this_element.dip()
		#
		this_y = y0
		#
		while this_y<y_max:
			this_x = x0
			y_index = int(this_y/dy)
			#
			while this_x<x_max:
				x_index = int(this_x/dx)
				#
				#disp_vector = (numpy.array(this_x, this_y) - pos_0[0:2])		# displacement vector (between a site and the epicenter).
				#
				# get okada displacement for this location:
				#positon_vector = quakelib.Vec3(*rotate_z(disp_vector, -(block_data['slip_phi'])))
				# use (??):
				disp_vector = get_okada_location(position_vec=[this_x, this_y, 0.], origin_vec=pos_0_okada, theta_CCW_from_x=block_data['slip_phi'])
				#print "disp_vector: ", disp_vector
				this_v = quakelib.Vec3(*disp_vector)
				#print this_v
				#return this_element
				#
				okada_disp = Okada_obj.calc_displacement_vector(this_v, c, this_element.dip(), L, W, US, UD, UT, block_data['lame_lambda'], block_data['lame_mu'])
				#
				#okada_disp = this_element.calc_displacement_vector(quakelib.Vec3(*disp_vector), slip_mag, block_data['lame_lambda'], block_data['lame_mu'])
				#
				#okada_disp = Okada_obj.calc_displacement_vector(positon_vector, c, block_data['dip_rad'], block_L, block_W, slip_mag*math.cos(block_data['rake_rad']), slip_mag*math.sin(block_data['rake_rad']), 0.0, block_data['lame_lambda'], block_data['lame_mu'])
				okada_disp = rotate_z(okada_disp, (+block_data['slip_phi']))
				#
				# now, update disp_field values:
				if not disp_field.has_key((x_index, y_index)): disp_field[(x_index, y_index)]={'xyz':[this_x, this_y, 0.], 'd_xyz':numpy.array([0., 0., 0.])}
				disp_field[(x_index, y_index)]['d_xyz'] += numpy.array(okada_disp)
				#
				#ax.plot(this_x, this_y, '.', color=this_color)
				#ax.arrow(this_x, this_y, arrow_factor*okada_disp[0], arrow_factor*okada_disp[1], head_width=.05, head_length=.1, fc=this_color, ec=this_color )
				#
				this_x += dx
			this_y += dy
		#
	plt.figure(fnum+1)
	plt.clf()
	for i, rw in disp_field.iteritems():
		x = rw['xyz']
		dx = rw['d_xyz']
		#
		plt.plot(x[0], x[1], 'b.')
		plt.arrow(x[0], x[1], arrow_factor*dx[0], arrow_factor*dx[1], 'b')
		
	'''
	if f_out!=None and isinstance(f_out, str):
		print "cPickle.dump() to: %s" % f_out
		with open(f_out, 'w') as fout:
			cPickle.dump(disp_field, fout)
		print "pickled..."
	'''
	#
	return disp_field
#
#####
# ripped off from pyvc (temporary...)
#def plot_backslip(sim_file, duration, section_filter=None, field_type='gravity',cutoff=None,padding=0.08,tag=None,fringes=False):
from pyvc import *
from pyvc import vcutils
from pyvc import vcplotutils
from pyvc import vcexceptions
from pyvc import vcanalysis
def get_field_vals(sim_file=default_sim_file, duration=10000, section_filter=None, field_type='gravity',cutoff=None,padding=0.08,tag=None,fringes=False):
	# just get the field values, using existing quakelib machinery.

	#
	output_directory       = 'backslip_only/'
	field_values_directory = '{}field_values/'.format(output_directory)
	
	if not os.path.exists(output_directory):
	    os.makedirs(output_directory)
	    
	if not os.path.exists(field_values_directory):
	    os.makedirs(field_values_directory)
	    
	
	
	# ----------------------------- Initializing --------------------------
	sys.stdout.write('Initializing plot :: ')
	sys.stdout.flush()
	    
	#---------------------------------------------------------------------------
	# Open the data file.
	#---------------------------------------------------------------------------
	with VCSimData() as sim_data:
	    # open the simulation data file
	    sim_data.open_file(sim_file)
	    
	    geometry    = VCGeometry(sim_data)
	    events      = VCEvents(sim_data)
	    
	    # Get global information about the simulations geometry
	    min_lat     = geometry.min_lat
	    max_lat     = geometry.max_lat
	    min_lon     = geometry.min_lon
	    max_lon     = geometry.max_lon
	    base_lat    = geometry.base_lat
	    base_lon    = geometry.base_lon
	    fault_traces= geometry.get_fault_traces()
	
	    # ------------------------------------------
	    # The blocks in slip_rates define the elements that are used for the plotting
	    slip_rates         = geometry.get_slip_rates(section_filter=section_filter,per_year=True)
	    # Set up the elements to evaluate Green's functions
	    ele_getter         = itemgetter(*slip_rates.keys())
	    element_data       = ele_getter(geometry)
	    
	    # Get event information, filter by section if specified
	    event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], section_filter=section_filter)

	    
	
	# -------------------------------------------------------------
	# Instantiate the field and the plotter
	# -------------------------------------------------------------
	if field_type == 'displacement':
	    EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
	    #EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
	    #EFP.calculate_look_angles(geometry[:])
	elif field_type == 'gravity':
	    EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
	    #EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
	
	
	# Apply backslip
	element_slips  = {bid:-1.0*duration*float(slip_rates[bid]) for bid in slip_rates.keys()}

	PRE = '{}{}_'.format(field_values_directory, int(duration)) 

	# Try and load the fields
	field_values_loaded = EF.load_field_values(PRE)
	if field_values_loaded:
	    sys.stdout.write('\nloaded {} years of backslip'.format(int(duration)))
	else:
	    # If they havent been saved then we need to calculate them
	    sys.stdout.write('\nprocessing {} elements :: '.format(len(element_slips.keys())))
	    sys.stdout.flush()
	        
	    EF.calculate_field_values(
	            element_data,
	            element_slips,
	            cutoff=cutoff,
	            save_file_prefix=PRE)

	# Make the plot and save it
	#generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,field_type='gravity')
	
	#generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,field_type=field_type)
	return EF


