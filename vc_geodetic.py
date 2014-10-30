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
		#mean_x = numpy.mean([rw['m_x_pt%d' % j] for j in [1,2,3,4]])
		#mean_y = numpy.mean([rw['m_y_pt%d' % j] for j in [1,2,3,4]])
		#mean_z = numpy.mean([rw['m_z_pt%d' % j] for j in [1,2,3,4]])
		# in a single list comprehension:
		
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
		fault_phi = math.atan(strike_vector[1]/strike_vector[0])		# strike angle (plus phase)
		strike_len = numpy.linalg.norm(strike_vector)
		strike_vector=numpy.array(strike_vector)/strike_len	# direction of strike
		#
		# get thrust vector:
		#dx_th, dy_th, dz_th 
		thrust_vector = [(-rw['m_%s_pt4' % xyz] - rw['m_%s_pt1' % xyz] + rw['m_%s_pt3' % xyz] + rw['m_%s_pt2' % xyz]) for xyz in ('x', 'y', 'z')]
		dx_th, dy_th, dz_th = thrust_vector
		fault_theta = math.acos(dz_th/math.sqrt(dx_th*dx_th + dy_th*dy_th + dz_th*dz_th))
		thrust_len = numpy.linalg.norm(thrust_vector)
		thrust_vector = numpy.array(thrust_vector)/thrust_len
		#
		# actual slip in the strike and thrust directions:
		#
		
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
			block_id=rw['block_id']
			event_number = rw['event_number']
			#
			event_time = events_data[event_number]
			#
			slip = rw['slip']*plot_factor
			#
			phi   = block_info[block_id]['slip_phi']		# more or less strike angle (maybe exactly strike angle); angle on surface about z^.
			#
			#x0=block_info[key]['positions'][-1][1]
			#y0=block_info[key]['positions'][-1][2]
			#z0=block_info[key]['positions'][-1][3]
			x0, y0, z0 = block_info[block_id]['positions'][-1][1:4]
			#
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
	#
	#return numpy.dot(rotation_matrix_x, vector)
	#return rotate_vector_general(vector=vector, axis=[1., 0., 0.], theta=theta)
	return numpy.dot(numpy.array([[1.0, 0., 0.], [0., math.cos(theta), -math.sin(theta)], [0., math.sin(theta), math.cos(theta)]]), vector)
#
def rotate_y(vector=None, theta=0.):
	# rotate about y axis (again, should hard-code these for speed. is there a numpy implementation?)
	#theta = float(theta)%(math.pi*2.0)
	#
	rotation_matrix_y = numpy.array([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-math.sin(theta), 0., math.cos(theta)]])
	#
	return numpy.dot(numpy.array([[math.cos(theta), 0., math.sin(theta)], [0., 1., 0.], [-math.sin(theta), 0., math.cos(theta)]]), vector)
	#return rotate_vector_general(vector=vector, axis=[0., 1., 0.], theta=theta)
#
def rotate_z(vector=None, theta=0.):
	# ... and we should probably pre-define these, for speed...
	#rotation_matrix_z = numpy.array([[math.cos(theta), -math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]])
	#
	return numpy.dot(numpy.array([[math.cos(theta), -math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]]), vector)
	#return rotate_vector_general(vector=vector, axis=[0., 0., 1.], theta=theta)
