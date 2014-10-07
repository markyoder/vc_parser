import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
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
import cStringIO
import sys
import json
import cPickle
import time
import os
#
import imp
import inspect
import multiprocessing as mpp

#pyvc = imp.load_source('pyvc', '../PyVC/pyvc')
#import pyvc as pyvc
#pvca = imp.load_source('pyvc.vcanalysis', '../PyVC/pyvc/vcanalysis.py')
#import pyvc
#import pyvc.vcanalysis as pvca
#

# sections filter for EMC related queries. (should be the set of fault sections used in the socal/EMC simulations). we can use this filter to mine the full
# AllCal data.
emc_section_filter = {'filter': (16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 69, 70, 73, 83, 84, 92, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 123, 124, 125, 126, 149)}
allcal_full_mks = '../ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'

class io_capture_list(list):
	# a little context handler to capture sdtout() 
	def __init__(self, delimiter=None):
		self.delimiter=delimiter
	def __enter__(self):
		self._stdout = sys.stdout
		sys.stdout = self._stringio = cStringIO.StringIO()
		return self
	#
	def __exit__(self, *args):
		#self.extend(self._stringio.getvalue().splitlines())
		self.extend(self._stringio.getvalue().split(self.delimiter))
		sys.stdout = self._stdout
#
class io_capture_strlike(object):
	# a little context handler to capture sdtout() 
	# it will act sort of like a string, but to make it act really like a string
	# i think we have to override all of its functions to use the str_val.
	#
	def __init__(self, str_val=''):
		self.str_val = str_val
	def __str__(self):
		return self.str_val
	def __repr__(self):
		return self.str_val
		
	def __enter__(self):
		self._stdout = sys.stdout
		sys.stdout = self._stringio = cStringIO.StringIO()
		return self
	#
	def __exit__(self, *args):
		self.str_val = self._stringio.getvalue()
		#
		sys.stdout = self._stdout	

#def get_h5_object(h5_in='../ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5', h5_mode='r'):
#	# return an h5 file object.
#	# ... but it looks like we really can't do this because if we leave the file object open, it creates a context problems.
#	# ... and it's just too easy to open the file (same one line of code).
#	return h5py.File(h5_in, h5_mode)

# example: to get column names:
#> f1=get_h5_object()
#> # sweeps talbe:
#> sweeps = fh5['event_sweep_table']
#> sweep_cols = sweeps.id.dtype
# returns:
# dtype([('event_number', '<u4'), ('sweep_num', '<u4'), ('block_id', '<u4'), ('slip', '<f8'), ('area', '<f8'), ('mu', '<f8'), ('shear_init', '<f8'), ('normal_init', '<f8'), ('shear_final', '<f8'), ('normal_final', '<f8')])
#
# then:
#> fields = sweep_cols.fields.keys()

def get_emc_events(sim_file=allcal_full_mks, sortby='event_magnitude', show=None, event_range=None, section_filter=emc_section_filter, magnitude_filter=None,return_data=True, print_data=False, min_mag=0.0):
#def get_emc_events(*args, **kwargs):
	# incidentally, this is how to pass args for this function directly to the next:
	# x = pvca.detailed_sim_info(**locals())
	# but note, of course, that the parameters must be an exact match (in name, expected value, and there cannot be extras... not without additional
	# processing of course.
	#
	return pvca.detailed_sim_info(**locals())
#
def get_event_sections(sim_file=allcal_full_mks, event_number=None):
	'''
	# get sections involved in a particular event. first, get sections, then blocks from that section.
	# this will probably all change for the better. for now, use some available tools from PyVC
	# returns a dictionary object like { <section id (int)>:<section name (str)> }
	# note the use of the stdout listener object io_capture_strlike().
	'''
	#
	# first, use the context handler to capture the event fault sections:
	# use the io_capture() class defined above.
	#with io_capture_list('\n') as sections_list:
	with io_capture_strlike() as sections_list:
		pvca.event_sections(sim_file, event_number)
	#
	sections_list = str(sections_list).strip()
	# now, split by line to get the event id (which of course we know... so we'll actually throw it away):
	sections_list = sections_list.split('\n')[1]
	#return sections_list
	# now, each section is split by a tab, then the section_id, section name are written like '{sec_number}{space(s)}{sec_name}'
	sections = {}
	for sec in sections_list.split('\t'):
		#print "sec: ", sec
		if len(sec)==0: continue
		num_str = ''
		i=0
		while sec[i].isdigit():
			i+=1
		if i==0: continue
		#
		sections[int(sec[0:i])]=sec[i:].strip()
	#
	#return sections_list
	return sections
#

def get_event_block_details(sim_file=allcal_full_mks, event_number=None, block_table_name='block_info_table'):
	'''
	# get blocks involved in a particular event. first, get sections, then blocks from that section.
	# this will probably all change for the better. for now, use some available tools from PyVC
	#
	# ... but have a closer look at this. i think it might be a bit redundant. event_sweep_table includes an event_number, so we should
	# be able to forego collecting all the sections, etc. unless we require detailed information about the blocks in the event... which this will give us.
	'''
	#
	#
	#event_sections = get_event_sections(**locals())
	event_sections = get_event_sections(sim_file=sim_file, event_number=event_number)
	
	section_ids = event_sections.keys()
	# now, sections is a dictionary of {sec_num: sec_name, ...} for event number event_number.
	# get these blocks from sim_file.
	# block data are in the block_info_table, uhh, table. i think just crack open the table and pull the data...
	#
	with h5py.File(sim_file, 'r') as vc_data:	
		block_info = vc_data[block_table_name]	# we should exception trap this of course...
		block_cols = block_info.id.dtype.fields	# %.dtype is a types container of some sort. %.dtype.fields is a "dictproxy" object... which acts (i think)
												# like a dictionary, so for example we can use .keys() to get just the col name.
												# entries look like {'block_id': (dtype('uint32'), 0), ...}
												# NOTE (IMPORTANT): if we fetch the list of column names, they may not appear in the same order
												# as the data in each row. columns should be specified by name, not index (without confirming position that is)...
												# though they do appear to be in the same order as in the hdf5 viewer tool... and note that %.id.dtype does
												# return the columns in the proper sequence.
		#
		#event_blocks = [x for x in block_info if x['section_id'] in section_ids]
		event_blocks = fetch_data_mpp(n_cpus=None, src_data=block_info, col_name='section_id', matching_vals=section_ids)
	#
	return [block_cols] + event_blocks
#
def get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, is_sorted=False):
	#
	# is_sorted: some sort based indexing in the find_it() (or whatever) function. this shortens the list over which
	# the data are searched, but at the expense of the list comprehension... so it actually runs slower by about 30% in preliminary tests.
	#
	
	# ... and what's weird is that this seems to run WAY faster using multiprocessing (non-linearly faster) -- 1 cpu running for a minute or so
	# and maybe a second or two using 4 processors... and in fact, even running 1 processor but using the multiprocessing structure is much,
	# much faster.
	#
	# allow for multi-processing:
	if n_cpus == None: n_cpus = mpp.cpu_count()
	#
	with h5py.File(sim_file,'r') as vc_data:
		#
		#
		view_ebs = vc_data['events_by_section']		# "events by section" view
		event_ids_tbl = view_ebs['section_%d' % section_id]
		#
		event_ids_ary = numpy.array(numpy.zeros(event_ids_tbl.len()), dtype='int64')	# array needs to be initialized with data type for hdf5.read_direct()
		#event_ids_ary.shape(numpy.zeros(event_ids_tbl.len(), len(event_ids_tbl[0])))
		event_ids_tbl.read_direct(event_ids_ary)
		#
		# now, fetch event data:
		# this is a choice, computationally. we can do it in a list comprehension [], or we can index the source list (event_id values
		# are sorted, so once we find one event in the master list, we can stop looking for it in the event_ids_ary list.
		# list comprehension is simpler...
		src_data = vc_data['event_table']
		col_names = cols_from_h5_dict(src_data.id.dtype.fields)
		col_name = 'event_number'
		#if n_cpus in (1, None):
		if n_cpus == None:
			# this should never be true any longer. see note above. strangely, using the mpp setup seems to make this run much, much
			# faster, even if only using the single CPU.
			#
			#section_events = [x for x in vc_data[src_data] if x[col_name] in event_ids_ary[:]]
			section_events = find_in(tbl_in=src_data, col_name=col_name, in_list=event_ids_ary, pipe_out=None) 
		else:
			my_pipes = []
			section_events = []
			my_processes = []
			sublist_len = min(len(src_data), (len(src_data)/n_cpus)+1)		# a little trick to ensure we get the full list and don't overrun...
			#
			for i in xrange(n_cpus):
				child_pipe, parent_pipe = mpp.Pipe()
				my_pipes+= [parent_pipe]
				#my_processes+=[mpp.Process(target=find_in, args=(src_data[i*sublist_len:(i+1)*sublist_len], col_name, event_ids_ary, child_pipe))]
				my_processes+=[mpp.Process(target=find_in, kwargs={'tbl_in':src_data[i*sublist_len:(i+1)*sublist_len], 'col_name':col_name, 'in_list':event_ids_ary, 'pipe_out':child_pipe, 'is_sorted':is_sorted})]
				my_processes[-1].start()
			#
			print "%d/%d processes started." % (len(my_pipes), len(my_processes))
			#
			# processes are now started and joined, and we'll wait until they are all finished (each join() will hold the calling
			# process until it's finished)
			for i, pipe in enumerate(my_pipes):
				in_my_pipe = pipe.recv()
				#print "adding %d new row..." % len(in_my_pipe)
				section_events += in_my_pipe
				my_pipes[i].close()
			#
			# and join, but note that the join() must follow pipe.recv(), or it hangs forever and ever because the sending
			# pipe never closes...
			for i, proc in enumerate(my_processes): 
				if my_processes[i].is_alive(): my_processes[i].join()					
	#
	#return event_ids_ary
	# columns we might care about:
	return [col_names] + section_events
#
def get_blocks_info_dict(sim_file=allcal_full_mks, block_ids=None):
	# if block_ids==None, get all of them...
	# and we could parallelize this, or we can use it as a single processor implementation that can be run in parallel...
	#
	with h5py.File(sim_file, 'r') as vc_data:
		blocks_table = vc_data['block_info_table']
		block_cols = cols_from_h5_dict(blocks_table.id.dtype.fields)
		#
		blockses = {}
		for block in blocks_table:
			if block_ids!=None and block['block_id'] not in block_ids: continue
			#
			block_dict = {}
			#
			[block_dict.__setitem__(key, block[key]) for key in block_cols]
			blockses[block['block_id']] = block_dict
		#
	#
	return blockses
#
# a suite of functions to get total CFF for an event by summing the CFF from all the blocks in the event.
# the variety of functions are to facilitate different approaches to multiprocessing...
#
def calc_final_CFF_from_block_data(block_data):
	return block_data['shear_final'] - block_data['mu']*block_data['normal_final']
#
def calc_CFF_from_block_data(block_data, initial_stress=True):
	# a row of block data (assume dict. or hdf5 format) (aka, data from 1 block).:
	if initial_stress==False:
		return block_data['shear_final'] - block_data['mu']*block_data['normal_final']
	else:
		return block_data['shear_init'] - block_data['mu']*block_data['normal_init']

def calc_CFFs_from_blocks(sweep_blocks, pipe=None):
	# a single-processor (spp) function to calculate total initial and final CFF from a set of "sweep" blocks:
	# note this is meant to be single-process (mpp.Process() ) friendly.
	#
	total_cff_init  = 0.
	total_cff_final = 0.
	#
	n_blocks = float(len(sweep_blocks))
	#
	for blk in sweep_blocks:
		total_cff_init  += blk['shear_init']  - blk['mu']*blk['normal_init']
		total_cff_final += blk['shear_final'] - blk['mu']*blk['normal_final']
	#
	if pipe!=None:
		try:
			pipe.send({'cff_init':total_cff_init, 'cff_final':total_cff_final})
		except:
			pipe=None
	if pipe==None: 
		return {'cff_init':total_cff_init, 'cff_final':total_cff_final, 'mean_cff_init':total_cff_init/n_blocks, 'mean_cff_final':total_cff_init/n_blocks}
	
def calc_CFF(shear_stress=0., normal_stress=0., mu=.8):
	return shear_stress - normal_stress*mu
#
def get_event_CFF(sweep_blocks, n_cpus=None, chunksize=100):
	# @sweep_blocks: blocks from event_sweep_table (for a given event(s))
	#
	# CFF = shear_stress - mu*normal_stress
	if n_cpus==None: n_cpus=mpp.cpu_count()
	#
	my_pool = mpp.Pool(processes=n_cpus)
	CFFs = my_pool.map_async(calc_CFF_from_block_data, sweep_blocks, chunksize=chunksize)	# imap returns an ordered iterable. map_async()
																		# returns a "results object"
																		# which is better? dunno. 
	my_pool.close()	# this is very very important. note a Pool() can be executed recursively, so this tells the Pool() object
					# that no more tasks are coming its way. otherwise, it remains open an eats up a bunch of memory and makes
					# a huge mess in the mem. stack.
					# don't know if the join() is necessary...
	my_pool.join()												
	#
	return sum(CFFs.get())
#
def get_final_event_CFF(sweep_blocks, n_cpus=None, chunksize=100):
	# aka, calc the "final", not the "initial" cff for an event (given the set of blocks for
	# that event).
	# @sweep_blocks: blocks from event_sweep_table (for a given event(s))
	# (see get_event_CFF() for additional notes)
	#
	# CFF = shear_stress - mu*normal_stress
	if n_cpus==None: n_cpus=mpp.cpu_count()
	#
	my_pool = mpp.Pool(processes=n_cpus)
	CFFs = my_pool.map_async(calc_final_CFF_from_block_data, sweep_blocks, chunksize=chunksize)	
	my_pool.close()	
	my_pool.join()												
	#
	return sum(CFFs.get())
#
def cff_dict_npy(dict_in):
	# convert an older style cff dict to a structured numpy array.
	# this is not tested, but it should work with minor tweaks at most.
	#
	dtype_names = []
	#dtype_types = []
	#
	lst_data = []
	dtype_names = [key for key in dict_in[0].keys()]
	#
	# since we're doing this dynamically, data will have to be accessed dynamically, but that shouldn't be a problem...
	#
	for rw in dict_in:
		lst_data += [[]]
		for key in keys:
			lst_data[-1]+=[rw[key]]
		#print lst_data[-1]
		#lst_rw = numpy.array(lst_rw, dtype=dtypes)
		
		#
	#
	
	# outputs = numpy.core.records.fromarrays(zip(*outputs), names=output_names, formats = [type(x).__name__ for x in outputs[0]])
	
	ary_data = numpy.core.records.fromarrays(zip(*lst_data), names=dtype_names, formats=[type(x).__name__ for x in lst_data[0]])
	#
	return lst_data
		
	#
#
def forecast_metric_1(ary_in='data/VC_CFF_timeseries_section_16.npy', m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False):
	'''
	#
	# forecast based on seismic acceleration. specifically, for a set of data (a fault-section time-series),
	# find the slopes of the inter-occurrence intervals. slope<0 implies acceleration, implies hazard.
	# the evaluation metric will be basically a ROC type:
	# 1) total time "at alert" serves as "false alarm"
	# 2) fraction of "predicted" events. this will vary. for starters:
	#     - slope at time t_i of event or time t_i-1 of the previous event is below b_0 (typically, b_0=0.0).
	#     - the "_i-1" condition is fair because 1) we don't know the slope is changin until the even happens
	#       and 2) in a real scenario, we'd have aftershocks, which woluld turn the post-seismic slope positive.
	'''
	#
	CFF = numpy.load(ary_in)
	#
	recurrence_data = mean_recurrence(ary_in=CFF, m0=m0)
	nyquist_len = max(int(nyquist_factor*recurrence_data['mean_dN_fault']), 2)
	nyquist_time = nyquist_factor*recurrence_data['mean_dT']
	#
	trend_data = get_trend_analysis(ary_in=CFF, nyquist_len = nyquist_len, nyquist_time=nyquist_time)
	#trend_data_dict = {trend_data['event_year']:x for x in trend_data}
	CFF_dict = {x['event_year']:x for x in CFF}
	print "trend lengths: ", len(trend_data), len(CFF), nyquist_len
	max_n = len(trend_data)
	#
	# trend_data includes columns: ['event_number', 'event_year', 'lin_fit_a', 'lin_fit_b', 'rb_ratio', 'interval_rate_ratios']
	#
	# first, just get the total time under alert:
	alert_time = 0.0
	alert_segments = [[]]		# a collection of lists...
	running_b_sequence = []
	#
	for i, rw in enumerate(trend_data):
		# when we find b<b_0, we issue an alert until the next earthquake -- unless this one was 'big',
		# in this case m>=7.
		#
		if i>=(max_n-1): break
		#
		# note the comment evidence of variations on this metric, primarily involving some sort of mean-slope averaging.
		# a more exhaustive, automated MC approach, analysis is necessary to be certain, but preliminary analysis suggests that
		# we don't gain anything from the averaging... in fact, we tend to loose measurable, at least on fault 16, where most
		# of the prelim examination was done.
		
		#mean_b = numpy.mean(trend_data['lin_fit_b'][max(0, i-nyquist_len) : i+1])
		#
		this_b = rw['lin_fit_b']
		#
		#this_b = mean_b
		#running_b_sequence = (running_b_sequence + [this_b])[max(0, len(running_b_sequence)-nyquist_len):]
		#mean_b = numpy.mean(running_b_sequence)
		#this_b = mean_b
		#
		if this_b >= b_0:
			# null case
			# if we've been collecting "alert" events, stop. if not, just troll along...
			if len(alert_segments[-1])>0:
				#
				alert_segments+= [[]]
		
		else:
			# accelerating case:
			#this_mag = CFF[-max_n:][i+1]['event_magnitude']
			this_mag = CFF_dict[rw['event_year']]['event_magnitude']
			#
			if len(alert_segments[-1])==0:
				alert_segments[-1]+=[[trend_data[i]['event_year'], trend_data[i]['lin_fit_b']]]
			#
			if this_mag<m0:
				# add the next event as the alert (aka, issue an alert until we have new data).
				#print "len alert_seg: %d" % len(alert_segments[-1])
				#
				# generalize language a bit:
				alert_segments[-1]+=[[trend_data[i+1]['event_year'], this_b]]
				pass
				#
			if this_mag>=m0:
				# this is "the" earthquake. add this entry (it's probably already there) from the previous entry.
				#alert_segments[-1]+=[[trend_data[i]['event_year'], trend_data[i]['lin_fit_b']]]
				alert_segments[-1]+=[[trend_data[i]['event_year'], this_b]]
				#running_b_sequence=[]
			#
		#
	#
	while len(alert_segments[-1])==0: alert_segments.pop()
	#
	# now, calculate total alert time:
	total_alert_time = 0.0
	total_total_time = CFF[-1]['event_year'] - CFF[0]['event_year']
	#
	for alert_segment in alert_segments:
		total_alert_time += (alert_segment[-1][0] - alert_segment[0][0])
	#
	# and prediction success:
	n_predicted = 0
	n_missed = 0
	#
	# was an alert issued at the time of an m>m0 event?
	# we should clean this up, but for now, make and use dictionaries...
	#alert_dict = {x[0]:x[1] for x in 
	alert_dict = {}
	for seg in alert_segments:
		#
		for rw in seg:
			alert_dict[rw[0]] = rw[1]
	#
	for i, rw in enumerate(CFF):
		if rw['event_magnitude']<m0: continue
		#
		# this is a really dumb way to do this, so fix it:
		#
		# otherwise, its a big one. we need two conditions:
		# 1) the alert was active during the event.
		# 2) the alert had already been issued.
		CFF_index = rw['event_number']
		prev_year = CFF[i-1]['event_year']
		#if alert_dict.has_key(rw['event_year']):
		#if alert_dict.has_key(rw['event_year']) and alert_dict.has_key(prev_year):
		if alert_dict.has_key(rw['event_year']):						
			n_predicted += 1
		else:
			n_missed += 1
	
	#
	#
	if do_plot:
		# diagnostic plots of forecast metric:
		plt.figure(0)
		plt.clf()
		#
		ax_metric = plt.gca()
		ax_mags = ax_metric.twinx()
		min_mag = min(CFF['event_magnitude'])
		#
		for segment in alert_segments:
			X,Y = zip(*segment)
			ax_metric.set_yscale('linear')
			ax_metric.fill_between(X,[-y+min_mag for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [y<0. for y in Y] )
			#
			ax_mags.fill_between(X, [min_mag for x in X], [m0 for x in X], zorder=5, alpha=.2, color='m')
	
	
		ax_mags.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='b', alpha=.7)
		X_big_mags, Y_big_mags = zip(*[[x['event_year'], x['event_magnitude']] for x in CFF if x['event_magnitude']>m0])
		ax_mags.vlines(X_big_mags, [min_mag for x in Y_big_mags], Y_big_mags, color='r', alpha=.9, lw=2.5)
	#
	#
	print "preliminary report:"
	print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, total_alert_time/total_total_time)
	print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, float(n_predicted)/(float(n_predicted)+n_missed))
	print "total: %f " % (float(n_predicted)/(float(n_predicted)+n_missed) - total_alert_time/total_total_time)
		
	
	#
	#return alert_segments
	# , 'alert_segments':alert_segments
	return {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'alert_segments':alert_segments, 'ary_in_name':ary_in, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
#
def plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=0.0, nyquist_factor=.5, do_spp=False, do_plot=False, do_clf=True, n_cpus=None):
	'''
	# scatter plot of hit_ratio vs alert_time_ratio for as many data as we throw at it.
	'''		
	# wrapper to convert a bunch of normal arrays or maybe lists to numpy structured arrays (numpy.recarray).
	G=glob.glob(file_profile)
	#
	X, Y = [], []
	#
	if do_spp:
		# SPP version (if we want to force this for some reason):
		for g in G:
			fc_data = forecast_metric_1(ary_in=g, m0=m0, b_0 = b_0, nyquist_factor=nyquist_factor)
			#
			Y+=[float(fc_data['n_predicted'])/(fc_data['n_predicted'] + float(fc_data['n_missed']))]
			X+=[fc_data['total_alert_time']/fc_data['total_time']]
	#
	else:
		# use MPP:
		if n_cpus==None: n_cpus = mpp.cpu_count()
		pool = mpp.Pool(n_cpus)
		if (m0, b_0, nyquist_factor)==forecast_metric_1.__defaults__[1:4]:
			print "defaults. use map_async()"
			# i'm guessing this is faster...
			result_set = pool.map_async(forecast_metric_1, G)
			pool.close()
			pool.join()
			#
			#print result_set
			resultses = result_set.get()	# will be a list of dictionary objects
			
		else:
			# add/"apply" each file to the pool.
			print "not default. use apply_async()"
			result_set_list = [pool.apply_async(forecast_metric_1, args = (g, m0, b_0, nyquist_factor)) for g in G]
			pool.close()
			pool.join()
			#
			#print "these prams: ", b_0, nyquist_factor
			#return result_set_list
			resultses = [x.get() for x in result_set_list]	# each entry in result_set_list is a dict; now we have a list of dicts.
			
		#
		# still MPP...
		#
		XY = [[x['total_alert_time']/x['total_time'] ,float(x['n_predicted'])/(x['n_predicted'] + float(x['n_missed']))] for x in resultses]
		X,Y = zip(*XY)
	# end MPP. now, either by SPP or MPP, we have X,Y
	#
	#
	#return fc_datas
	#
	if do_plot:
		#
		plt.figure(0)
		if do_clf: plt.clf()
		plt.plot(X, Y, '.')
		plt.plot([0., 1.], [0., 1.], '-')
		plt.ylabel('n_predicted/n_total')
		plt.xlabel('percent alert time')
	#
	if resultses[0].has_key('alert_segments'): [x.pop('alert_segments') for x in resultses]
	return resultses
#
def plot_aggregate_metric(scores_in, n_top=None):
	# make a pretty (3d) plot of the aggregate type optimizer solutions.
	if isinstance(scores_in, str):
		scores_in = numpy.load(scores_in)
	#scores_in.sort(order=('b_0', 'nyquist_factor'))
	#
	f2 = plt.figure(3)
	f2.clf()
	ax3d = f2.add_subplot(111, projection='3d')
	ax3d.plot(scores_in['b_0'], scores_in['nyquist_factor'], scores_in['score'], '.')
	ax3d.set_xlabel('threshold slope, $b_0$')
	ax3d.set_ylabel('nyquist_factor')
	ax3d.set_zlabel('ROC metric, (Percent Predicted) - (False alarm rate)')
	#
	if n_top==None:
		n_top = min(100, int(len(scores_in)/100))
	#
	# and plot the top n_top values separately:
	scores_in.sort(order='score')
	ax3d.plot(scores_in['b_0'][-n_top:], scores_in['nyquist_factor'][-n_top:], scores_in['score'][-n_top:], 'yo')
	#
	#scores_in.sort(order=('nyquist_factor', 'b_0'))
	#ax3d.plot(scores_in['b_0'], scores_in['nyquist_factor'], scores_in['score'], '.')
	#
	# best score?
	best_row = scores_in[scores_in['score'].tolist().index(max(scores_in['score']))]
	ax3d.plot(scores_in['b_0'][-1:], scores_in['nyquist_factor'][-1:], scores_in['score'][-1:], 'r*', ms=15, label='best: $b_0 = %.3f, nf=%.3f, score=%.3f \pm %.3f$' % (best_row['b_0'], best_row['nyquist_factor'], best_row['score'], best_row['stdev']) )
	#
	plt.legend(loc=0, numpoints=1)
	#
	print scores_in.dtype.names
	print best_row
	return scores_in
#
def optimize_metric_full_aggregate(b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01, nits=None, do_plot=False, n_cpus=None, data_in=None):
	# run a whole bunch of metric_1 and get the best nyquist_factor, b_0 combination.
	# this runs a fully composite optimization, which produces some not super believable... or optimal
	# results. it seems that we would do better to optimize each fault independently.
	# data_in: start with, and append, an existing data file.	

	#
	R_b   = random.Random()
	R_nyq = random.Random()
	delta_b = b_max-b_min
	delta_nyq = nyquist_max - nyquist_min
	# and let's do this right and randomly sample...
	if nits==None: nits = 1 + int(abs((delta_b/d_b)*((delta_nyq)/d_nyquist)))	# safe-guard for n->0
	#
	total_scores = []	# cumulative total score, like [[b_0, nyquist_factor, mean_score, score_stdev]]
	if data_in!=None: 
		if isinstance(data_in, str):
			total_scores = numpy.load(data_in).tolist()
		else:
			total_scores = data_in
	#
	for i in xrange(nits):
		this_b   = b_min       + delta_b*R_b.random()
		this_nyq = nyquist_min + delta_nyq*R_nyq.random()
		#
		# quick mod to remove the apparent "fobidden" zone:
		while this_b < (-.3143 + .2286*this_nyq):
			# we've determined that we don't get valid results from this domain.
			this_b   = b_min       + delta_b*R_b.random()
			this_nyq = nyquist_min + delta_nyq*R_nyq.random()
		print "************\n*************\n***************\n*************\n"
		
		try:
			datas = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=this_b, nyquist_factor=this_nyq, do_spp=False, do_plot=False, n_cpus=n_cpus)	# note: this will fully multiprocess.
		except:
			print "ERROR!!! datas would not assimilate. probably bogus prams: b=%f, nq_fact=%f" % (this_b, this_nyq)
			continue
		#
		# now, aggregate for this pram-set.
		# datas is a list of dictionary objects like: {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'ary_in_name':ary_in, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
		#
		#score_row = [x['b'], x['nyquist_factor'], x['n_predicted']/(x['n_predicted']+x['n_missed']), x['total_alert_time']/x['total_time'], (x['n_predicted']/(x['n_predicted']+x['n_missed']) - x['total_alert_time']/x['total_time']) for x in datas]
		scores = [float(x['n_predicted'])/(float(x['n_predicted'])+float(x['n_missed'])) - x['total_alert_time']/x['total_time'] for x in datas]
		mean_score = numpy.mean(scores)
		score_stdev = numpy.std(scores)
		#
		total_scores += [[this_b, this_nyq, mean_score, score_stdev]]
		#
		if i%100==0:
			# intermediate dump...
			dump_object = numpy.core.records.fromarrays(zip(*total_scores), names=['b_0', 'nyquist_factor', 'score', 'stdev'], formats = [type(x).__name__ for x in total_scores[0]])
			try:
				dump_object.dump('dumps/aggregate_optimize_n_%d.npy' % nits)
			except:
				print "failed to dump array..."
		#
	total_scores.sort(key = lambda x: x[2])
	#
	total_scores = numpy.core.records.fromarrays(zip(*total_scores), names=['b_0', 'nyquist_factor', 'score', 'stdev'], formats = [type(x).__name__ for x in total_scores[0]])
	#
	try:
		total_scores.dump('dumps/aggregate_optimize_n_%d.npy' % nits)
	except:
		print "failed to dump array..."
	#
	A=total_scores	# for shorthand...
	best_fit_row = A[A['score'].tolist().index(max(A['score']))]
	print best_fit_row
	if do_plot:
		plt.ion()
		f=plt.figure()
		ax3d = f.add_subplot(111, projection='3d')
		ax3d.plot(A['b_0'], A['nyquist_factor'], A['score'], '.')
	#
	return total_scores

#
def optimize_metric_faultwise(b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01, nits=None, dump_file='dumps/optimize_faultwise'):
	# run a whole bunch of metric_1 and get the best nyquist_factor, b_0 combination.
	# this runs a fully composite optimization, which produces some not super believable... or optimal
	# results. it seems that we would do better to optimize each fault independently.
	#
	R_b   = random.Random()
	R_nyq = random.Random()
	delta_b = b_max-b_min
	delta_nyq = nyquist_max - nyquist_min
	# and let's do this right and randomly sample...
	if nits==None: nits = 1 + int(abs((delta_b/d_b)*((delta_nyq)/d_nyquist)))	# safe-guard for n->0
	#
	fault_scores={}
	#
	for i in xrange(nits):
		this_b   = b_min       + delta_b*R_b.random()
		this_nyq = nyquist_min + delta_nyq*R_nyq.random()
		print "************\n*************\n***************\n*************\n"
		
		try:
			datas = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=this_b, nyquist_factor=this_nyq, do_spp=False, do_plot=False, n_cpus=None)	# note: this will fully multiprocess.
		except:
			print "ERROR!!! datas would not assimilate. probably bogus prams: b=%f, nq_fact=%f" % (this_b, this_nyq)
			continue
		#
		# now, aggregate for this pram-set.
		# datas is a list of dictionary objects like: {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'ary_in_name':ary_in, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
		for rw in datas:
			#
			if fault_scores.has_key(rw['ary_in_name'])==False:
				fault_scores[rw['ary_in_name']] = {}
			this_score = float(rw['n_predicted'])/(float(rw['n_predicted'])+rw['n_missed']) - rw['total_alert_time']/rw['total_time']
			fault_scores[rw['ary_in_name']][this_score] = rw	# index each row by the score for fast sorting.
	#
	# now, get the max score set for each row:
	#return fault_scores
	
	scores_out = []
	for key in fault_scores.keys():
		rw = fault_scores[key]
		i0 = key.index('section_') + len('section_')
		fault_number = int(key[i0:key.index('.', i0)])
		#
		best_score = max(rw.keys())
		best_set = rw[best_score]
		#scores_out += [[key, best_set['b'], best_set['nyquist_factor'], best_score]]
		
		# controlled:
		scores_out += [[fault_number, best_set['b'], best_set['nyquist_factor'], best_set['n_predicted'], best_set['n_missed'], best_set['total_alert_time'], best_set['total_time'], best_score]]
		# automated:
		#scores_out += [[fault_number] + [best_set[key] for key in best_set] + [best_score]]
	#
	#
	
	#print scores_out[0]
	#print len(scores_out), len(scores_out[0])
	scores_out = numpy.core.records.fromarrays(zip(*scores_out), names=['fault_id', 'b_0', 'nyquist_factor', 'n_predicted', 'n_missed', 'total_alert_time', 'total_time', 'score'], formats = [type(x).__name__ for x in scores_out[0]])
	#
	scores_out.dump('%s_best_scores.npy' % dump_file)
	numpy.array(fault_scores).dump('%s_fault_scores.npy' % dump_file)
	#
	return scores_out

def plot_best_opt_prams(scores_in):
	'''
	# plots from faultwise optimization
	#
	'''
	if isinstance(scores_in, str):
		scores_in = numpy.load(scores_in)
	#
	b_col = scores_in['b_0']
	nq_col = scores_in['nyquist_factor']
	#
	# information gain:
	plt.figure(0)
	X = [scores_in['total_alert_time'][i]/t for i, t in enumerate(scores_in['total_time'])]
	Y = [float(N)/(float(N)+scores_in['n_missed'][i]) for i, N in enumerate(scores_in['n_predicted'])]
	plt.clf()
	plt.plot(X,Y, '.')
	plt.plot([0., 1.], [0.,1.], '-')
	plt.xlabel('percent alert time')
	plt.ylabel('percent predicted')
	#
	# best-fit parameters:
	plt.figure(1)
	plt.clf()
	ax=plt.gca()
	X=scores_in['b_0'].copy()
	lst_scores=scores_in['score'].copy()
	XY = zip(*[X, lst_scores])
	XY.sort(key=lambda x:x[0])
	X,Y = zip(*XY)
	#ax.plot(scores_in['b_0'], scores_in['score'], 'b.-')
	ax.plot(X,Y, 'b.-')
	ax.set_xlabel('b_0')
	ax.set_ylabel('score')
	#
	#ax2 = ax.twiny()
	plt.figure(2)
	plt.clf()
	ax2=plt.gca()
	XY = zip(*[scores_in['nyquist_factor'].copy(), scores_in['score'].copy()])
	XY.sort(key=lambda x: x[0])
	X,Y = zip(*XY)
	ax2.plot(X, Y, 'g.-')
	#ax2.plot(scores_in['nyquist_factor'], scores_in['score'], 'g.-')
	ax2.set_xlabel('nyquist factor')
	#
	f2 = plt.figure(3)
	f2.clf()
	ax3d = f2.add_subplot(111, projection='3d')
	ax3d.plot(scores_in['b_0'], scores_in['nyquist_factor'], scores_in['score'], '.')
	#
	mean_b = numpy.mean(scores_in['b_0'])
	std_b  = numpy.std(scores_in['score'])
	#
	mean_nf = numpy.mean(scores_in['nyquist_factor'])
	std_nf = numpy.std(scores_in['nyquist_factor'])
	#
	print "mean_b: %f +/- %f" % (mean_b, std_b)
	print "mean_nf:    %f +/- %f" % (mean_nf, std_nf)

#
def bunch_o_metric_1(nyquist_factor=0.425556):
	# a bunch of metric 1 data...
	my_b = .1
	all_the_data = {}
	while my_b>-.1:
		datas = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=my_b, nyquist_factor=nyquist_factor, do_spp=False, do_plot=False)
		# datas will be a list of dictionary objects. we want to construct independent "time" series for each fault (dict object):
		for data in datas:
			key = data['ary_in_name']
			if all_the_data.has_key(key) == False: all_the_data[key] = []
			#if all_the_data.has_key(key) == False: all_the_data[key] = {'prams':[], 'data':[]}
			#
			all_the_data[key]+=[[data['total_alert_time']/data['total_time'] ,float(data['n_predicted'])/(data['n_predicted'] + float(data['n_missed'])), data['b'], data['m0'], data['nyquist_factor']]]
			
		my_b-=.02
	#
	# wrap these up into recarrays:
	for key in all_the_data.keys():
		all_the_data[key] = numpy.core.records.fromarrays(zip(*all_the_data[key]), names=['false_alarm', 'hit_rate', 'b', 'm0', 'nyquist_factor'], formats = [type(x).__name__ for x in all_the_data[key][0]])
		#
	return all_the_data

def plot_bunch_o_metric_1(all_the_data=None):
	if all_the_data==None: all_the_data = bunch_o_metric()
	#
	plt.figure(0)
	plt.clf()
	plt.xlabel('percent alert time (false alarm)')
	plt.ylabel('percent predicted')
	plt.plot([0., 1.], [0.,1.], '--')
	#
	plt.figure(1)
	plt.clf()
	plt.xlabel('critical slobe $b_0$')
	plt.ylabel('total score: hit_rate - false_alarm_rate')
	plt.ion()
	#
	best_bs = []
	for key in all_the_data.keys():
		#X,Y = zip(*all_the_data[key])[0:2]
		X,Y = all_the_data[key]['false_alarm'], all_the_data[key]['hit_rate']
		bs = all_the_data[key]['b']
		total_scores = [x['hit_rate'] - x['false_alarm'] for x in all_the_data[key]]
		plt.figure(0)
		plt.plot(X,Y, '.-')
		plt.figure(1)
		plt.plot(bs, total_scores, '.-')
		max_score = max(total_scores)
		max_index = total_scores.index(max_score)
		best_bs+=[bs[max_index]]
		plt.plot([bs[max_index]], [max_score], '*')
		#
		print "max predictor: ", X[max_index], Y[max_index], bs[max_index], max_score
	mean_best_b = numpy.mean(best_bs)
	std_best_b = numpy.std(best_bs)
	print "best b: mean: %f, stdev: %f" % (mean_best_b, std_best_b)
	#
	return all_the_data
#
def bunch_o_metric_nyquist(b_0=0.009333, d_nyq = .05):
	# a bunch of metric 1 data...
	my_nyq = .25
	all_the_data = {}
	while my_nyq < .75:
		datas = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=b_0, nyquist_factor=my_nyq, do_spp=False, do_plot=False)
		# datas will be a list of dictionary objects. we want to construct independent "time" series for each fault (dict object):
		for data in datas:
			key = data['ary_in_name']
			if all_the_data.has_key(key) == False: all_the_data[key] = []
			#if all_the_data.has_key(key) == False: all_the_data[key] = {'prams':[], 'data':[]}
			#
			all_the_data[key]+=[[data['total_alert_time']/data['total_time'] ,float(data['n_predicted'])/(float(data['n_predicted']) + float(data['n_missed'])), data['b'], data['m0'], data['nyquist_factor']]]
			
		my_nyq+=d_nyq
	#
	# wrap these up into recarrays:
	for key in all_the_data.keys():
		all_the_data[key] = numpy.core.records.fromarrays(zip(*all_the_data[key]), names=['false_alarm', 'hit_rate', 'b', 'm0', 'nyquist_factor'], formats = [type(x).__name__ for x in all_the_data[key][0]])
		#
	return all_the_data

def plot_bunch_o_metric_nyquist(all_the_data=None):
	if all_the_data==None: all_the_data = bunch_o_metric_nyquist()
	#
	plt.figure(0)
	plt.clf()
	plt.xlabel('percent alert time (false alarm)')
	plt.ylabel('percent predicted')
	plt.plot([0., 1.], [0.,1.], '--')
	#
	plt.figure(1)
	plt.clf()
	plt.xlabel('nyquist factor')
	plt.ylabel('total score: hit_rate - false_alarm_rate')
	plt.ion()
	#
	best_nyquists = []
	for key in all_the_data.keys():
		#X,Y = zip(*all_the_data[key])[0:2]
		X,Y = all_the_data[key]['false_alarm'], all_the_data[key]['hit_rate']
		nyquists = all_the_data[key]['nyquist_factor']
		total_scores = [x['hit_rate'] - x['false_alarm'] for x in all_the_data[key]]
		plt.figure(0)
		plt.plot(X,Y, '.-')
		plt.figure(1)
		plt.plot(nyquists, total_scores, '.-')
		max_score = max(total_scores)
		max_index = total_scores.index(max_score)
		best_nyquists+=[nyquists[max_index]]
		plt.plot([nyquists[max_index]], [max_score], '*', ms=10)
		#
		print "max predictor: ", X[max_index], Y[max_index], nyquists[max_index], max_score
	mean_best_nyquists = numpy.mean(best_nyquists)
	std_best_nyquists = numpy.std(best_nyquists)
	print "best nyquist: mean: %f, stdev: %f" % (mean_best_nyquists, std_best_nyquists)
	#
	return all_the_data
#
def mean_recurrenceses(section_ids=[], m0=7.0, file_path_format='data/VC_CFF_timeseries_section_%d.npy'):
	# wrapper script to plot a whole bunch of mean_recurrence data onto a single figure.
	if section_ids == None or (hasattr(section_ids, '__len__' and len(section_ids)==0)):
		section_ids = emc_section_filter['filter']
	#
	for sec_id in section_ids:
		do_clf=False
		if sec_id==section_ids[0]:
			do_clf=True
		#
		z=mean_recurrence(ary_in=file_path_format % sec_id, m0=m0, do_plots=True, do_clf=do_clf)
	#
	return None
#
def f_weibull2(x=None, x0=None, chi=1.1, beta=1.0):
	# prams re-ordered for different fitting algorithms.
	return f_weibull(x=x, x0=x0, chi=chi, beta=beta)
#
def f_weibull(x=None, chi=1.0, beta=1.0, x0=None):
	'''
	# weibull distribution (for fitting).
	'''
	if x0==None: x0=0.0		# ... but we give it None so the fitting algorithm does not try to fit...
	#
	return 1.0 - numpy.exp(((x0/chi)**beta) - ((x/chi)**beta))
#
def recurrence_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, keep_figs=False, output_dir='CDF_figs'):
	'''
	# for each section_id, fetch the mean_recurrence data.
	# 1) plot each cumulative probability
	# 2) including a weibull fit.
	# 3) plot each cum. prob onto a single figure
	# 4) a cumulative figure
	'''
	#
	if section_ids in (None, [], ()): section_ids = list(emc_section_filter['filter'])
	plt.ion()
	#if (-1) not in section_ids: section_ids+=[-1]		# add this special case which we will use for a composite mean fit.
	if section_ids[-1]!=-1: section_ids+=[-1]		# ... and it needs to be the final entry.
	section_ids += [-2]								# use this for the faultwise composite
	#
	# make sure we have a valid output dir:
	i=0
	new_output_dir = output_dir
	if os.path.isdir(new_output_dir)==False:
		while glob.glob(new_output_dir)!=[]:
			# the dir does not exist or ther is a file where our dir should be...
			new_output_dir = '%s_%d' % (output_dir, i)
			if i>=1000:
				print "can't find a valid output dir..."
				break
		output_dir = new_output_dir
		os.mkdir(output_dir)
	new_output_dir = None
	#
	# pressing on...
	#
	for i in xrange(len(section_ids)+1):
		# clean up a bit:
		#
		#f=plt.figure(i)
		#f.clf()
		if i==0:
			plt.figure(i)
			plt.clf()
		#
		if i>0: plt.close()
	#
	# initialize some working variables and figures...
	plt.figure(0)
	max_x = 0.
	min_x = 0.
	dT_composite_faultwise = []
	full_catalog = None
	# control color cycling:
	colors_ =  mpl.rcParams['axes.color_cycle']
	sections ={'all_cdf':{'fig':0}}
	best_fit_array = []		# ...and we'll cast this as a recarray later.
	for j, sec_id in enumerate(section_ids):
		this_color = colors_[j%len(colors_)]
		i=j+1
		sections[sec_id] = {'fig':i}
		plt.figure(i)
		#
		# set up some variables we'll need later:
		# line-widths for plotting:
		lw=2.5
		ms=5.
		#
		# fetch the data from the file. we're doing the composite wrong -- fitting basically the mean intervals on all the faults
		# in aggregate. what we want is to fit the whole collection of events together, so fetch the data and append
		# to a full catalog.
		#
		if sec_id == -1:
			# use the full_catalog...
			# but note that (i think, some events might be duplicated -- maybe when an earthquake occurs on multiple fault_segments?)
			# remove duplicates:
			full_catalog = numpy.array(list(set(full_catalog.tolist())), dtype=full_catalog.dtype)
			#
			ary_in = full_catalog
			# and we'll probably sort it in mean_recurrence(), but just in case (since we're definitely compiling
			# in an unsorted fashion):
			ary_in.sort(order='event_year')
			#
			mean_rec_data = mean_recurrence(ary_in=ary_in, m0=m0)
			X = mean_rec_data['dists']['dT'].tolist()
			lw=5.
			ms=10.
			this_lbl = '(composite mean)'
			this_lbl_composite = 'all faults (combined)'
		#
		elif sec_id == -2:
			# "faultwise" composite:
			X = dT_composite_faultwise
			lw=5.
			ms=10.
			this_lbl = '(faultwise mean)'
			this_lbl_composite = 'all faults (faultwise)'
		#
		else:
			ary_in = numpy.load(file_path_pattern % sec_id)
			#
			if full_catalog in (None, []) or len(full_catalog)==0:
				full_catalog = ary_in
			else:
				full_catalog = numpy.append(full_catalog, ary_in)
			#
			mean_rec_data = mean_recurrence(ary_in=ary_in, m0=m0)
			X = mean_rec_data['dists']['dT'].tolist()
			dT_composite_faultwise += X
			this_lbl = 'section %d' % sec_id
			this_lbl_composite = None
		#
		#
		#if sec_id != -1:
		#	mean_rec_data = mean_recurrence(ary_in=file_path_pattern % sec_id, m0=m0)
		#	#N=float(len(mean_rec_data['dists']['dT']))
		#	#X=[x for x in mean_rec_data['dists']['dT']]
		#	X = mean_rec_data['dists']['dT'].tolist()
		#	dT_composite += X
		#else:
		#	X=dT_composite
		#
		N = float(len(X))	
		X.sort()
		#
		# keep track of min/max X values for later global plot... or we could use them here as well.
		max_x = max(max_x, max(X))
		min_x = min(min_x, min(X))
		#
		Y = [x/N for x in range(1, len(X)+1)]
		#mean_dT = mean_rec_data['mean_dT']
		mean_dT = numpy.mean(X)
		#
		# get a fit:
		#. and it looks like this, by default, multiprocesses; we get all cpus to about 50% without doing any native mpp.
		fit_prams, fit_cov = spo.curve_fit(f_weibull, xdata=numpy.array(X), ydata=numpy.array(Y), p0=numpy.array([mean_dT, 1.5]))
		pram_sigmas = numpy.sqrt(numpy.diag(fit_cov))
		mean_chi_sqr = numpy.mean([(f_weibull(x=X[k], chi=fit_prams[0], beta=fit_prams[1], x0=0.0)-Y[k])**2. for k, xx in enumerate(X)]) # in xrange(len(X))])
		stdErr = mean_chi_sqr/math.sqrt(N)
		#
		print fit_cov
		print "fit_prams(%d/%d): %s, %f, %f" % (i, sec_id, str(fit_prams), mean_chi_sqr, stdErr)
		best_fit_array += [[sec_id, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr]]
		X_fit = numpy.arange(min(X), max(X)*1.5, (max(X)-min(X))/500.)
		#
		#ax=figs['all_cdf'].gca()
		f = plt.figure(0)
		plt.plot(X,Y, '.-', color=this_color)
		plt.plot(X_fit, [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=0.) for x in X_fit], '--', color=this_color, lw=lw, ms=ms, label=this_lbl_composite)
		#
		f = plt.figure(i)
		plt.plot(X,Y, '.-', color=this_color, label=this_lbl)
		plt.plot(X_fit, [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=0.) for x in X_fit], '--', color=this_color, label='$\\beta=%.3f, \\tau=%.3f, \\chi ^2=%.3e/%.3e$' % (fit_prams[1], fit_prams[0], mean_chi_sqr, stdErr), lw=lw, ms=ms)
		plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t$ (years)' % m0)
		plt.ylabel('Probability $P(t)$')
		
		sec_name = str(sec_id)
		if sec_id==-1:
			sec_name = '(mean composite)'
			#plt.gca().set_xlim([0., 4500.])
			plt.show()
		plt.gca().set_ylim([0., 1.1])
		plt.title('CDF for m>7 on fault section %s' % sec_name)
		plt.legend(loc=0, numpoints=1)
		plt.savefig('%s/VC_CDF_m%s_section_%d.png' % (output_dir, str(m0).replace('.', ''), sec_id))
		if keep_figs==False: plt.close(i)
	#
	# composite/mean fits. instead of recoding all of this, maybe use a special section_id=-1, which we add to the section_id list
	# above and do "if" on the fetch-data part...
	'''
	dT_composite.sort()
	N=float(len(dT_composite))
	Y_composite = [x/N for x in range(1, int(N)+1)]
	fit_prams, fit_cov = spo.curve_fit(f_weibull, xdata=numpy.array(dT_composite), ydata=numpy.array(Y_composite), p0=numpy.array([mean_dT, 1.5]))
	plt.figure()
	plt.plot(dT_composite, Y_composite, '.-', color='b')
	X_fit = numpy.arange(min(dT_composite), 1.2*max(dT_composite), (min(dT_composite)-1.2*max(dT_composite))/500.)
	plt.plot(X_fit, [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=0.) for x in X_fit], '-', color='b', label='$\\beta=%.3f, \\tau=%.3f$' % (fit_prams[1], fit_prams[0]))
	plt.legend(loc=0, numpoints=1)
	plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t$' % m0)
	plt.ylabel('Probability $P(t)')
	plt.savefig('CDF_figs/VC_CDF_m%s_section_mean_composite.png' % (str(m0).replace('.', '')))
	best_fit_array += [[-1, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr]]
	#
	'''
	#
	plt.figure(0)
	#X_fit = numpy.arange(min_x, max_x*1.2, (max_x-min_x)/500.)
	#plt.plot(X_fit, [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=0.) for x in X_fit], '-', lw=5, alpha=.7, zorder=7, label='mean best fit')
	#mean_tau = numpy.mean([x[1] for x in best_fit_array])
	#mean_beta = numpy.mean([x[2] for x in best_fit_array])
	#
	plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t$ (years)' % m0)
	plt.ylabel('Probability $P(t)$')
	plt.gca().set_xlim([0., 5000.])
	plt.legend(loc=0, numpoints=1)
	plt.gca().set_xscale('linear')
	plt.gca().set_yscale('linear')
	plt.gca().set_ylim([0., 1.1])
	plt.show()
	plt.savefig('%s/VC_CDF_m%s_section_composite.png' % (output_dir, str(m0).replace('.', '')))
	#
	best_fit_array = numpy.core.records.fromarrays(zip(*best_fit_array), names = ['section_id', 'tau', 'beta', 'sigma_tau', 'sigma_beta', 'mean_chi_sqr'], formats = [type(x).__name__ for x in best_fit_array[0]])
	best_fit_array.dump('%s/VC_CDF_Weibull_fits_dump.npy' % output_dir)
	#
	# and some interesting weibull pram plots are like:
	# (obviously other scalings are relevant; this seems to be the most linear)
	#plt.loglog(best_fit_array['beta'], best_fit_array['tau'], '.')
	# and of course, the distributions of beta, tau (dists. of beta are not terribly interesting,
	# but the relationship between beta and tau might be).
	'''
	plt.figure()
	ax1 = plt.gca()
	Xbeta = best_fit_array['beta'].tolist()
	Xbeta.sort()
	ax1.plot(Xbeta, range(len(best_fit_array['beta'])), 'g.-')
	ax2 = plt.gca().twiny()
	Xtau = best_fit_array['tau'].tolist()
	Xtau.sort()
	ax2.plot(Xtau, range(len(best_fit_array['tau'])), 'r.-')
	'''
	#
	return best_fit_array
	#return full_catalog
#
def waiting_time_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], keep_figs=False, output_dir='VC_CDF_WT_figs'):
	'''
	# for each section_id, fetch the mean_recurrence data.
	# 1) plot each cumulative probability
	# 2) including a weibull fit.
	# 3) do this for a set of waiting times:0, <t>/2, <t>, 1.5<t>, 2<t>, 2.5<t>
	# 4) a cumulative figure
	'''
	#
	if section_ids in (None, [], ()): section_ids = list(emc_section_filter['filter'])
	plt.ion()
	section_ids+=[-1, -2]		# ... assume they're not there and put them at the end...
	#							# (special cases for aggregate)
	#
	# make sure we have a valid output dir:
	i=0
	new_output_dir = output_dir
	if os.path.isdir(new_output_dir)==False:
		while glob.glob(new_output_dir)!=[]:
			# the dir does not exist or ther is a file where our dir should be...
			new_output_dir = '%s_%d' % (output_dir, i)
			if i>=1000:
				print "can't find a valid output dir..."
				break
		output_dir = new_output_dir
		os.mkdir(output_dir)
	new_output_dir = None
	#
	for i in xrange(len(section_ids)+1):
		# clean up a bit:
		#
		#f=plt.figure(i)
		#f.clf()
		if i==0:
			plt.figure(i)
			plt.clf()
		#
		if i>0: plt.close()
	#
	lw=2.5
	ms=5.
	max_x = 0.
	min_x = 0.
	dT_composite_faultwise = []
	full_catalog = []
	# control color cycling:
	colors_ =  mpl.rcParams['axes.color_cycle']
	sections ={'all_cdf':{'fig':0}}
	#best_fit_array = []		# ...and we'll cast this as a recarray later.
	best_fit_dict = {}
	for j, sec_id in enumerate(section_ids):
		this_color = colors_[j%len(colors_)]
		#i=j+1
		i=j
		sections[sec_id] = {'fig':i}
		plt.figure(i)
		plt.clf()
		#
		if sec_id == -1:
			# use the full_catalog...
			# but note that (i think, some events might be duplicated -- maybe when an earthquake occurs on multiple fault_segments?)
			# remove duplicates:
			full_catalog = numpy.array(list(set(full_catalog.tolist())), dtype=full_catalog.dtype)
			#
			ary_in = full_catalog
			# and we'll probably sort it in mean_recurrence(), but just in case (since we're definitely compiling
			# in an unsorted fashion):
			ary_in.sort(order='event_year')
			#
			mean_rec_data = mean_recurrence(ary_in=ary_in, m0=m0)
			X = mean_rec_data['dists']['dT'].tolist()
			lw=5.
			ms=10.
			this_lbl = '(composite mean)'
			this_lbl_composite = 'all faults (combined)'
			sec_name = 'combined mean'
		#
		elif sec_id == -2:
			# "faultwise" composite:
			X = dT_composite_faultwise
			lw=5.
			ms=10.
			this_lbl = '(faultwise mean)'
			this_lbl_composite = 'all faults (faultwise)'
			sec_name = 'faultwise mean'
		#
		else:
			ary_in = numpy.load(file_path_pattern % sec_id)
			#
			if full_catalog in (None, []) or len(full_catalog)==0:
				full_catalog = ary_in
			else:
				full_catalog = numpy.append(full_catalog, ary_in)
			#
			mean_rec_data = mean_recurrence(ary_in=ary_in, m0=m0)
			X = mean_rec_data['dists']['dT'].tolist()
			dT_composite_faultwise += X
			this_lbl = 'section %d' % sec_id
			this_lbl_composite = None
			sec_name = this_lbl
		N = float(len(X))	
		X.sort()
		#
		# keep track of min/max X values for later global plot... or we could use them here as well.
		max_x = max(max_x, max(X))
		min_x = min(min_x, min(X))
		#
		#Y = [x/N for x in range(1, len(X)+1)]
		#mean_dT = mean_rec_data['mean_dT']
		mean_dT = numpy.mean(X)
		this_t0s = [mean_dT * x for x in t0_factors]
		#
		best_fit_dict[sec_id]={}
		for t0 in this_t0s:
			max_x = 0.
			min_x = 0.
			this_X = [x-t0 for x in X if (x-t0)>=0.]
			#
			# skip it if there aren't very many data:
			if len(this_X)<5: continue
			#
			N = float(len(this_X))
			this_X.sort()
			Y = [j/N for j in range(1, len(this_X)+1)]
			max_x = max(max_x, max(X))
			min_x = min(min_x, min(X))
			#
			plt.plot([t0 + x for x in this_X], Y, '.-', color = this_color)
			# this tends to break:
			#try:
			#fit_prams, fit_cov = spo.curve_fit(f_weibull2, xdata=numpy.array([numpy.array(this_X), numpy.array([t0 for x in this_X])]), ydata=numpy.array(Y), p0=numpy.array([mean_dT, 1.5]))
			#
			# f_weibull(x=None, chi=1.0, beta=1.0, x0=None)
			fit_prams, fit_cov = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=t0), xdata=numpy.array(this_X), ydata=numpy.array(Y), p0=numpy.array([mean_dT, 1.5])) 
			
			#fit_prams, fit_cov = spo.curve_fit(f_weibull2, xdata=numpy.array([numpy.array(this_X), numpy.array([t0])]), ydata=numpy.array(Y), p0=numpy.array([mean_dT, 1.5]))
			#fit_prams, fit_cov = spo.curve_fit(f_weibull2, xdata=numpy.array(this_X), ydata=numpy.array(Y), p0=numpy.array([mean_dT, 1.5]))
			pram_sigmas = numpy.sqrt(numpy.diag(fit_cov))
			mean_chi_sqr = numpy.mean([(f_weibull(x=X[k], chi=fit_prams[0], beta=fit_prams[1], x0=0.0)-Y[k])**2. for k, xx in enumerate(this_X)]) # in xrange(len(X))])
			stdErr = mean_chi_sqr/math.sqrt(N)
			print fit_cov
			print "fit_prams(%d/%d): %s" % (i, sec_id, str(fit_prams))
			#
			best_fit_dict[sec_id][t0] = [sec_id, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr]
			X_fit = numpy.arange(min(this_X), max(this_X)*1.5, (max(this_X)-min(this_X))/500.)
			#
			
			plt.plot([t0 + x for x in X_fit], [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=t0) for x in X_fit], '--', color=this_color, lw=lw, ms=ms, label=this_lbl)
			this_chi_0  = best_fit_dict[sec_id][0.][1]
			this_beta_0 =best_fit_dict[sec_id][0.][2]
			plt.plot(X_fit, [f_weibull(x=x, chi=this_chi_0, beta=this_beta_0, x0=0.0) for x in X_fit], '--', color=this_color, lw=lw, ms=ms, label=this_lbl + ' ($t_0 = 0$)')
			#
			#except:
			#	print "fit failed. move on..."
			#	#print std.err
		#
		plt.savefig('%s/VC_CDF_WT_m%s_section_%d.png' % (output_dir, str(m0).replace('.', ''), sec_id))
		#
		plt.gca().set_ylim([0., 1.1])
		plt.title('CDF for m>7 on fault section %s' % sec_name)
		plt.legend(loc=0, numpoints=1)
		plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t$ (years)' % m0)
		plt.ylabel('Probability $P(t)$')
		
		if keep_figs==False: plt.close(i)
	#
	#
	# and some interesting weibull pram plots are like:
	# (obviously other scalings are relevant; this seems to be the most linear)
	#plt.loglog(best_fit_array['beta'], best_fit_array['tau'], '.')
	# and of course, the distributions of beta, tau (dists. of beta are not terribly interesting,
	# but the relationship between beta and tau might be).
	'''
	plt.figure()
	ax1 = plt.gca()
	Xbeta = best_fit_array['beta'].tolist()
	Xbeta.sort()
	ax1.plot(Xbeta, range(len(best_fit_array['beta'])), 'g.-')
	ax2 = plt.gca().twiny()
	Xtau = best_fit_array['tau'].tolist()
	Xtau.sort()
	ax2.plot(Xtau, range(len(best_fit_array['tau'])), 'r.-')
	'''
	#
	#return best_fit_array
	return best_fit_dict
	
#		
def mean_recurrence(ary_in='data/VC_CFF_timeseries_section_123.npy', m0=7.0, do_plots=False, do_clf=True):
	# find mean, stdev Delta_t, N between m>m0 events in ary_in.
	# for now, assume ary_in is a structured array, h5, or name of a structured array file.
	#
	if isinstance(ary_in, str)==True:
		ary_in = numpy.load(ary_in)
	# ... and for good measure, let's always sort these (just in case)... or at least try to.
	try:
		ary_in.sort(order='event_year')
	except:
		print "ary_in sorting failed, but continuing anyway, hoping for not nonsense..."
	#
	Ns, Js, Ts, Ms = zip(*[[x['event_number'], j, x['event_year'], x['event_magnitude']] for j, x in enumerate(ary_in) if float(x['event_magnitude'])>m0])
	Ns_total, Js_total, Ts_total = zip(*[[x['event_number'], j, x['event_year']] for j, x in enumerate(ary_in) ])
	#
	# just large events:
	dNs = [Ns[i]-Ns[i-1] for i, n in enumerate(Ns[1:])][1:]	# total event numbers (of the simulation, including other faults)
	dJs = [Js[i]-Js[i-1] for i, n in enumerate(Js[1:])][1:]	# sequence number along this fault.
	#dTs = [(Ts[i]-Ts[i-1])*10.**(4.5-m0) for i, t in enumerate(Ts[1:])][1:]
	dTs = [Ts[i]-Ts[i-1] for i, t in enumerate(Ts[1:])][1:]		# at some point, we were correcting this for magnitude...
	#															# but i don't recall why... so let's not...
	#
	# all events:
	stress_drop_total = [(x['cff_initial']-x['cff_final'])**2. for x in ary_in[2:]]
	stress_drop = [(x['cff_initial']-x['cff_final'])**2. for x in ary_in[2:] if float(x['event_magnitude'])>m0]	# large m stress-drop
	dNs_total = [Ns_total[i]-Ns_total[i-1] for i, n in enumerate(Ns_total[1:])][1:]
	dJs_total = [Js_total[i]-Js_total[i-1] for i, n in enumerate(Js_total[1:])][1:]
	dTs_total = [Ts_total[i]-Ts_total[i-1] for i, t in enumerate(Ts_total[1:])][1:]
	#
	recurrence_dist_array = numpy.array(zip(*[dJs, dTs]), dtype = [('dN', int), ('dT', float)])
	#
	r_dict = {'mean_dN': numpy.mean(dNs), 'stdev_dN':numpy.std(dNs), 'mean_dT':numpy.mean(dTs), 'stdev_dT':numpy.std(dTs), 'mean_dN_fault':numpy.mean(dJs), 'stdev_dN_fault':numpy.std(dJs), 'dists':recurrence_dist_array}
	#
	#
	if do_plots:
		#print len(dTs_total), len(dJs_total), len(stress_drop_total)
		#print len(dTs), len(dJs), len(stress_drop)
		plt.figure(1)
		if do_clf: plt.clf()
		ax = plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		plt.plot(dTs_total, stress_drop_total, '.', zorder=4)
		plt.plot( dTs, stress_drop[-len(dNs):],'.', zorder=5)
		plt.xlabel('interval $\\Delta t$')
		plt.ylabel('stress drop $(CFF_{final} - CFF_{initial})^2$')
		#
		plt.figure(2)
		if do_clf: plt.clf()
		ax = plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		plt.plot(dNs_total, stress_drop_total, '.')
		plt.plot(dNs, stress_drop[-len(dNs):], '.')
		#plt.plot( stress_drop[2:], dJs,'.')
		plt.xlabel('event-interval, $\\Delta N$')
		plt.ylabel('stress drop $(CFF_{final} - CFF_{initial})^2$')
		#
		#print r_dict
	#
	return r_dict
#
def get_trend_analysis(ary_in=None, nyquist_len=10, nyquist_time=None):
		if isinstance(ary_in, str)==True:
			ary_in = numpy.load(ary_in)
		#
		Xs = ary_in['event_year'][1:]
		Ns = ary_in['event_number'][1:]
		intervals = numpy.array(map(math.log10, ary_in['event_year'])[1:]) - numpy.array(map(math.log10, ary_in['event_year'])[:-1])
		mean_interval = numpy.mean(intervals)
		#
		X_numpy = numpy.array([[x, 1.0] for x in Xs])
		
		# first, get a fixed length line-fit:
		#for i in xrange(nyqquist_len, len(ary_in)):
		#	rw = ary_in[i]
		fitses = [numpy.linalg.lstsq(X_numpy[i-nyquist_len:i], intervals[i-nyquist_len:i])[0] for i in xrange(nyquist_len, len(ary_in))]
		#
		output_names = ['event_number', 'event_year', 'lin_fit_a', 'lin_fit_b']
		#
		# get record-breaking intervals
		nrbs=[]
		output_names += ['rb_ratio']
		for i in xrange(nyquist_len, len(ary_in)):
			#
			nrbs_up=[intervals[i-nyquist_len]]
			nrbs_dn=[intervals[i-nyquist_len]]
			#
			for j, interval in enumerate(intervals[i-nyquist_len:i]):
				if interval > nrbs_up[-1]: nrbs_up+=[interval]
				if interval < nrbs_dn[-1]: nrbs_dn+=[interval]
			#
			nrbs += [math.log10(float(len(nrbs_up))/float(len(nrbs_dn)))]
			#rb_ratio = float(len(nrbs_up))/float(len(nrbs_dn))
			#outputs[i]+=[rb_ratio]
			#
		#nrbs = [None for i in xrange(nyquist_len)] + nrbs
		#
		outputs = [[Ns[i], Xs[i], fitses[i][0], fitses[i][1], nrbs[i], intervals[i]/mean_interval] for i in xrange(len(fitses))]
		output_names += ['interval_rate_ratio']
		print "lens:: ", len(nrbs), len(outputs)
		#
		#CFF = numpy.core.records.fromarrays(CFF.transpose(), names=['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final'], formats=[type(x).__name__ for x in CFF[0]])
		outputs = numpy.core.records.fromarrays(zip(*outputs), names=output_names, formats = [type(x).__name__ for x in outputs[0]])
		#
		return outputs			
#
def plot_CFF_ary(ary_in='data/VC_CFF_timeseries_section_125.ary', fnum=0, nyquist_factor=.5):
	# this script is for some of the earlier CFF numpy array types. newer arrays will require different scripting.
	# for older data sets (most likely):
	#	# 2 cols: event_number, CFF_initial
	#   # 3 cols: event_number, event_year, CFF_initial
	#   # 4 cols: event_number, event_year, CFF_initial, CFF_final
	#  later versions will include column names.
	#
	CFF = numpy.load(ary_in)
	#
	recurrence_data = mean_recurrence(ary_in=CFF, m0=7.0)
	nyquist_len = int(nyquist_factor*recurrence_data['mean_dN_fault'])
	nyquist_time = nyquist_factor*recurrence_data['mean_dT']
	#
	trend_data = get_trend_analysis(ary_in=CFF, nyquist_len = nyquist_len, nyquist_time=nyquist_time)
	#
	# what kind of array did we get?
	# this is not very efficient, in that we rewrite the whole enchilada for the new types, but i don't
	# expect that we'll be returning to these unformatted (unstructured) array types.
	# so, focus on newer structured arrays. leave the old thing in out of principle.
	#	
	if isinstance(CFF, numpy.recarray)==True:
		# it's a structured array.
		#cols = map(operator.itemgetter(0), CFF.dtype.descr)
		col = CFF.dtype.names
		# cols should be like: ['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'event_area']
		#
		f=plt.figure(fnum)
		f.clf()
		#
		# create two axes:
		# magnitudes plot
		ax_mag = f.add_axes([.1, .05, .85, .25])
		# CFF plot.
		ax_CFF = f.add_axes([.1, .35, .85, .25], sharex=ax_mag)
		ax_dCFF = ax_CFF.twinx()	# over-plot stress (CFF) drop...
		ax_ints = f.add_axes([.1, .65, .85, .25], sharex=ax_mag)
		ax_mag2 = ax_ints.twinx()
		ax_trend=ax_mag.twinx()
		#
		X_init = CFF['event_year']
		X_finals = [x+.01 for x in X_init]
		#
		Y0 = -1*CFF['cff_initial']
		Y_final = -1*(CFF['cff_final'])
		
		X = list(X_init) + list(X_finals)
		X.sort()
		#
		intervals = X_init[1:] - X_init[:-1]
		
		big_mags = zip(*[[x['event_year'], x['event_magnitude']] for x in CFF if x['event_magnitude']>7.0])
		#
		Y = []
		CFF_drops = [(x['cff_initial'] - x['cff_final'])**2. for x in CFF]
		for i, y in enumerate(Y0):
			Y += [Y0[i]]
			Y += [Y_final[i]]
		#
		# use "peak" values to cut through some noise.
		peaks = get_peaks(zip(*[X,Y]), col=1, peak_type='upper')
		X_peaks, Y_peaks = zip(*peaks)
		#
		# CFF Plot:
		#ax = plt.gca()
		ax_CFF.set_xscale('linear')
		ax_CFF.set_yscale('log')
		ax_CFF.set_ylabel('CFF')
		ax_dCFF.set_xscale('linear')
		ax_dCFF.set_yscale('log')

		# first, raw CFF (initial):
		ax_CFF.plot(X, Y, '.-', color='b', alpha=.2, zorder=4)
		ax_CFF.fill_between(X, Y, y2=min(Y), color='b', alpha=.2, zorder=4)
		ax_CFF.plot(X_peaks, Y_peaks, '-', zorder=5)
		ax_dCFF.plot(CFF['event_year'], CFF_drops, 'g.-', zorder=7, alpha=.9)
		#
		# Magnitude plot (by itself):
		ax_mag.set_xscale('linear')
		ax_mag.set_yscale('linear')
		ax_mag.set_ylabel('event magnitude $m$')
		ax_mag.set_xlabel('event year $t$')
		min_mag = min(CFF['event_magnitude']) - .5
		ax_mag.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='b', alpha=.9)
		#
		ax_trend.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], 'r-', zorder=5, alpha=.8)
		ax_trend.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='r', zorder=2, alpha=.8)
		ax_trend.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		ax_trend.set_ylabel('(log) interval slope $b$')
		#
		ax_trend2 = ax_ints.twinx()
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], 'r-', zorder=5, alpha=.8)
		#
		ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b']  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='m', zorder=1, alpha=.5)
		ax_trend2.fill_between([x['event_year'] for x in trend_data], [1.  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0.0 for x in trend_data], color='m', zorder=1, alpha=.25)
		#ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b'] / x['interval_rate_ratio'] for x in trend_data], y2=[0.0 for x in trend_data], where=[(x['lin_fit_b'] / x['interval_rate_ratio'])<0. for x in trend_data], color='m', zorder=1, alpha=.5)
		
		ax_trend2.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		#
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['rb_ratio'] for x in trend_data], 'c--')
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] + x['rb_ratio'] for x in trend_data], 'c-')
		
		#
		# Interval Plots:
		ax_ints.set_xscale('linear')
		ax_ints.set_yscale('log')
		ax_ints.set_ylabel('intervals $\\Delta t$')
		ax_ints.plot(X_init[1:], intervals, '.-', alpha=.9)
		ax_mag2.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='g', alpha=.9, lw=2)
		ax_mag2.vlines(big_mags[0], [min_mag for x in big_mags[1]], big_mags[1], color='r', lw=2.5, alpha=.9)
		ax_mag2.set_ylabel('magnitude $m$')
	#	
	if isinstance(CFF, numpy.recarray)==False:
		# a regular, old-style, numpy.ndarray -- aka, no columns. guess the column structure from what we know...
		#
		# assume either [event_id/num, year, CFF] or [event_id/num, CFF]
		if len(CFF[0])==2:
			y_col=1
		if len(CFF[0])>=3:
			x_col=1
			y_col=2
		if len(CFF[0])>=5:
			x_col=1
			y_col=3
		#CFF_peaks = get_peaks(data_in=CFF, col=y_col, peak_type='lower')
		#
		zCFF = zip(*CFF)
		X=zCFF[x_col]
		Y=zCFF[y_col]
		Y=[-1*y for y in Y]
		CFF_peaks = get_peaks(data_in = zip(*[X,Y]), col=1, peak_type='upper')
		#
		zCFF_peaks = zip(*CFF_peaks)
		X_peaks = zCFF_peaks[0]
		Y_peaks = zCFF_peaks[1]
		#
		'''
		if len(CFF[0])==2:
			X=zCFF[0]
			Y=zCFF[1]
		if len(CFF[0])>=3:
			X=zCFF[1]
			Y=zCFF[2]
		'''
		#
		plt.figure(fnum)
		plt.clf()
		#
		ax = plt.gca()
		ax.set_xscale('linear')
		ax.set_yscale('log')
		plt.fill_between(X, Y, y2=min(Y), color='b', alpha=.2, zorder=4)
		plt.plot(X_peaks, Y_peaks, '.-', zorder=5)				
	#
	return CFF
#
# end CFF calculators and helpers...	
#
def get_EMC_CFF(sections=None):
	if sections==None: sections = [16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 69, 70, 73, 83, 84, 92, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 123, 124, 125, 126, 149]
	#
	event_sweeps_dict = None
	#with h5py.File(sim_file) as vc_data:
	#	# pre-calc a dictionary of event-sweeps. each event_number will include a list of sweep events.
	#	event_sweeps_dict = make_h5_indexed_dict_spp(h5_in=vc_data['event_sweep_table'], index_col='event_number')
	#
	for section in sections:
		X=get_CFF_on_section(section_id=section, event_sweeps_dict=event_sweeps_dict)
		X=numpy.array(X)
		X.dump('data/VC_CFF_timeseries_section_%d.npy' % section)
		#
		# ... and we may have redefined the offending variable types by using the structured
		# numpy array, so let's have a go at json as well:
		#with open('data/VC_CFF_timeseries_section_%d.json' % section, 'w') as f_out:
		#	#... well crap, these won't json.dump(s)() because the 'integer' values in the data are not
		#	# proper integers; they numpy.uint32 (i presume  numpy unsigned integers, 32 bit), and json
		#	# does not know what to do with them. we could try to handle them, or we can forego json...
		#	json.dump(X, f_out, indent=1)
		#
	return None
#
def get_CFF_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, event_sweeps_dict=None):
	# CFF = shear_stress - mu*normal_stress
	#  blocks: dictionary of blocks indexed by event number: {<event_number>:[block_data, block_data, etc.], <>:[]}
	# for runs of multiple instances, this can be pre-calculated and shared.
	#
	# 1) get events for the section (use "events_by_section")
	# 2) for each event, gather all the blocks (use the event_sweep_table).
	# 3) for each block-sweep in the event, add the local CFF (aka, add the shear_init + mu*normal_init for each entry in the event.
	# 4) time (year) of the event can be taken from the event_table.
	#
	time_start = time.time()
	print "starting CFF(t) for section: %d (%s)" % (section_id, time.ctime())
	if n_cpus==None: n_cpus = mpp.cpu_count()
	CFF=[]
	events = get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=section_id, n_cpus=n_cpus)
	t0=time.time()
	#
	# not using this any longer...
	#col_dict = {}
	#map(col_dict.__setitem__, events[0], range(len(events[0])))
	#col_dict = {event:i for i, event in enumerate(events[0])}
	#
	# pre-load blocks data:
	# (actually, we don't need this all the required data are in "sweeps").
	#blocks_data = get_blocks_info_dict(sim_file=sim_file, block_ids=None)
	#
	# so, events is a time series of events (earthquakes) on a particular fault section.
	# for each event, the normal/shear stresses are calculated.
	#
	# now, fetch block level data. for now, let's return a simpler data set than we initially read.
	#
	with h5py.File(sim_file) as vc_data:
		sweep_data = vc_data['event_sweep_table']
		max_index = len(events)-1
		#
		# pre-calc a dictionary of event-sweeps. each event_number will include a list of sweep events.
		#if event_sweeps_dict==None: 
		#	print "setting up sweeps dictionary..."
		#	event_sweeps_dict = make_h5_indexed_dict_spp(h5_in=sweep_data, index_col='event_number')
		#
		# make some mpp queues:
		jobs     = []
		#pipes    = []
		#cff_outs = []
		for ev_count, event in enumerate(events[1:]):
			# get_event_blocks(sim_file=allcal_full_mks, event_number=None, block_table_name='block_info_table')
			#
			# "events" is a list type, but each row is some sort of numpy recarray (returns numpy.void), but acts like rec_array
			event_id=event['event_number']
			event_year = event['event_year']
			event_area = event['event_area']
			# this needs to be sped up maybe -- perhaps by maintaining the hdf5 context?
			#
			#blocks = fetch_h5_data(sim_file=sim_file, n_cpus=n_cpus, table_name='event_sweep_table', col_name='event_number', matching_vals=[event_id])
			#
			# this should be faster, since we don't have to open and reopen the file... but i'm not sure how much
			# faster it is. certainly for smaller jobs, the sipler syntax is probably fine.
			# ... but still slow, so use a pre-indexed list of sweep events...
			if event_sweeps_dict ==None:
				blocks = fetch_h5_table_data(h5_table=sweep_data, n_cpus=None, col_name='event_number', matching_vals=[event_id])
			else:
				# ... though i'm starting to think that this dict. object is not any faster (in fact likely slower)
				# than (repeatedly) accessing the h5 table directly.
				blocks = event_sweeps_dict[event_id]
			#
			if ev_count%250==0:
				#print "(%d) blocks (%d) fetched: %d" % (ev_count, event_id, len(blocks))
				t1=time.time()
				print "fetching CFF events: %d/%d, %d blocks fetched, dt=%s" % (ev_count, max_index, len(blocks), str(t1-t0))
				t0=t1
			#
			# ... and now get the block info for each block in the event (from the blocks_data dict).
			# ... but more specifically, let's calculate the cumulative CFF for this event (summing the local CFF for each block).
			#
			#CFF+=[[event_id, event_year, event['event_magnitude'], get_event_CFF(blocks), get_final_event_CFF(blocks)]]
			#
			# ... and this approach appears to be killing us. i think mpp is simply not efficient with this generalized approach
			#  (aka, doing the mpp on each set of blocks.). it may be faster and more memory efficient
			# to do each set of blocks as a process, either by mapping (which i think will kill us in memory)
			# or by using processes (and for now, just suck up that larger chunks are going to slow down smaller chunks
			# inthe join).
			#
			# ... and instead of returning a list of dicts, let's return a structured numpy.array().
			# we can define the column names and types; see below.
			#
			pipe1, pipe2 = mpp.Pipe()
			P=mpp.Process(target=calc_CFFs_from_blocks, args=[blocks, pipe1])
			P.start()
			#jobs+=[{'process':P, 'pipe':pipe2, 'cff_out':{'event_id':event_id, 'event_year':event_year,'event_magnitude':event['event_magnitude']}}]
			jobs+=[{'process':P, 'pipe':pipe2, 'cff_out':[event_id, event_year, event['event_magnitude']]}]
			#pipes+=[pipe2]
			#cff_outs += [{'event_id':event_id, 'event_year':event_year,'event_magnitude':event['event_magnitude']}]
			#
			# we're spinning through this, so we should have n>1 processes running. let's slow them down...
			if ev_count%n_cpus==(n_cpus-1) or ev_count==max_index:
				#print "do joins..."
				for j, job in enumerate(jobs):
					#print "from pipe[%d]: " % (j), pipes[j].recv()
					cff_vals = job['pipe'].recv()
					job['pipe'].close()
					#job['cff_out'].update(cff_vals)		# cff_vals like: {'cff_init':total_cff_init, 'cff_final':total_cff_final}
					#job['cff_out'] += [cff_vals['cff_init'], cff_vals['cff_final'], cff_vals['mean_cff_init'], cff_vals['mean_cff_final']]
					# ... but instead of the "block-mean" value returned by calc_CFFs_from_blocks, let's use the event_area:
					job['cff_out'] += [cff_vals['cff_init']/event_area, cff_vals['cff_final']/event_area, event_area]
					#
					CFF+=[job['cff_out']]
				#
				#print "waiting on pipes..."
				# and just to be sure, join until we're all done with this set:
				for j, job in enumerate(jobs):
					job['process'].join()
				#print "joins complete"
				jobs=[]
			#
			##CFF+=[{'event_id':event_id, 'event_year':event_year,'event_magnitude':event['event_magnitude'], 'CFF':get_event_CFF(blocks), 'CF_final':get_final_event_CFF(blocks)}]
			#CFF+=[{'event_id':event_id, 'event_year':event_year,'event_magnitude':event['event_magnitude']}]
			#CFF['CFF']      = get_event_CFF(blocks, chunksize=128)
			#CFF['CF_final'] = get_final_event_CFF(blocks, chunksize=128)		# process quasi-serially...
			#
			# calc_CFFs_from_blocks(sweep_blocks)
			#
			blocks=None
		#
	# convert CFF to a structured array (basically, we're after named columns; in this case, defining the data types is
	# easy as well).
	# this is what we'll want to do, but for some reason it doesn't work worth a damn. just send back the list for now...
	#CFF = numpy.array(CFF, dtype=[('event_number', 'int32'), ('event_year', 'int32'), ('event_magnitude', 'float64'), ('CFF_init', 'float64'), ('CFF_final', 'float64')])
	#
	print "finished CFF(t) for section: %d (%s: %f)" % (section_id, time.ctime(), time.time()-time_start)
	#
	# convert to structured array with named cols (note this syntax, because this is not as easy as it should be):
	CFF = numpy.core.records.fromarrays(CFF.transpose(), names=['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'event_area'], formats=[type(x).__name__ for x in CFF[0]])
	#
	return CFF
#
def make_structured_arrays(file_profile = 'data/VC_CFF_timeseries_section_*.npy'):
	# wrapper to convert a bunch of normal arrays or maybe lists to numpy structured arrays (numpy.recarray).
	G=glob.glob(file_profile)
	#
	for g in G:
		print "fixing file: %s" % g
		try:
			z=make_structured_array_from_file(fname_in=g, fname_out=g)
			print "successful..."
		except:
			print "failed to 'fix' file %s. it might have already been converted..." % g

def fix_CFF_in_struct_arrays(file_profile = 'data/VC_CFF_timeseries_section_*.npy', h5file = allcal_full_mks):
	# ALMOST CERTAINLY A ONE-TIME JOB SCRIPT...
	#
	# when the CFF time series were initially compiled, we summed the CFF functions but did not provide a mean value,
	# so the CFF is sort of nonsense. <CFF> can be calculated from the event_table['event_area'], or
	# it can be caluclated from the length of the sweep-subset -- which is probably how it would be caluclated in real-time,
	# so for now, let's just do that (recognizing that it will take a bit longer).
	G=glob.glob(file_profile)
	#
	with h5py.File(h5file, 'r') as f:
		for g in G:
			print "fixing file: %s" % g
			#
			with open(g, 'r') as f_read:
				datas = numpy.load(f_read)
			new_data=[]
			names, types = zip(*datas.dtype.descr)
			#
			# now, for each event in the data:
			# get the event_number and corresponding event data
			# update the CFF as CFF/Area
			# add area to the data
			# then, output a new structured array.
			#
			for i, rw in enumerate(datas):
				event_number = rw['event_number']
				event_data = f['event_table'][event_number]	# noting that event_number will be the index of this table (i=event_number)
				event_area = event_data['event_area']
				#
				new_data_row = rw.copy()
				new_data_row['cff_initial']/=event_area
				new_data_row['cff_final']/=event_area
				#
				new_data += [list(new_data_row.tolist())]
				#print "new data: ", new_data[-1]
				new_data[-1] += [event_area]
			#
			# make structured array from new_data
			new_data = numpy.array(new_data)
			new_data = numpy.core.records.fromarrays(new_data.transpose(), names=['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'event_area'], formats=[type(x).__name__ for x in new_data[0]])
			#
			with open(g, 'w') as f2:
				new_data.dump(f2)
			#
			print "fixed (hopefully) file: ", g
		
def make_structured_array_from_file(fname_in, col_names = ['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'mean_cff_init', 'mean_cff_final'], col_formats = None, fname_out=None):
	# note: this yields a numpy.recraray, as opposed to the standard numpy.ndarray type/instance. note that, given
	# recarray A and ndarray B:
	#
	# isinstance(A, numpy.ndarray):  True
	# isinstance(B, numpy.ndarray):  True
	# isinstance(A, numpy.recarray): True
	# isinstance(B, numpy.recarray): False
	#
	# so, recarray inherits ndarray.
	#
	with open(fname_in, 'r') as f:
		X  = numpy.load(f)
		X=numpy.array(X.tolist())		# if we want to modify an existing array, we'll probably need to convert it first.
		#
		if col_formats == None: col_formats = [type(x).__name__ for x in X[0]]
		#
		Xs = numpy.core.records.fromarrays(X.transpose(), names=col_names, formats=col_formats)
		#
	#
	if fname_out==None:
		dot_index = fname_in.rfind('.')
		fname_out = fname_in[0:dot_index] + '_struct' + fname_in[dot_index:]
	#
	with open(fname_out, 'w') as f:
		Xs.dump(f)
	#
	return Xs
	
		
#
def h5_index_event_number(h5_in=None, index_col='event_number'):
	# this is a little helper funciton to facilitate map() and other mpp scripts.
	# ... though we might be able to use initialize() pool kwarg. for now...
	return make_h5_indexed_dict_spp(h5_in=h5_in, index_col=index_col)
#
def make_h5_indexed_dict_spp(h5_in=None, index_col=None, pipe_out=None, data_set_name=None):
	# worker for make_h5_indexed_dict (or can be used standalone).
	# group rows by unique value(s) in index_col
	# @data_set_name: use only if h5_in is a string defining the name of the hdf5 filename.
	#
	# we might want to call this outside an HDF5 context (maybe we want to do this always...)
	# so, if the h5_in data type is a string, assume it's a file-name, and assume we get a not-None data_set_name
	# open the file
	if isinstance(h5_in, str):
		# we've (probably) been passed the name of the hdf5 source file. open it and call this function recursively:
		with h5py.File(h5_in, 'r') as f:
			return make_h5_indexed_dict_spp(h5_in=f[data_set_name], index_col=index_col, pipe_out=pipe_out)
	#
	#
	print "setting up keys:"
	#[return_dict.__setitem__(key, []) for key in key_vals]
	return_dict = {key:[] for key in set(h5_in[index_col])}
	print "keys configured. assign rows:"
	[return_dict[rw[index_col]].append(rw) for rw in h5_in]
	#
	if pipe_out!=None:
		print "f_piping out...(%d)" % len(return_dict)
		pipe_out.send(return_dict)
		print "piped. now close..."
		pipe_out.close()
	else:
		return return_dict

	
#
# this class could be used to (re)-index a table, aka grouping into a dictionary based on common values in index_col.
# this approach would use individual Process() objects, and the source data would be partitioned (evenly) to 
# n_cpu processes, aka: h5_indexed_dict_worker(data_in=[data:N/n_cpu]), ...
class h5_indexed_dict_worker(mpp.Process):
	#
	def __init__(self, data_in, index_col=None, pipe_out=None):
		mpp.Process.__init__(self)
		#
		self.data_in=data_in
		self.index_col=index_col
		#
		self.pipe_out = pipe_out
		#
		self.run = self.index_data	# but now this can be changed if we want...
		print "h5 indexer initted..."
	#
	def index_data(self, data_in=None, index_col=None, pipe_out=None):
		if index_col == None: index_col = self.index_col
		if data_in == None: data_in = self.data_in
		if pipe_out == None: pipe_out = self.pipe_out
		#
		#return_dict={}
		#
		# this assumes hdf5 type input...
		key_vals = set(data_in[index_col])
		print "setting up keys (%d)..." % len(key_vals)
		#[return_dict.__setitem__(key, []) for key in key_vals]
		return_dict = {key:[] for key in key_vals}
		print "keys set. asign..."
		#[return_dict[rw[index_col]].append(rw) for rw in data_in]
		#for rw in data_in:
		#	return_dict[rw[index_col]] += [rw]
		#print return_dict
		#
		if pipe_out==None:
			return return_dict
		else:
			try:
				print "piping out... (%d)" % len(return_dict), type(pipe_out)
				pipe_out.send(return_dict)
				#pipe_out.send({'a':1, 'b':2})
				print "... and close..."
				pipe_out.close()
				print "and pipe out and closed..."
				#
			except:
				print 'excepting...'
				return return_dict
	#
#
def make_h5_indexed_pool(h5_table=None, n_cpus=None, chunksize=None):
	# incidentals:
	index_col='event_number'
	if n_cpus==None: n_cpus = mpp.cpu_count()
	data_len = len(h5_table)
	job_indices = [i*data_len/n_cpus for i in xrange(n_cpus)] + [data_len]
	if chunksize==None: chunksize = len(h5_table)/n_cpus
	#
	print "job_indices: ", job_indices
	#
	jobs = []
	pipes = []
	output_dict={key:[] for key in set(h5_table[index_col])}
	#
	pool =  mpp.Pool(processes=n_cpus)
	#
	for i in xrange(n_cpus):
		jobs += [h5_table[job_indices[i]:job_indices[i+1]]]
	#
	results = pool.map_async(h5_index_event_number, jobs, chunksize=chunksize)
	pool.close()
	pool.join()
	#
	# [output_dict[key].extend(tmp_dict[key]) for key in tmp_dict.keys()]
	returns = results.get()
	for dct in returns:
		[output_dict[key].extend(dct[key]) for key in dct.keys()]
	#
	return output_dict
	
	
#
def make_h5_indexed_dict(h5_table=None, n_cpus=None, index_col=None, use_obj=False):
	'''
	# (mpp implied)
	#
	# (re)-index a data set. for example, for event_sweep_table, take a set like
	# {event_id, blah, blah}:[ [0, blah0], [0, blah1], [0, blah2], [1, blah0], [1, blah1]...]
	# and make: [{0:[[0, blah0], [0, blah1], [0, blah2]]}, {1:[[1, blah0], [1, blah1]], etc.  ]
	# now, events can be extracted easily and indexidly.
	#
	# and we find this weird behavior where out-sourcing the process seems to always improve speed performance, so let's
	# always mpp... we should bench this against Pool() as well (original test were with direct Process() implementation.
	'''
	#
	# incidentals:
	if n_cpus==None: n_cpus = mpp.cpu_count()
	data_len = len(h5_table)
	job_indices = [i*data_len/n_cpus for i in xrange(n_cpus)] + [data_len]
	#
	print "job_indices: ", job_indices
	#
	jobs = []
	pipes = []
	output_dict={key:[] for key in set(h5_table[index_col])}
	#
	for i in xrange(n_cpus):
		p1, p2 = mpp.Pipe()
		pipes += [p1]
		#print "for indices: ", job_indices[i], job_indices[i+1]
		#
		if use_obj:
			jobs  += [h5_indexed_dict_worker(data_in=h5_table[job_indices[i]:job_indices[i+1]], index_col=index_col, pipe_out=p2)]
		else:
			jobs += [mpp.Process(target=make_h5_indexed_dict_spp, kwargs={'h5_in':h5_table[job_indices[i]:job_indices[i+1]], 'index_col':index_col, 'pipe_out':p2})]
		jobs[-1].start()	
	#jobs = [h5_indexed_dict_worker(data_in=h5_table[job_indices[i]:job_indices[i+1]], index_col=index_col) for i in xrange(n_cpus)]
	#
	#
	# and wait until these jobs finish...
	# and then fetch the results.
	# note: we've possibly split up some results, so we may need to append. nominally, these will be at the begining/end
	# of subsets. we could handle it in processor assignment, or we can do it by pre-indexing the dictionary (like we do
	# in the worker(s) ).
	#output_dict = {key:[] for key in set(h5_table[index_col])}	# gives an initialized dict. with empty value lists.
	#
	for i, p in enumerate(pipes):
		#output_dict.update(pipes[i].recv())
		# this doesn't appear to be too time-costly... but, at lest once the system finds its indices, 
		#the spp seems to be much faster...; the fast majority of the 'cost' is in piping back the data.
		# pool.map_async() might be faster
		tmp_dict = pipes[i].recv()
		[output_dict[key].extend(tmp_dict[key]) for key in tmp_dict.keys()]
		pipes[i].close()
	#
	for i, job in enumerate(jobs):
		print "joining: %d" % i
		jobs[i].join()
				
	return output_dict
#
def index_dict_test(N=10**6):
	# test mpp index maker...
	# and the results... it's much, much faster to do dictionary assignment using a single processor.
	# mpp can be marginally faster to assign values over 10**6 values or so, but the cost of piping back
	# the data is HUGE. basically, this test takes about 2 seconds in spp. this process takes about
	# 1.5-1.7 seconds in at least two of the mpp modes, but piping back the data puts the whole process at
	# ~ 7 seconds, 10 seconds using the pool.map_async() method.
	#
	with h5py.File(allcal_full_mks,'r') as f:
		sweeps = f['event_sweep_table']
		if N==None: N=len(sweeps)
		times=[]
		times+=[time.time()]
		#
		sweeps1 = make_h5_indexed_dict_spp(h5_in=sweeps[0:N], index_col='event_number')
		times+=[time.time()]
		#
		sweeps2 = make_h5_indexed_dict(h5_table=sweeps[0:N], index_col='event_number', use_obj=False)
		times+=[time.time()]

		sweeps3 = make_h5_indexed_dict(h5_table=sweeps[0:N], index_col='event_number', use_obj=False)
		times+=[time.time()]
		#
		sweeps4 = make_h5_indexed_pool(h5_table=sweeps[0:N], n_cpus=None)
		times+=[time.time()]
		#
		for i in xrange(1, len(times)):
			print "time_%d: %f" % (i, times[i]-times[i-1])
		print "lens (1,2,3, 4): %d, %d, %d, %d" % (len(sweeps1), len(sweeps2), len(sweeps3), len(sweeps4))
		equal_tests = [[sweeps1, sweeps2], [sweeps2,sweeps3], [sweeps3,sweeps1], [sweeps1,sweeps4], [sweeps2,sweeps4], [sweeps3,sweeps4]]
		#
		
		# process() MPP eqality test:
		# (this approach is definitely more memory efficient and tractable than the pool().map_async() method.
		# it might also be worth looking into a Queue() method, but generally it seems that 
		#
		jobs=[]
		pipes=[]
		n_procs=mpp.cpu_count()
		print "equalities (in some order):"
		for i, rw in enumerate(equal_tests):
			if len(rw)<2: continue
			#P=mpp.Process(target=operator.eq, args=[rw[0], rw[1]])
			pipe1, pipe2 = mpp.Pipe()
			P=mpp.Process(target=print_eq, args=[rw[0], rw[1], pipe2])
			P.start()
			jobs+=[P]
			pipes+=[pipe1]
			# ... and of course, we could be cleaning these up as well, and probably monitoring them and spawning new processes
			# as they finish (using the pipes() maybe?
			# .. but we need better tracking of how many processes we started.
			if i%n_procs==(n_procs-1) or i==len(equal_tests)-1:
				print "do joins..."
				for j, job in enumerate(jobs):
					#
					print "from pipe[%d]: " % (j), pipes[j].recv()
					pipes[j].close()
				#
				for job in jobs: job.join()
				jobs=[]
				pipes=[]
				#
			#
		#
		equal_tests=None
		#
		'''
		# testing pool method for evaluating equality. it is super memory intensive and problematic. maybe try direct processes?
		3
		#print "equalsies in: ", [[operator.eq] + [S] for S in equal_tests]
		arg_list = [[operator.eq] + [S] for S in equal_tests]
		for l in arg_list: print l[0], len(l), len(l[1])
		#
		eq_pool = mpp.Pool(mpp.cpu_count())
		equalses = eq_pool.map_async(map_helper, arg_list)
		#print "equalities: ", sweeps1==sweeps2, sweeps2==sweeps3, sweeps3==sweeps1, sweeps1==sweeps4, sweeps2==sweeps4, sweeps3==sweeps4
		results = equalses.get()
		eq_pool.close()
		print "equalities: ", results
		'''
		#
	return None
#
def print_eq(a,b, pipe=None):
	truth = (a==b)
	if pipe!=None:
		pipe.send(truth)
		pipe.close()
	print truth
#
def get_stress_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, fignum=0):
	# ... and "time_series" is implied.
	# get time series of stress on a block. plot 2 axes on one plot: normal_stress, shear_stress.
	# note: the CFF, what we'll ultimately want, is:
	# CFF_j = shear_stress_j(t) - mu_j*normal_stress_j(t)
	#
	events = get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=section_id, n_cpus=n_cpus)
	col_dict = {}
	map(col_dict.__setitem__, events[0], range(len(events[0])))
	#
	# each time step has an initial + final stress value.
	ts_single = []
	mags = []
	ts = []
	shear_stress  = []
	normal_stress = []
	ave_slip = []
	for rw in events[1:]:
		ts += [rw[col_dict['event_year']], rw[col_dict['event_year']]]
		ts_single += [rw[col_dict['event_year']]]
		mags += [rw[col_dict['event_magnitude']]]
		shear_stress += [rw[col_dict['event_shear_init']], rw[col_dict['event_shear_final']]]
		normal_stress += [rw[col_dict['event_normal_init']], rw[col_dict['event_normal_final']]]
		ave_slip += ['event_average_slip']
	#
	#
	# plotting bits:
	myfig = plt.figure(fignum)
	myfig.clf()
	#
	myfig.add_axes([.05, .05, .9, .4], label='shear')
	myfig.add_axes([.05, .5, .9, .4], sharex=myfig.axes[0])
	#
	#ax = plt.gca()
	for ax in myfig.axes:
		ax.set_xscale('linear')
		ax.set_yscale('log')
		ax.set_ylabel('stress (VC units)', size=14)
		#
	#
	# let's look at peak values...
	upper_shear_init = get_peaks(events[1:], col=col_dict['event_shear_init'], peak_type='upper')
	upper_shear_final = get_peaks(events[1:], col=col_dict['event_shear_final'], peak_type='upper')
	#print "lens: ", len(upper_shear_init), len(upper_shear_final)
	#print upper_shear_init[-5:]
	#print upper_shear_final[-5:]
	lower_shear_init = get_peaks(events[1:], col=col_dict['event_shear_init'], peak_type='lower')
	lower_shear_final = get_peaks(events[1:], col=col_dict['event_shear_final'], peak_type='lower')
	#
	upper_normal_init = get_peaks(events[1:], col=col_dict['event_normal_init'], peak_type='upper')
	upper_normal_final = get_peaks(events[1:], col=col_dict['event_normal_final'], peak_type='upper')
	lower_normal_init = get_peaks(events[1:], col=col_dict['event_normal_init'], peak_type='lower')
	lower_normal_final = get_peaks(events[1:], col=col_dict['event_normal_final'], peak_type='lower')
	#
	ax = myfig.axes[0]
	ax.plot(ts, shear_stress, '.-', label='shear stress')
	for tbl, col,lbl in [(upper_shear_init, 'event_shear_init', 'init_upper'), (upper_shear_final, 'event_shear_final', 'final_upper'), (lower_shear_init, 'event_shear_init', 'init_lower'), (lower_shear_final, 'event_shear_final', 'final_lower')]:
		#ax.plot(map(operator.itemgetter(col_dict['event_year'], upper_shear_init)), map(operator.itemgetter(col_dict['event_shear_init'], upper_shear_init)),
		#print tbl[0]
		X = map(operator.itemgetter(col_dict['event_year']), tbl)
		Y = map(operator.itemgetter(col_dict[col]), tbl)
		ax.plot(X, Y, '--.', label=lbl)
	#
	ax2 = ax.twinx()
	ax2.set_yscale('linear')
	ax2.plot(ts_single, mags, 'o')
	#
	ax.set_xlabel('time (year)', size=14)
	#
	ax = myfig.axes[1]
	ax.plot(ts, normal_stress, '.-', label='normal stress')
	for tbl, col, lbl in [(upper_normal_init, 'event_normal_init', 'init_upper'), (upper_normal_final, 'event_normal_final', 'final_upper'), (lower_normal_init, 'event_normal_init', 'init_lower'), (lower_normal_final, 'event_normal_final', 'final_lower')]:
		X = map(operator.itemgetter(col_dict['event_year']), tbl)
		Y = map(operator.itemgetter(col_dict[col]), tbl)
		ax.plot(X, Y, '--.', label=lbl)
	#
	for ax in myfig.axes:
		ax.legend(loc=0, numpoints=1)
	#	
#
# helper functions:
#
def tmp_test():
	A=[[1,2], [2,3], [4,5], [6,7]]
	P=mpp.Pool(mpp.cpu_count())
	X=P.map_async(map_helper, [[operator.eq]+a for a in A])
	#
	return X.get()


def null_funct(args=[], kwargs={}):
	pass
#
def map_helper(args_in = [null_funct, [], {}]):
	# helper function for pool.map_async(). pass data as a list(-like object):
	# [function, [args], {kwargs}] (though we'll allow for some mistakes).
	#
	funct = args_in[0]
	#
	# allow for different formatting options:
	if not (isinstance(args_in[1], list) or isinstance(args_in[1], tuple) or isinstance(args_in[1], dict)):
		# probably passed a list of parameters. just use them:
		args = args_in[1:]
		#
		return funct(*args)
	#
	# if the args are "properly" formatted:
	args=[]
	kwargs = {}
	for arg in args_in[1:]:
		if isinstance(arg, list) or isinstance(arg, tuple): args += arg
		if isinstance(arg, dict): kwargs.update(arg)
	return funct(*args, **kwargs)
#
def get_peaks(data_in=[], col=0, peak_type='upper'):
	peaks_out = []
	for i, rw in enumerate(data_in[1:-1]):
		if hasattr(rw, '__iter__')==False and hasattr(rw, '__len__')==False:
			#print "making iterable...", rw
			rw=[rw]
			col=0	
		#
		#print data_in[i-1][col], data_in[i][col], data_in[i+1][col], data_in[i-1][col]<data_in[i][col], data_in[i+1][col]<data_in[i][col]
		if peak_type == 'upper':
			if data_in[i][col]>data_in[i-1][col] and data_in[i][col]>data_in[i+1][col]: peaks_out += [data_in[i]]
		if peak_type == 'lower':
			if data_in[i][col]<data_in[i-1][col] and data_in[i][col]<data_in[i+1][col]: peaks_out += [data_in[i]]
		#
	#
	return peaks_out
#
def fetch_data_mpp(n_cpus=None, src_data=[], col_name='', matching_vals=[], is_sorted=False):
	#
	# (there's probably a simpler way to do this usnig pool.apply_async() or something...
	# wrap up the whole multi-processing fetch data process. aka, use this to to a multi-processing search for values in a table.
	# n_cpus: number of processors to use. pass None to get all cpus.
	# src_data: data to be searched.
	# col_name: column in src_data to be searched
	# matching_vals: values to match.
	#
	# basically, it looks like this should aways be used. there appears to be some latency that Process() (or something in the mpp framework)
	# overcomes, so even for a single CPU, this runs WAY faster than a straight list comprehension statement... maybe because it runs compiled?
	#
	if n_cpus == None: n_cpus = mpp.cpu_count()
	#
	my_pipes = []
	output = []
	my_processes = []
	sublist_len = min(len(src_data), (len(src_data)/n_cpus)+1)		# a little trick to ensure we get the full list and don't overrun...
	#
	for i in xrange(n_cpus):
		child_pipe, parent_pipe = mpp.Pipe()
		my_pipes+= [parent_pipe]
		#my_processes+=[mpp.Process(target=find_in, args=(src_data[i*sublist_len:(i+1)*sublist_len], col_name, event_ids_ary, child_pipe))]
		my_processes+=[mpp.Process(target=find_in, kwargs={'tbl_in':src_data[i*sublist_len:(i+1)*sublist_len], 'col_name':col_name, 'in_list':matching_vals, 'pipe_out':child_pipe, 'is_sorted':is_sorted})]
		my_processes[-1].start()
	#
	#print "%d/%d processes started." % (len(my_pipes), len(my_processes))
	#
	for i, proc in enumerate(my_processes): my_processes[i].join()
	#
	# processes are now started and joined, and we'll wait until they are all finished (each join() will hold the calling
	# process until it's finished)
	for i, pipe in enumerate(my_pipes):
		#in_my_pipe = pipe.recv()
		#
		output += pipe.recv()
		my_pipes[i].close()
	#
	return output

def find_in(tbl_in=[], col_name='', in_list=[], pipe_out=None, is_sorted=False):
	# a "find" utility function. look for qualifying rows (presuming at this point a hdf5 type layoud). return the
	# out_list via "return" or a pipe in the case of mpp.
	# (see notes below regarding @is_sorted)
	#
	if is_sorted != True: 
		output = [x for x in tbl_in if x[col_name] in in_list]
	if is_sorted == True:
		# if the data are sorted, we can use an indexing of sorts to improve performance... but we break list-comprehension syntax,
		# and we have to find the max/min values, so in the end this is much more costly and should not be used except perhaps in cases
		# where the soruce data are HUGE and the "in_list" data are quite small. we might also (see below) add the condition that the
		# in_list data are sorted so we can quickly extract the min/max values.
		min_val = min(in_list)
		max_val = max(in_list)
		output = []
		for i, rw in enumerate(tbl_in):
			if rw[col_name]<min_val: continue
			if rw[col_name]>max_val: break
			#
			if rw[col_name] in in_list: output+=[rw]
		#
		# or a variation on the list comprehension?
		#output = [x for x in tbl_in if x[col_name]>=min_val and x[col_name]<=max_val and x[col_name] in in_list]
		#
		#it seems that the sorted method is consistently slower, probably  because it requreis first finding the min/max values (which
		# involves spinning the list. if the in_list were sorted, we could get something out of this by pulling mthe max/min
		# vals directly from that list. also, the list-comprehension syntax does not appear to avoid spinning the whole list, so the
		# standard "for" syntax is preferable.
			
	#
	if pipe_out!=None:
		pipe_out.send(output)
		pipe_out.close()
	else:
		return output
#		
def fetch_h5_data(sim_file=allcal_full_mks, n_cpus=None, table_name=None, col_name=None, matching_vals=[]):
	with h5py.File(sim_file, 'r') as h5_data:
		this_table = h5_data[table_name]
		table_cols = cols_from_h5_dict(this_table.id.dtype.fields)
		#
		output_data = fetch_data_mpp(n_cpus=n_cpus, src_data=this_table, col_name=col_name, matching_vals=matching_vals)
	#
	return output_data
#		
def fetch_h5_table_data(h5_table=None, n_cpus=None, col_name=None, matching_vals=[]):
	# fetch hdf5 data, but with the hdf5 file context specified in the calling function, aka -- skip to the table.
	# (see also fetch_h5_data() for context-not-specified).
	#
	#table_cols = cols_from_h5_dict(h5_table.id.dtype.fields)
	#
	output_data = fetch_data_mpp(n_cpus=n_cpus, src_data=h5_table, col_name=col_name, matching_vals=matching_vals)
	#
	return output_data
#
def cols_from_h5_dict(cols_dict):
	# return column mames from an h5 columns dictionary.
	# dict(proxy) will be like: {str{col_name}:(dtype, width_index),...}
	# SO, we want to return them in order o width_index
	#
	col_names=[]
	for key in cols_dict.keys():
		col_names += [[key, cols_dict[key][1]]]
	col_names.sort(key=lambda x:x[1])
	#
	return [x[0] for x in col_names]

def _worker_in_table(src_table, test_table, test_col, target_Q):
	#
	# note we must be passed a 1-D vector test_table. 
	target_Q.put([x for x in src_table if x[test_col] in test_table])
	#return [x for x in src_table if src_table[test_col] in test_table]
		
	
	
			

