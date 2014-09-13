import matplotlib.pyplot as plt
plt.ion()
#
import math
import h5py
import numpy
import scipy
import operator
#
import cStringIO
import sys
#
import imp
import inspect
import multiprocessing as mpp

#pyvc = imp.load_source('pyvc', '../PyVC/pyvc')
import pyvc as pyvc
pvca = imp.load_source('pyvc.vcanalysis', '../PyVC/pyvc/vcanalysis.py')
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

def get_event_blocks(sim_file=allcal_full_mks, event_number=None, block_table_name='block_info_table'):
	'''
	# get blocks involved in a particular event. first, get sections, then blocks from that section.
	# this will probably all change for the better. for now, use some available tools from PyVC
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
	

def get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None):
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
				my_processes+=[mpp.Process(target=find_in, kwargs={'tbl_in':src_data[i*sublist_len:(i+1)*sublist_len], 'col_name':col_name, 'in_list':event_ids_ary, 'pipe_out':child_pipe})]
				my_processes[-1].start()
			#
			print "%d/%d processes started." % (len(my_pipes), len(my_processes))
			#
			for i, proc in enumerate(my_processes): my_processes[i].join()
			#
			# processes are now started and joined, and we'll wait until they are all finished (each join() will hold the calling
			# process until it's finished)
			for i, pipe in enumerate(my_pipes):
				in_my_pipe = pipe.recv()
				print "adding %d new row..." % len(in_my_pipe)
				section_events += in_my_pipe
				my_pipes[i].close()
					
	#
	#return event_ids_ary
	# columns we might care about:
	return [col_names] + section_events
#
def get_CFF_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, fignum=0):
	# CFF = shear_stress - mu*normal_stress
	# 1) get events for the section (use "events_by_section")
	# 2) for each event, gather all the blocks (use the event_sweep_table).
	# 3) for each block-sweep in the event, add the local CFF (aka, add the shear_init + mu*normal_init for each entry in the event.
	# 4) time (year) of the event can be taken from the event_table.
	#
	CFF=[]
	events = get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=section_id, n_cpus=n_cpus)
	col_dict = {}
	map(col_dict.__setitem__, events[0], range(len(events[0])))
	#
	# so, events is a time series of events (earthquakes) on a particular fault section.
	# for each event, the normal/shear stresses are calculated.
	#
	# now, fetch block level data. for now, let's return a simpler data set than we initially read.
	#
	for event in events:
		# get_event_blocks(sim_file=allcal_full_mks, event_number=None, block_table_name='block_info_table')
		blocks = get_event_blocks(sim_file=sim_file, event_number=event[col_dict['event_number']])
		
#
	
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
def fetch_data_mpp(n_cpus=None, src_data=[], col_name='', matching_vals=[]):
	#
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
		my_processes+=[mpp.Process(target=find_in, kwargs={'tbl_in':src_data[i*sublist_len:(i+1)*sublist_len], 'col_name':col_name, 'in_list':matching_vals, 'pipe_out':child_pipe})]
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
		#print "adding %d new row..." % len(in_my_pipe)
		output += pipe.recv()
		my_pipes[i].close()
	#
	return output

def find_in(tbl_in=[], col_name='', in_list=[], pipe_out=None):
	# a "find" utility function. look for qualifying rows (presuming at this point a hdf5 type layoud). return the
	# out_list via "return" or a pipe in the case of mpp.
	#
	output = [x for x in tbl_in if x[col_name] in in_list]
	#
	if pipe_out!=None:
		pipe_out.send(output)
		pipe_out.close()
	else:
		return output


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
		
	
	
			

