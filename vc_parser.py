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
import json
import cPickle
import time
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
def calc_final_CFF_from_block_data(block_data):
	return block_data['shear_final'] - block_data['mu']*block_data['normal_final']
#
def calc_CFF_from_block_data(block_data, initial_stress=True):
	# a row of block data (assume dict. or hdf5 format):
	return block_data['shear_init'] - block_data['mu']*block_data['normal_init']
#
def calc_CFF(shear_stress=0., normal_stress=0., mu=.8):
	return shear_stress - normal_stress*mu
#
def get_event_CFF(sweep_blocks, n_cpus=None):
	# @sweep_blocks: blocks from event_sweep_table (for a given event(s))
	#
	# CFF = shear_stress - mu*normal_stress
	if n_cpus==None: n_cpus=mpp.cpu_count()
	#
	my_pool = mpp.Pool(processes=n_cpus)
	CFFs = my_pool.map_async(calc_CFF_from_block_data, sweep_blocks)	# imap returns an ordered iterable. map_async()
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
def get_final_event_CFF(sweep_blocks, n_cpus=None):
	# aka, calc the "final", not the "initial" cff for an event (given the set of blocks for
	# that event).
	# @sweep_blocks: blocks from event_sweep_table (for a given event(s))
	# (see get_event_CFF() for additional notes)
	#
	# CFF = shear_stress - mu*normal_stress
	if n_cpus==None: n_cpus=mpp.cpu_count()
	#
	my_pool = mpp.Pool(processes=n_cpus)
	CFFs = my_pool.map_async(calc_final_CFF_from_block_data, sweep_blocks)	
	my_pool.close()	
	my_pool.join()												
	#
	return sum(CFFs.get())
#
def plot_CFF_ary(ary_in='data/VC_CFF_section_125.ary', fnum=0):
	CFF = numpy.load(ary_in)
	#
	# assume either [event_id/num, year, CFF] or [event_id/num, CFF]
	if len(CFF[0])==2:
		y_col=1
	if len(CFF[0])>=3:
		y_col=2
	CFF_peaks = get_peaks(data_in=CFF, col=y_col, peak_type='lower')
	#
	zCFF = zip(*CFF)
	X=zCFF[0]
	Y=zCFF[y_col]
	#
	zCFF_peaks = zip(*CFF_peaks)
	X_peaks = zCFF_peaks[0]
	Y_peaks = zCFF_peaks[y_col]
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
	#plt.semilogy(X, [-1*y for y in Y], '.-')
	plt.semilogy(X_peaks, [-1*y for y in Y_peaks], '.-')
	
#
def get_EMC_CFF():
	sections = [16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 69, 70, 73, 83, 84, 92, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 123, 124, 125, 126, 149]
	#
	for section in sections:
		X=get_CFF_on_section(section_id=section)
		#X=numpy.array(X)
		#X.dump('data/EMC_CFF_timeseries_section_%d.npy' % section)
		with open('data/EMC_CFF_timeseries_section_%d.npy' % section, 'w') as f_out:
			json.dump(X, f_out, indent=1)
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
	CFF=[]
	events = get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=section_id, n_cpus=n_cpus)
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
		# pre-calc a dictionary of event-sweeps. each event_number will include a list of sweep events.
		if event_sweeps_dict==None: event_sweeps_dict = make_h5_indexed_dict_spp(h5_in=sweep_data, index_col='event_number')
		#
		for ev_num, event in enumerate(events[1:]):
			# get_event_blocks(sim_file=allcal_full_mks, event_number=None, block_table_name='block_info_table')
			#
			#event_id=event[col_dict['event_number']]
			event_id=event['event_number']
			event_year = event['event_year']
			# this needs to be sped up maybe -- perhaps by maintaining the hdf5 context?
			#
			#blocks = fetch_h5_data(sim_file=sim_file, n_cpus=n_cpus, table_name='event_sweep_table', col_name='event_number', matching_vals=[event_id])
			#
			# this should be faster, since we don't have to open and reopen the file... but i'm not sure how much
			# faster it is. certainly for smaller jobs, the sipler syntax is probably fine.
			# ... but still slow, so use a pre-indexed list of sweep events...
			#blocks = fetch_h5_table_data(h5_table=sweep_data, n_cpus=None, col_name='event_number', matching_vals=[event_id])
			blocks = event_sweeps_dict[event_id]
			#
			if ev_num%250==0: print "(%d) blocks (%d) fetched: %d" % (ev_num, event_id, len(blocks))
			#
			# ... and now get the block info for each block in the event (from the blocks_data dict).
			# ... but more specifically, let's calculate the cumulative CFF for this event (summing the local CFF for each block).
			#
			#CFF+=[[event_id, event_year, get_event_CFF(blocks), get_final_event_CFF(blocks)]]
			#
			CFF+=[{'event_id':event_id, 'event_year':event_year,'event_magnitude':event['event_magnitude'], 'CFF':get_event_CFF(blocks), 'CFF_final':get_final_event_CFF(blocks)}]
			#print "CFF calculated: %s" % str(CFF[-1])
			#if len(blocks)>10: return blocks
			blocks=None
		#
	return CFF
#
def h5_index_event_number(h5_in=None, index_col='event_number'):
	# ... though we might be able to use initialize() pool kwarg. for now...
	return make_h5_indexed_dict_spp(h5_in=h5_in, index_col=index_col)
#
def make_h5_indexed_dict_spp(h5_in=None, index_col=None, pipe_out=None):
	# worker for make_h5_indexed_dict (or can be used standalone).
	# group rows by unique value(s) in index_col
	#
	# guessing a bit about how to properly optimize this...
	#return_dict = {}
	#key_vals = set(h5_in[index_col])
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
def make_h5_indexed_pool(h5_table=None, n_cpus=None):
	# incidentals:
	index_col='event_number'
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
	pool =  mpp.Pool(processes=n_cpus)
	#
	for i in xrange(n_cpus):
		jobs += [h5_table[job_indices[i]:job_indices[i+1]]]
	#
	results = pool.map_async(h5_index_event_number, jobs)
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
		print "equalities: ", sweeps1==sweeps2, sweeps2==sweeps3, sweeps3==sweeps1, sweeps1==sweeps4, sweeps2==sweeps4, sweeps3==sweeps4
		#
	return None
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
		#print "adding %d new row..." % len(in_my_pipe)
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
		
	
	
			

