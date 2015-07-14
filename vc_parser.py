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
#
import numpy
import scipy
import scipy.optimize as spo
from scipy.spatial import KDTree
from scipy.spatial import cKDTree
#
import operator
import glob
import random
import datetime as dtm
import pytz
#
import cStringIO
import sys
import json
import cPickle
import time
import os
import datetime as dtm
import pytz
#
import imp
import inspect
import multiprocessing as mpp
#
import ANSStools
import pca_tools
try:
	import BASScast

except:
	pass

try:
	import quakelib
	import pyvc
	import pyvc.vcanalysis as pvca
except:
	print "exception loading quakelib. quakelib functions not available."

#pyvc = imp.load_source('pyvc', '../PyVC/pyvc')
#import pyvc as pyvc
#pvca = imp.load_source('pyvc.vcanalysis', '../PyVC/pyvc/vcanalysis.py')
#import pyvc
#import pyvc.vcanalysis as pvca
#

# sections filter for EMC related queries. (should be the set of fault sections used in the socal/EMC simulations). we can use this filter to mine the full
# AllCal data.
napa_region_section_filter = {'filter':set([45, 50, 172, 171, 170, 169, 44, 168, 167, 139, 40, 142, 41, 46])}
napa_sections = list(napa_region_section_filter['filter'])

emc_event = {'lat': 32.128, 'lon':-115.303, 'mag':7.2, 'event_time':dtm.datetime(2010,4,4,3,40,41,tzinfo=pytz.timezone('US/Pacific'))}
emc_section_filter = {'filter': (16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 69, 70, 73, 83, 84, 92, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 123, 124, 125, 126, 149)}
allcal_full_mks = '../ALLCAL2_1-7-11_no-creep_dyn-05_st-20.h5'
default_sim_file = allcal_full_mks
emc_sections = list(emc_section_filter['filter'])

# all sections:
with h5py.File(default_sim_file, 'r') as F:
	all_section_ids = list(set(F['block_info_table']['section_id'][()]))

class getter(object):
	# a simple container class to emulate objects that return values. for example,
	# emulate Random.random() when a known value is always returned.
	def __init__(self, rand_val=0.0, get_val=0.0):
		self.rand_val = rand_val
		self.get_val = get_val
	#
	def get(self, val=0.0):
		if val==None: val=self.get_val
		#
		return val
	#
	def random(self, val=0.0):
		if val==None: val=self.rand_val
		#
		return val
		

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
def get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, is_sorted=False, recarray=False):
	#
	# TODO: there is a little bit of an inconsistenty. when we make the pre-calc data sets, we include some derived fields, namely CFF calcs.
	# so for now, let's include these in the recarray, but to avoid compatibility problems, not in the list-type returns (though since we'd be
	# tacking onto the ends of each row, we could probably get away with doing that... and anyway, we'll want to move in the direction of
	# including those information -- aka, so the dynamic/pre-calc returns are as identical as possible.
	#
	# note in comments and code that this is optimized using the h5py.File()[{index}][()]  (specifically, the
	# [()] ) syntax ), which immediately copies the hdf5 content to a numpy array. note also that this is not
	# always a good idea for large data sets, or when we only want to fetch a small piece of a data set
	# (see inline comments).
	#
	# note on multiprocessing: the basic trend seems to be that 1) there is overhead in using mpp (aka, use mpp
	# for 1 cpu; it takes a bit longer than the "raw" SPP process (for overhead, pickling, etc.).
	# BUT, all and all, multiprocessing seems to be faster, even for just 2 processors (though perhaps
	# marginally so for a single section).
	#
	if n_cpus == None: n_cpus = mpp.cpu_count()
	#
	with h5py.File(sim_file,'r') as vc_data:
		#
		#
		# note: vc_data['events_by_section'] is a group, not an array. each entry is an array, so we can't
		# just use [()] syntax. but we can make a dict. of the entries. if we wanted to:
		# ebs_index = {key:tbl for key,tbl in ebs.iteritems()}
		# which returns immediately and returns a dict. of pointers/iterators to the hdf5 entries. this would be,
		# probably, fastest for a scenario in which we want to return a few arrays from a large data set, particularly
		# if data volume is an issue. note, however, that the hdf5 File object must remain open.
		# alternatively:
		# ebs_index = {key:tbl[()] for key,tbl in ebs.iteritems()}
		# returns a dict with the tables automatically loaded into nd arrays, so type(key)=string, type(value)=ndarray
		# this syntax loads a bit more slowly but is then faster if we want to look at many or all of the arrays
		# note also that we can close the h5py.File() object, and of course that this approach is memory intensive.
		view_ebs = vc_data['events_by_section']		# "events by section" group
		#
		#event_ids_tbl = view_ebs['section_%d' % section_id]
		event_ids_ary = view_ebs['section_%d' % section_id][()]
		#
		#event_ids_ary = numpy.array(numpy.zeros(event_ids_tbl.len()), dtype='int64')	# array needs to be initialized with data type for hdf5.read_direct()
		##event_ids_ary.shape(numpy.zeros(event_ids_tbl.len(), len(event_ids_tbl[0])))
		#event_ids_tbl.read_direct(event_ids_ary)		# so we initialize the array, then read_direct() into it. can we insteadn just x=XX[()] ??
		
		#print "types: ", type(event_ids_ary2), type(event_ids_tbl), type(event_ids_ary)
		#print "dtypes: ", event_ids_ary2.dtype, event_ids_tbl.dtype, event_ids_ary.dtype
		
		#
		# now, fetch event data:
		# this is a choice, computationally. we can do it in a list comprehension [], or we can index the source list (event_id values
		# are sorted, so once we find one event in the master list, we can stop looking for it in the event_ids_ary list.
		# list comprehension is simpler...
		src_data = vc_data['event_table'][()]
		#col_names = cols_from_h5_dict(src_data.id.dtype.fields)		# or this can be done from .dtype.names, .dtype.formats (or whatever that is).
		col_names = src_data.dtype.names
		col_name = 'event_number'
		#if n_cpus in (1, None):
		if n_cpus == None or n_cpus==1:
			# this should never be true any longer. see note above.
			print "running in SPP mode"
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
	# add 'cff_initial', 'cff_final' columns. try for both list and recarray return. if backwards compatiblity is a problem, well, handle it...??
	#
	if not recarray:
		return [col_names] + section_events
	else:
		return numpy.core.records.fromarrays(zip(*section_events), names=col_names, formats=[type(x).__name__ for x in section_events[0]])
#
def get_blocks_info_dict(sim_file=allcal_full_mks, block_ids=None):
	# if block_ids==None, get all of them...
	# and we could parallelize this, or we can use it as a single processor implementation that can be run in parallel...
	#
	if isinstance(block_ids,float) or isinstance(block_ids,int):
		block_ids = [int(block_ids)]
	#
	with h5py.File(sim_file, 'r') as vc_data:
		blocks_table = vc_data['block_info_table'][()]		# "blocks" should be pretty small, and this runs about
															# 10,000 times faster (literally maybe... at least) with [()]

		#block_cols = vc_data['block_info_table'].dtype.names
		#
		'''
		blockses = {}
		for block in blocks_table:
			if block_ids!=None and block['block_id'] not in block_ids: continue
			#
			#block_dict = {}
			block_dict = {key:block[key] for key in block_cols}
			#
			#[block_dict.__setitem__(key, block[key]) for key in block_cols]
			blockses[block['block_id']] = block_dict
		'''
		#
		#blockses = {block['block_id']:{key:block[key] for key in block_cols} for block in blocks_table if not (block_ids!=None and block['block_id'] not in block_ids)}
		#blockses = {block['block_id']:{key:block[key] for key in vc_data['block_info_table'].dtype.names} for block in blocks_table if not (block_ids!=None and block['block_id'] not in block_ids)}
	#
	return {block['block_id']:{key:block[key] for key in blocks_table.dtype.names} for block in blocks_table if not (block_ids!=None and block['block_id'] not in block_ids)}
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
	total_cff_init = sum(sweep_blocks['shear_init'] - sweep_blocks['mu']*sweep_blocks['normal_init'])
	total_cff_final = sum(sweep_blocks['shear_final'] - sweep_blocks['mu']*sweep_blocks['normal_final'])
	
	#
	'''
	for blk in sweep_blocks:
		total_cff_init  += blk['shear_init']  - blk['mu']*blk['normal_init']
		total_cff_final += blk['shear_final'] - blk['mu']*blk['normal_final']
	if total_cff_init2!=total_cff_init: print "total_init error"
	if total_cff_final2!=total_cff_final: print "total final error"
	'''
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
def combine_section_CFFs(sections=[], ary_in_format='data/VC_CFF_timeseries_section_%d.npy', start_year=0., end_year=None, sim_file=default_sim_file, create_missing_catalogs=True):
	# concatenate (then sort by year) a set of section catalogs. assume proper recarray() format, etc.
	# ... and this should definitely be used for multi-fault combined catalogs; there are ~50 duplicate events in
	# the full EMC set.
	# create_missing_catalogs: true: if the pre-processed array does not exist, then create it. this can be time consuming,
	# so it may be desirable to set this to False and let the script fail.
	#
	# handle some input variations:
	if sections==None: return sections
	if not hasattr(sections, '__len__'): sections = [sections]
	if len(sections)==0:
		if sim_file==None: return sections
		#
		# try to get all the sections. otherwise, return None/sections in defeat:
		try:
			with h5py.File(sim_file) as vc_data:
				sections = list(set(vc_data['block_info_table']['section_id']))
		except:
			return sections
	#
	# let's just do this here... check all section_ids to see that their pre-processed catalogs have been made. if not, make them.
	for i,sec_id in enumerate(sections):
			# get_CFF_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, event_sweeps_dict=None)
			# get_EMC_CFF(sections=None, file_out_root='data/VC_CFF_timeseries_section')
			ary_file_name = ary_in_format % sec_id
			if ary_file_name not in glob.glob(ary_file_name):
				print "creating section catalog: %s" % ary_file_name
				#ary_data = get_CFF_on_section(sim_file=default_sim_file, section_id=sec_id)	
				ary_data = get_EMC_CFF(sections=[sec_id], file_out_root=ary_in_format.split('_%d')[0])
	
	#
	if len(sections)>=1:
		combined_catalog = numpy.load(ary_in_format % sections[0])
	#
	# ... and let's be sure we're not getting duplicate events as per section partitioning...
	for i,sec_id in enumerate(sections[1:]):
		# (or we can load and append in one fell swoop, but this should be fine).
		try:
			#
			combined_catalog = numpy.append(combined_catalog, numpy.load(ary_in_format % sec_id))
		except:
			print "failed to add sub-catalog for %d:: %s, %d" % (sec_id, ary_in_format % sec_id, len(combined_catalog))
	#
	#l_cat=[rw.tolist() for rw in combined_catalog if rw['event_year']>=start_year]
	#print "type: ", type(l_cat)
	# only unique rows (re-create the recarray from a list-ofthe-set-ofthe-list-ofthe-recarray):
	#combined_catalog=numpy.rec.array(list(set(combined_catalog.tolist())), dtype=combined_catalog.dtype)
	#
	end_year = (end_year or max(combined_catalog['event_year']))
	#combined_catalog=numpy.rec.array(list(set([rw.tolist() for rw in combined_catalog if rw['event_year']>=start_year])), dtype=combined_catalog.dtype)
	combined_catalog=numpy.rec.array(list(set([rw.tolist() for rw in combined_catalog if (rw['event_year']>=start_year and rw['event_year']<=end_year)])), dtype=combined_catalog.dtype)
	
	combined_catalog.sort(order='event_year')
	#
	return combined_catalog
#
##########################################
# ROC and predictability bits.
#
def get_random_forecast_set(CFF=None, section_id=16, m0=7.0, P_min=0., P_max=1.0, set_name=None, nits=1000, format='recarray', do_roc_plot=False, fignum=0):
	'''
	# get a set of random forecasts, for comparison to real forecasts. return as list of dicts.
	'''
	R=random.Random()
	P_min = max(0.0, P_min)
	P_max = min(1.0, P_max)
	#
	
	if (CFF==None or len(CFF)==0) and section_id!=None:
		CFF = combine_section_CFFs(sections=section_id)
	if isinstance(CFF, str):
		CFF = numpy.load(ary_in)
	#
	roc_prams = []
	for j in xrange(nits):
		this_P = P_min + (P_max-P_min)*R.random()
		print "this_P: ", this_P
		fc = forecast_random_1(ary_in=CFF, m0=m0, do_plot=False, set_name=set_name, section_id=None, P=this_P)
		fc_eval = evaluate_alert_segments(alert_segments=fc, CFF=CFF, do_plot=False)
		fc_eval['P']=this_P
		roc_prams+=[fc_eval]
	#
	return_object = roc_prams
	#
	if format in ('array', 'rec', 'recarray') or do_roc_plot:
		lst_prams = [[val for val in rw.itervalues()] for rw in roc_prams]
		col_names = roc_prams[0].keys()
		col_formats = [type(x).__name__ for x in lst_prams[0]]
		#
		roc_prams_ary = numpy.core.records.fromarrays(zip(*lst_prams), names=roc_prams[0].keys(), formats=[type(x).__name__ for x in lst_prams[0]])
		return_object = roc_prams_ary
	#
	if do_roc_plot:
		plt.figure(fignum)
		plt.clf()
		plt.plot(roc_prams_ary['F'], roc_prams_ary['H'], '.', zorder=2)
		plt.plot([0., 1.], [0., 1.], '-', lw=2, zorder=1)
		#
		f=plt.figure(fignum+1)
		plt.clf()
		ax3d = f.add_axes([.1, .1, .8, .8], projection='3d')
		ax3d.plot(roc_prams_ary['P'], roc_prams_ary['F'], roc_prams_ary['H'], '.')
		ax3d.set_xlabel('Probability $P$')
		ax3d.set_ylabel('False Alarm Rate $F$')
		ax3d.set_zlabel('Hit Rate $H$')
		ax3d.plot([0., 1.], [0., 1.], [0., 1.], '-', lw=2)
		
	#
	return return_object
#
def forecast_random_1(ary_in=None, m0=7.0, do_plot=False, fnum=0, set_name=None, section_id=16, P=.5, R=None):
	'''
	# a random forecast. for now, let's leave the array alone. we'll just use it to get the correct length, etc.
	# assign a random forecast by randomly choosing to turn off/on the forecast alert... i suppose with some probability.
	#
	# R: random number generator. can be None and we'll get one...
	#
	# basically: 
	#  - keep track of alert on/off statis
	#  - at each step:
	#    - if alert is on: turn off with probability Q
	#    - if alert if off: turn on with probability P
	# this should be equivalent to (right?) at each step, alert status is ON with probability P, OFF with probability Q.
	# then, use: evaluate_alert_segments() to get basic hit/false alarm stats.
	# 
	#
	'''
	#
	if P==None: P=.5
	if R==None: R=random.Random()
	print "from fcr1: P=%f" % P
	#Q = 1.0-P
	#
	# handle catalog inputs:
	#
	if (ary_in==None or len(ary_in)==0) and section_id!=None:
		#ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id
		ary_in = combine_section_CFFs(sections=section_id)
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		# otherwise, assume we've been passed a proper CFF object:
		CFF = ary_in
		# give the set a name  so we don't return the whole input data object...
		if set_name==None: set_name='unnamed_CFF'
	#
	if set_name==None:
		# still none? must be a string type ary_in...
		set_name=ary_in
	#
	alert_segments = []
	alert_status=False
	this_b = None		# holder variable to make these alert_segments look more like a real alert metric.
	for i,rw in enumerate(CFF):
		r = R.random()
		#
		#print "diag: ", r, alert_status
		#
		if alert_status==True and r>P:
			# alert is on and we cancel it...
			alert_status=False
			alert_segments[-1]
			alert_segments[-1]+=[[rw['event_year'], this_b]]
			#alert_segments += [[]]
		if alert_status==False and r<=P:
			alert_status=True
			alert_segments +=[ [[rw['event_year'], this_b]] ]
		#
	if alert_status==True:
		alert_segments[-1] += [[rw['event_year'], this_b]]
	while len(alert_segments)>0 and len(alert_segments[-1])<2: alert_segments.pop(-1)
	#
	#return alert_segments
	#
	# plot it up:
	if do_plot:
		plt.figure(0)
		plt.clf()
		for rw in alert_segments:
			plt.plot([rw[0][0], rw[1][0]], [0.,0.], 'o-')
	#
	return alert_segments	
#
def evaluate_alert_segments(alert_segments=None, CFF=None, section_id=None, m0=7.0, do_plot=True, fnum=0):
	# generic alert analyzier:
	# this should be used by all alert/forecast algorithms. basically, produce a forecast as alert segments
	# like [ [t_0_0, t_0_1, t_0_2, t_0_3], [t_1_0, t_1_1, t_1_2], etc. ] so each alert exists for t>t_j_0 and t<=t_j_-1
	# given a set of alert segments and a CFF, evaluate the alerts:
	# ... and let's just do this the hard way; for each earthquake see if its in an alert bin. going forward, we'll want to index this.
	if CFF==None and section_id!=None:
		CFF=combine_section_CFFs(sections=section_id)
	#
	total_alert_time = 0.0
	total_total_time = CFF[-1]['event_year'] - CFF[0]['event_year']
	#
	# first, did we get any segments? it's possible (especially for a short catalog), to return zero segments.
	if len(alert_segments)==0:
		return {'total_alert_time': 0., 'total_time':0., 'n_predicted':0., 'n_missed':0., 'H':None, 'F':None}
	#
	# initial versions of this use an alert_segments format like [ [ [t_1, alert_val_1], [t2, alert_val_2] ]...]
	# for most purposes, we want to simplify this to : [ [t1, t2], ...]
	if hasattr(alert_segments[0][0], '__len__'):
		# we have the original list-list format.
		for i,rw in enumerate(alert_segments):
			alert_segments[i][0] = alert_segments[i][0][0]
			alert_segments[i][1] = alert_segments[i][1][0]
	#
	for alert_segment in alert_segments:
		#total_alert_time += (alert_segment[-1][0] - alert_segment[0][0])
		total_alert_time += (alert_segment[-1] - alert_segment[0])
	#
	# and prediction success:
	n_predicted = 0
	n_missed = 0
	n_total = 0
	#
	j_alert_start = 0
	for i, rw in enumerate(CFF):
		if rw['event_magnitude']<m0: continue
		n_total+=1
		#
		#
		#while rw['event_year']>alert_segments[j_alert_start][1][0] and j_alert_start<(len(alert_segments)-1): j_alert_start+=1
		while rw['event_year']>alert_segments[j_alert_start][-1] and j_alert_start<(len(alert_segments)-1): j_alert_start+=1
		#
		for alert in alert_segments[j_alert_start:]:
			#if rw['event_year'] > alert[0][0] and rw['event_year'] <= alert[1][0]:
			if rw['event_year'] > alert[0] and rw['event_year'] <= alert[-1]:
				# note: using alert[-1] instead of alert[1] allows for longer alert segments, aka, [start, middle, middle, mittle, end]
				# instead of just [start, end]
				n_predicted +=1
				#n_missed -=1
				break

	n_missed = n_total-n_predicted
	
	#
	#
	if do_plot:
		# diagnostic plots of forecast metric:
		plt.figure(fnum)
		plt.clf()
		plt.title('General Forecast Metric')
		#
		ax_mags = plt.gca()
		ax_metric = ax_mags.twinx()		
		#
		ax_ints = ax_metric.twinx()
		ax_ints.set_yscale('log')
		ax_metric.set_yscale('linear')
		#
		min_mag = min(CFF['event_magnitude'])
		max_mag = max(CFF['event_magnitude'])
		intervals = [x-CFF['event_year'][i] for i,x in enumerate(CFF['event_year'][1:])]
		#
		#metric_pad_factor = min_mag
		#metric_pad_factor = 0.
		#min_metric = 0.
		#min_metric = alert_segments[0][0][1]
		#max_metric = alert_segments[0][0][1]
		min_metric=0.
		max_metric=1.
		#
		# do a quick spin to get min/max values and other useful stats:
		#for segment in alert_segments:
		#	X,Y = zip(*segment)
		#	min_metric = min(min(Y), min_metric)
		#	max_matric = max(max(Y), max_metric)
		#
		min_x = None
		max_x = None
		for segment in alert_segments:
			#X,Y = zip(*segment)
			#min_metric = min(min(Y), min_metric)
			X = segment
			Y = [1. for x in X]
			#
			# ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b']  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='m', zorder=1, alpha=.5)
			#ax_trend2.fill_between([x['event_year'] for x in trend_data], [1.  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0.0 for x in trend_data], color='m', zorder=1, alpha=.25)
			#
			# show the metric value:
			#ln_metric = ax_metric.fill_between(X,[y for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [y<0. for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			#ln_metric = ax_metric.fill_between(X,[y for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [f_gt_lt(y,0.) for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			# later, add the option for 'detail' metrics (aka, [ [[t,val], [t,val], ...] ]instead of [ [t,t,t,t], [t,t,t,] ]
			#
			#ax_metric.plot([x for i,x in enumerate(X) if Y[i]<0] , [0.0 for y in Y if y<0], 'm-')
			#ax_metric.plot([x for i,x in enumerate(X) if Y[i]<0], [y for y in Y if y<0], 'm-')
			
			ax_metric.plot(X,[0. for x in X], 'mo-')
			ax_mags.fill_between(X, [max_mag for x in X], [min_mag for x in X], color='m', alpha=.2)
			#
			min_x = min((min_x or X[0]), min(X))
			max_x = max((max_x or X[0]), max(X))
			
			#
			# and just an "alert!" box:
			#ln_metric = ax_metric.fill_between(X,[min_metric for y in Y],y2=[0.0 for y in Y], color='m', alpha=.15, where = [y<0. for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			#ln_metric = ax_metric.fill_between(X,[min_metric for y in Y],y2=[0.0 for y in Y], color='m', alpha=.15, where = [f_gt_lt(y,0.) for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			#
			#ax_mags.fill_between(X, [min_mag for x in X], [m0 for x in X], zorder=5, alpha=.2, color='c')
		#
		ax_metric.plot([min_x, max_x], [0., 0.], 'k-')
		#
		ln_ints = ax_ints.plot(CFF['event_year'][1:], intervals, 'b.-', lw=2, ms=7, label='Intervals $\Delta t = t_i - t_{i-1}$')
		#ax_metric.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		ln_mags = ax_mags.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='g', alpha=.7, lw=1.5, label='magnitudes')
		X_big_mags, Y_big_mags = zip(*[[x['event_year'], x['event_magnitude']] for x in CFF if x['event_magnitude']>m0])
		ax_mags.vlines(X_big_mags, [min_mag for x in Y_big_mags], Y_big_mags, color='r', alpha=.9, lw=2.75)
		#
		# cosmetics:
		ax_mags.set_xlabel('Event Year $t$')
		ax_ints.set_ylabel('Inter-event interval $\Delta t$')
		ax_metric.yaxis.set_ticks([])
		ax_mags.set_ylabel('Earthquake magnitude $m$')
		#
		ax_metric.set_ylim(ymin=1.15*min_metric, ymax = -min_metric/10.)
	#
	print "preliminary report:"
	H = float(n_predicted)/(float(n_predicted)+n_missed)
	F = total_alert_time/total_total_time
	#print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, total_alert_time/total_total_time)
	#print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, float(n_predicted)/(float(n_predicted)+n_missed))
	#print "score: H-F: %f " % (float(n_predicted)/(float(n_predicted)+n_missed) - total_alert_time/total_total_time)
	print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, F)
	print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, H)
	print "score: H-F, H/F: %f / %f " % (H-F, H/F)
	#
	return {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'H':H, 'F':F}
	#, 'ary_in_name':set_name, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
	# ... and we might add some sort of returned parameters (from the metric component), but for now, keep it simple.
#
def ROC_single_prams(section_ids=emc_sections, b_0=0., nyquist_factor=.5, m0=7.0, fignum=None):
	'''
	# produce a simple, single parameterization ROC plot (and data set) for a bunch of section_ids.
	'''
	#
	roc_prams={sec_id:{} for sec_id in section_ids}
	#
	for sec_id in section_ids:
		alert_segs=psa_forecast_1(section_id=sec_id, m0=m0, b_0=b_0, nyquist_factor=nyquist_factor, fnum=fignum, do_plot=False, f_gt_lt=operator.lt, detail=False)
		roc_prams[sec_id] = evaluate_alert_segments(alert_segments=alert_segs, CFF=None, section_id=sec_id, m0=m0, do_plot=False, fnum=fignum)
		#
	#
	if fignum!=None:
		plt.figure(fignum)
		plt.clf()
		plt.plot(range(2), range(2), 'r-', lw=2)
		plt.plot([rw['F'] for rw in roc_prams.itervalues()], [rw['H'] for rw in roc_prams.itervalues()], 'o')
		#
	return roc_prams
#
def psa_forecast_1(ary_in=None, m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False, fnum=0, set_name=None, f_gt_lt=operator.lt, section_id=16, detail=False):
	'''
	# a simplified version of forecast_metric_1.
	# output only the alert segments; then use evaluate_alert_segments() to get ROC data.
	# note that, at least in its native form, this outputs a binary alert: alert is on/off; there is no 'value' (aka, slope value b)
	# or degree associated. note that another model for this sort of work is to just output a full 1:1 alert set, aka
	# [ [time, alert_value]
	'''
	#
	# copied (mostly) from the first half of forecast_metric_1() (which can probably be removed at this point).:
	#
	if (ary_in==None or len(ary_in)==0) and section_id!=None:
		#ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id		
		ary_in=combine_section_CFFs(section_id)		# ... and this function knows to use the string format above...
		#											# ... and note we do CFF=ary_in in the if-else clause below...
		#
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		# otherwise, assume we've been passed a proper CFF object:
		CFF = ary_in
		# give the set a name  so we don't return the whole input data object...
		if set_name==None: set_name='unnamed_CFF'
	#
	if set_name==None:
		# still none? must be a string type ary_in...
		set_name=str(ary_in)
	#
	#print "get recurrence...** %d **" % len(CFF)
	recurrence_data = mean_recurrence(ary_in=CFF, m0=m0)
	nyquist_len = max(int(nyquist_factor*recurrence_data['mean_dN_fault']), 2)
	nyquist_time = nyquist_factor*recurrence_data['mean_dT']
	#
	#print "recurrenced. get trend data..."
	#
	trend_data = get_trend_analysis(ary_in=CFF, nyquist_len = nyquist_len, nyquist_time=nyquist_time)
	#trend_data_dict = {trend_data['event_year']:x for x in trend_data}
	#CFF_dict = {x['event_year']:x for x in CFF}
	CFF_dict = {x['event_number']:x for x in CFF}
	
	print "trend lengths: ", len(trend_data), len(CFF), nyquist_len
	max_n = len(trend_data)
	#
	# trend_data includes columns: ['event_number', 'event_year', 'lin_fit_a', 'lin_fit_b', 'rb_ratio', 'interval_rate_ratios']
	#
	# first, just get the total time under alert:
	alert_segments = [[]]		# a collection of lists...
	#alert_segments = []
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
		#
		#if this_b >= b_0:
		if not f_gt_lt(this_b, b_0):		# abstracting the gt/lt option...
			# null case; does not meet "alert" criteria.
			# if we've been collecting "alert" events, stop. if not, just troll along...
			if len(alert_segments[-1])>0:
				#
				alert_segments+= [[]]
		
		else:
			# f_gt_lt(this_b, b_0) is True
			# accelerating case (if f_gt_lt = operator.lt):
			#this_mag = CFF[-max_n:][i+1]['event_magnitude']
			this_mag = CFF_dict[rw['event_number']]['event_magnitude']
			#
			if len(alert_segments[-1])==0:
				# this is the first entry into this alert segment.
				alert_year = trend_data[i]['event_year']
				if detail:
					alert_segments[-1]+=[[alert_year, trend_data[i]['lin_fit_b']]]
				else:
					alert_segments[-1]+=[alert_year]
			#
			# these two cases are additional entries into this alert segment. if it's a small event, extend alert to the
			# next event (aka, in real live, we'd turn "on" the alert and reevaluate at the next event. if it's a big event,
			# it's the one we've been lookng for. we have a rule to turn off an alert immediately after a big event... though we might
			# reevaluate this since we've started to favor the operator.gt() version of this metric.
			if this_mag<m0:
				# add the *next* event as the alert (aka, issue an alert until we have new data).
				#print "len alert_seg: %d" % len(alert_segments[-1])
				#
				# generalize language a bit:
				alert_year = trend_data[i+1]['event_year']
				if detail:
					alert_segments[-1]+=[[alert_year, this_b]]
				else:
					alert_segments[-1]+=[alert_year]
				#pass
				#
			if this_mag>=m0:
				# this is "the" earthquake. add this entry (it's probably already there) from the previous entry.
				#
				alert_year = trend_data[i]['event_year']
				if detail:
					alert_segments[-1]+=[[alert_year, this_b]]
				else:
					alert_segments[-1]+=[alert_year]
			#
		#
	#
	while len(alert_segments)>0 and len(alert_segments[-1])==0: alert_segments.pop()
	#
	return alert_segments
#
def plot_best_interval_roc(n_best=1, fnum=0):
	# parameterize this later...
	my_path = 'interval_metric_output_mpp'
	#my_path = 'interval_metric_output'
	#
	fls = glob.glob('%s/roc_interval_metric_*_LT.npy' % my_path)
	#return fls
	for j in xrange(3):
		plt.figure(fnum+j)
		plt.clf()
	plt.figure(fnum)
	plt.plot(range(2), range(2), 'r-', lw=2.5)
	for fl in fls:
		# calc. the score and somehow return the correct rows. let's just calc. an index...
		ary = numpy.load(fl)
		scores = [[j, rw['H']-rw['F']] for j,rw in enumerate(ary)]
		scores.sort(key=lambda x: x[1])
		#
		best_ary = [ary[rw[0]] for rw in scores[-n_best:]]
		#
		H=[rw['H'] for rw in best_ary]
		F=[rw['F'] for rw in best_ary]
		dt_0s = [rw['dt_0'] for rw in best_ary]
		dt_m0 = best_ary[0]['dt_m0']
		#
		plt.figure(fnum)
		plt.plot(F,H, 'o-', lw=2)
		#
		plt.figure(fnum+1)
		ax1=plt.gca()
		ax1.plot(dt_0s, [h-f for h,f in zip(H,F)], 'o-')
		plt.figure(fnum+2)
		ax2=plt.gca()
		#ax2=ax1.twinx()
		ax2.set_xscale('log')
		ax2.plot([dt/dt_m0 for dt in dt_0s], [h-f for h,f in zip(H,F)], 'o-')
		
		
#
def plot_interval_roc(roc_data=16, fignum_0=0):
	# do some plots for the interval roc metric.
	# roc_data can be the return from interval_forecast_mc(), the output file there-from (a numpy array; we should trap for list-o-dicts.)
	# if roc_data is an integer, we'll infer it's supposed to be a section_id and grab a default file format.
	#
	if isinstance(roc_data, int):
		roc_data = 'interval_metric_output/roc_interval_metric_%d_LT.npy' % roc_data
	if isinstance(roc_data, str):
		roc_data=numpy.load(roc_data)
	#
	max_score = max(roc_data['H']-roc_data['F'])
	scores = roc_data['H']-roc_data['F']
	max_row = [rw for rw in roc_data if rw['H']-rw['F']==max_score][0]
	
	# and we might trap for dict-->ary here (and we should write a short script to do this).
	#
	figs=[]
	# start with a basic ROC:
	figs+=[plt.figure(fignum_0+len(figs))]
	plt.clf()
	ax1 = figs[-1].gca()
	ax2 = ax1.twinx()
	ax1.plot(range(2), range(2), 'r-', lw=2)
	ax1.plot(roc_data['F'], roc_data['H'], 'b.', zorder=3, label='ROC')
	ax1.plot(max_row['F'], max_row['H'], 'r*', ms=15, zorder=5, label='$max(H-F)$')
	ax1.plot(max_row['F'], max_row['H'], 'k*', ms=18, zorder=4)
	ax1.set_xlabel('False alarm rate $F$')
	ax1.set_ylabel('Hit rate $H$')
	ax1.legend(loc='center right', numpoints=1)
	#
	scoreses = [[x,y] for x,y in zip(roc_data['F'], scores)]
	scoreses.sort(key=lambda x: x[0])
	scoreses = zip(*scoreses)
	ax2.plot([x for x in scoreses[0]], scoreses[1], 'g.', alpha=.7, label='$score=H-F$', zorder=2)
	ax2.legend(loc='lower right', numpoints=1)
	#
	figs+=[plt.figure(fignum_0+len(figs))]
	figs[-1].clf()
	ax3d = figs[-1].add_subplot(111, projection='3d')
	ax3d.plot(roc_data['nyquist_factor'], roc_data['dt_0'], roc_data['H']-roc_data['F'], '.')
	ax3d.set_xlabel('nyquist_factor')
	ax3d.set_ylabel('dt_0')
	ax3d.set_zlabel('score')
	#
	if 'dt_m0' in roc_data.dtype.names:
		figs+=[plt.figure(fignum_0+len(figs))]
		figs[-1].clf()
		ax3d = figs[-1].add_subplot(111, projection='3d')
		ax3d.plot(roc_data['nyquist_factor'], roc_data['dt_0']/roc_data['dt_m0'], roc_data['H']-roc_data['F'], '.')
		ax3d.set_xlabel('nyquist_factor')
		ax3d.set_ylabel('dt_0/dt_m0')
		ax3d.set_zlabel('score')
	#
	return scoreses
#
def interval_forecast_mc_set_mpp(section_ids=emc_sections, nits=1000, dt_min=None, dt_max=None, nyquist_min=.2, nyquist_max=1.2, m0=7.0, f_gt_lt=operator.lt, percentile_lower=.1, percentile_upper=.8, output_dir='interval_metric_output', fignum=0, fignum_1=1, n_cpus=None):
	#
	# churn out a bunch of MC rate forecast analyses.
	# mpp version of this (which seems to work pretty well). we should merge this with the spp version, basically split the spp version... but it would be nice, from a huge amounts of memory perspective, to be able to do these one at a time (aka, not have to collect the whole results collection at one time. of course, the tricky part is to handle pyplot (tricky modal state).
	#
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	metric_names={operator.lt:'LT', operator.gt:'GT', operator.le:'Leq', operator.ge:'Geq'}
	#
	roc_maxes = []
	#
	if n_cpus==None: n_cpus=mpp.cpu_count()
	
	mpp_loop_index = n_cpus
	while (mpp_loop_index-n_cpus)<len(section_ids):
		P=mpp.Pool(n_cpus)
		print "***************\n new loop: ", mpp_loop_index
		#
		# 'dump_file':'%s/roc_interval_metric_%d_%s.npy'% (output_dir, sec_id, metric_names.get(f_gt_lt, 'other'))
		secs = section_ids[max(mpp_loop_index-n_cpus,0):min(mpp_loop_index, len(section_ids))]
		print "secs: ", secs
		p_results = [P.apply_async(interval_forecast_mc, kwds={'nits':nits, 'dt_min':dt_min, 'dt_max':dt_max, 'nyquist_min':nyquist_min, 'nyquist_max':nyquist_max, 'CFF':None, 'm0':m0, 'do_plot':False, 'fnum':0, 'set_name':None, 'f_gt_lt':f_gt_lt, 'section_id':sec_id, 'detail':False, 'percentile_lower':percentile_lower, 'percentile_upper':percentile_upper, 'dump_file':None}) for sec_id in  secs]
		P.close()
		P.join()
		#
		#for sec_id in section_ids:
		for sec_id, res in zip(secs, p_results):
		#for sec_id in section_ids:
			roc_datas = res.get()
			#roc_datas = interval_forecast_mc(nits=nits, dt_min=dt_min, dt_max=dt_max, nyquist_min=nyquist_min, nyquist_max=nyquist_max, CFF=None, m0=m0, do_plot=False, fnum=0, set_name=None, f_gt_lt=f_gt_lt, section_id=sec_id, detail=False, percentile_lower=percentile_lower, percentile_upper=percentile_upper)
			#
			#roc_datas = numpy.load('%s/roc_interval_metric_%d_%s.npy'% (output_dir, sec_id, metric_names.get(f_gt_lt, 'other')))
			# note, we can also (nominally) use: plot_interval_roc(roc_data=16, fignum_0=0) and then do the save bit on the first fignum.
			#
			# might not get any data back:
			if roc_datas==None or len(roc_datas)==0: continue
			#
			plt.figure(fignum)
			plt.clf()
			plt.plot(range(2), range(2), 'r-', lw=2)
			plt.plot(roc_datas['F'], roc_datas['H'], '.', zorder=1)
			#rw_max = [rw for rw in roc_datas if rw['H']-rw['F']==max(roc_datas['H']-roc_datas['F'])]
			#print rw_max
			max_val = max(roc_datas['H']-roc_datas['F'])
			for rw in roc_datas:
				if rw['H']-rw['F']==max_val:
					rw_max = rw
					break
			#
			roc_maxes += [rw_max]
			plt.plot(rw_max['F'], rw_max['H'], 'r*', ms=15, zorder=3)
			plt.plot(rw_max['F'], rw_max['H'], 'k*', ms=18, zorder=2)
			plt.xlabel('False alarm rate $F$')
			plt.ylabel('Hit rate $H$')
			plt.title('Section ID %d, %s metric' % (sec_id, metric_names.get(f_gt_lt, 'other')))
			plt.savefig('%s/roc_interval_metric_%d_%s.png'% (output_dir, sec_id, metric_names.get(f_gt_lt, 'other')))
			roc_datas.dump('%s/roc_interval_metric_%d_%s.npy'% (output_dir, sec_id, metric_names.get(f_gt_lt, 'other')))
			#
		mpp_loop_index+=1
	#
	#rw_maxes = numpy.core.records.fromarrays(zip(*roc_maxes), names='
	F,H = zip(*[[rw['F'], rw['H']] for rw in roc_maxes])
	plt.figure(fignum)
	plt.clf()
	plt.plot(range(2), range(2), 'r-', lw=2)
	plt.plot(F,H, 'o')
	plt.xlabel('False alarm rate $F$')
	plt.ylabel('Hit rate $H$')
#
def interval_forecast_mc_set(section_ids=emc_sections, nits=1000, dt_min=None, dt_max=None, nyquist_min=.2, nyquist_max=1.2, m0=7.0, f_gt_lt=operator.lt, percentile_lower=.1, percentile_upper=.8, output_dir='interval_metric_output', fignum=0, fignum_1=1):
	#
	# churn out a bunch of MC rate forecast analyses.
	#
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	metric_names={operator.lt:'LT', operator.gt:'GT', operator.le:'Leq', operator.ge:'Geq'}
	#
	roc_maxes = []
	#
	for sec_id in section_ids:
		roc_datas = interval_forecast_mc(nits=nits, dt_min=dt_min, dt_max=dt_max, nyquist_min=nyquist_min, nyquist_max=nyquist_max, CFF=None, m0=m0, do_plot=False, fnum=0, set_name=None, f_gt_lt=f_gt_lt, section_id=sec_id, detail=False, percentile_lower=percentile_lower, percentile_upper=percentile_upper)
		#
		# might not get any data back:
		if roc_datas==None or len(roc_datas)==0: continue
		#
		plt.figure(fignum)
		plt.clf()
		plt.plot(range(2), range(2), 'r-', lw=2)
		plt.plot(roc_datas['F'], roc_datas['H'], '.', zorder=1)
		#rw_max = [rw for rw in roc_datas if rw['H']-rw['F']==max(roc_datas['H']-roc_datas['F'])]
		#print rw_max
		max_val = max(roc_datas['H']-roc_datas['F'])
		for rw in roc_datas:
			if rw['H']-rw['F']==max_val:
				rw_max = rw
				break
		#
		roc_maxes += [rw_max]
		plt.plot(rw_max['F'], rw_max['H'], 'r*', ms=15, zorder=3)
		plt.plot(rw_max['F'], rw_max['H'], 'k*', ms=18, zorder=2)
		plt.xlabel('False alarm rate $F$')
		plt.ylabel('Hit rate $H$')
		plt.title('Section ID %d, %s metric' % (sec_id, metric_names.get(f_gt_lt, 'other')))
		plt.savefig('%s/roc_interval_metric_%d_%s.png'% (output_dir, sec_id, metric_names.get(f_gt_lt, 'other')))
		roc_datas.dump('%s/roc_interval_metric_%d_%s.npy'% (output_dir, sec_id, metric_names.get(f_gt_lt, 'other')))
		#
	#
	#rw_maxes = numpy.core.records.fromarrays(zip(*roc_maxes), names='
	F,H = zip(*[[rw['F'], rw['H']] for rw in roc_maxes])
	plt.figure(fignum)
	plt.clf()
	plt.plot(range(2), range(2), 'r-', lw=2)
	plt.plot(F,H, 'o')
	plt.xlabel('False alarm rate $F$')
	plt.ylabel('Hit rate $H$')
		
	
#
def interval_forecast_mc(nits=1000, dt_min=None, dt_max=None, nyquist_min=.2, nyquist_max=1.2, CFF=None, m0=7.0, do_plot=False, fnum=0, set_name=None, f_gt_lt=operator.lt, section_id=16, detail=False, percentile_lower=.1, percentile_upper=.8, dump_file=None):
	#
	# quick MC evaluation of a rate based metric, interval_forecast_1().
	#
	if (CFF==None or len(ary_in)==0) and section_id!=None:
		#ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id		
		CFF=combine_section_CFFs(section_id)		# ... and this function knows to use the string format above...
		#											# ... and note we do CFF=ary_in in the if-else clause below...
		#
	if isinstance(CFF, str):
		CFF = numpy.load(ary_in)
	if set_name==None: set_name='unnamed_CFF'
	#
	if dt_min==None or dt_max==None:
		delta_ts = [t-CFF['event_year'][j] for j,t in enumerate(CFF['event_year'][1:])]
		delta_ts.sort()
		if dt_min==None: dt_min=delta_ts[int(percentile_lower*len(delta_ts))]
		if dt_max==None: dt_max=delta_ts[int(percentile_upper*len(delta_ts))]
		#
		print "sampling intervals between %f, %f" % (dt_min, dt_max)
	#
	ROC_prams = []
	#
	# get mean m0 interval (this all needs to be reorganized and put into a class structure for efficiency...)
	m0s = [rw['event_year'] for rw in CFF if rw['event_magnitude']>=m0]
	#return m0s
	mean_dt_m0 = numpy.mean([x-m0s[j] for j,x in enumerate(m0s[1:])])
	#
	R_alpha = random.Random()
	R_dt = random.Random()
	delta_alpha = nyquist_max-nyquist_min
	delta_t0 = dt_max-dt_min
	#
	#
	for j in xrange(nits):
		dt_0 = dt_min + delta_t0*R_dt.random()
		alpha = nyquist_min + delta_alpha*R_alpha.random()
		#
		alert_segments = interval_forecast_1(ary_in=CFF, m0=m0, dt_0=dt_0, nyquist_factor=alpha, set_name=set_name, f_gt_lt=f_gt_lt, section_id=section_id, detail=False)
		#
		roc_prams = evaluate_alert_segments(alert_segments=alert_segments, CFF=CFF, section_id=section_id, m0=m0, do_plot=False)
		roc_prams.update({'dt_0':dt_0, 'nyquist_factor':alpha, 'm0':m0, 'dt_m0':mean_dt_m0})
		ROC_prams+=[roc_prams]
		#
	#
	# and let's make a rec array:
	cols = roc_prams.keys()
	output = [[rw[key] for key in cols] for rw in ROC_prams if not (rw==None or None in rw.values())]
	if not (output==None or len(output)==0):
		output = numpy.core.records.fromarrays(zip(*output), names=cols, formats = [type(x).__name__ for x in output[0]])
	if not dump_file==None: output.dump(dump_file)
	#
	return output
#
def interval_forecast_1(ary_in=None, m0=7.0, dt_0 = 0.0, nyquist_factor=.5, set_name=None, f_gt_lt=operator.lt, section_id=16, detail=False):
	'''
	# a rate (or interval) based forecast. same idea -- sample nyquist_factor, critical_rate for but gt/lt.
	#
	# note: lower interval (aka, dt<dt_0) produces information gain. dt>dt_0 yields information loss.
	#
	# output only the alert segments; then use evaluate_alert_segments() to get ROC data.
	# note that, at least in its native form, this outputs a binary alert: alert is on/off; there is no 'value' (aka, slope value b)
	# or degree associated. note that another model for this sort of work is to just output a full 1:1 alert set, aka
	# [ [time, alert_value]
	'''
	#
	# copied (mostly) from the first half of psa_metric_1():
	#
	# first, get the catalog data.
	if (ary_in==None or len(ary_in)==0) and section_id!=None:
		#ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id		
		ary_in=combine_section_CFFs(section_id)		# ... and this function knows to use the string format above...
		#											# ... and note we do CFF=ary_in in the if-else clause below...
		#
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		# otherwise, assume we've been passed a proper CFF object:
		CFF = ary_in
		# give the set a name  so we don't return the whole input data object...
		if set_name==None: set_name='unnamed_CFF'
	#
	if set_name==None:
		# still none? must be a string type ary_in...
		set_name=str(ary_in)
	#
	#print "get recurrence...** %d **" % len(CFF)
	recurrence_data = mean_recurrence(ary_in=CFF, m0=m0)
	nyquist_len = max(int(nyquist_factor*recurrence_data['mean_dN_fault']), 2)
	nyquist_time = nyquist_factor*recurrence_data['mean_dT']
	#
	#print "nyquist_len: ", nyquist_len
	#
	#print "recurrenced. get trend data..."
	#
	#trend_data = get_trend_analysis(ary_in=CFF, nyquist_len = nyquist_len, nyquist_time=nyquist_time)
	# note, mean intervals over nyquist_len events, including the very first event (because we can calc. mean intervals directly from time data).
	interval_data = [[j, rw[0], rw[1], (rw[1]-CFF['event_year'][j])/float(nyquist_len)] for j, rw in enumerate(zip(CFF['event_number'], CFF['event_year'])[nyquist_len-1:])]
	interval_data = numpy.core.records.fromarrays(zip(*interval_data), names=['index', 'event_number', 'event_year', 'interval'], formats=[type(x).__name__ for x in interval_data[0]])
	#
	# an index for fast reference...
	CFF_dict = {x['event_number']:x for x in CFF}
	
	print "trend lengths: ", len(interval_data), len(CFF), nyquist_len
	max_n = len(interval_data)
	#
	# trend_data includes columns: ['event_number', 'event_year', 'lin_fit_a', 'lin_fit_b', 'rb_ratio', 'interval_rate_ratios']
	#
	# first, just get the total time under alert:
	alert_segments = [[]]		# a collection of lists...
	#alert_segments = []
	#
	for i, rw in enumerate(interval_data):
		# when we find dt</>dt_0, we issue an alert until the next earthquake -- unless this one was 'big',
		# in this case m>=7.
		#
		if i>=(max_n-1): break
		#
		#
		this_dt = rw['interval']
		#
		#
		#if this_b >= b_0:
		if not f_gt_lt(this_dt, dt_0):		# abstracting the gt/lt option...
			# null case; does not meet "alert" criteria.
			# if we've been collecting "alert" events, stop. if not, just troll along...
			if len(alert_segments[-1])>0:
				#
				alert_segments+= [[]]
		
		else:
			# f_gt_lt(this_b, b_0) is True
			# accelerating case (if f_gt_lt = operator.lt):
			#this_mag = CFF[-max_n:][i+1]['event_magnitude']
			this_mag = CFF_dict[rw['event_number']]['event_magnitude']
			#
			alert_year = interval_data[i]['event_year']
			if len(alert_segments[-1])==0:
				# this is the first entry into this alert segment.
				#alert_year = interval_data[i]['event_year']
				if detail:
					alert_segments[-1]+=[[alert_year, interval_data[i]['interval']]]
				else:
					alert_segments[-1]+=[alert_year]
			#
			# these two cases are additional entries into this alert segment. if it's a small event, extend alert to the
			# next event (aka, in real live, we'd turn "on" the alert and reevaluate at the next event. if it's a big event,
			# it's the one we've been lookng for. we have a rule to turn off an alert immediately after a big event... though we might
			# reevaluate this since we've started to favor the operator.gt() version of this metric.
			if this_mag<m0:
				# add the *next* event as the alert (aka, issue an alert until we have new data).
				#print "len alert_seg: %d" % len(alert_segments[-1])
				#
				# generalize language a bit:
				#alert_year = interval_data[i+1]['event_year']
				if detail:
					alert_segments[-1]+=[[alert_year, this_dt]]
				else:
					alert_segments[-1]+=[alert_year]
				#pass
				#
			if this_mag>=m0:
				# this is "the" earthquake. add this entry (it's probably already there) from the previous entry.
				#
				#alert_year = trend_data[i]['event_year']
				if detail:
					alert_segments[-1]+=[[alert_year, this_dt]]
				else:
					alert_segments[-1]+=[alert_year]
			#
		#
	#
	while len(alert_segments)>0 and len(alert_segments[-1])==0: alert_segments.pop()
	#
	return alert_segments
#
def check_psa_metric(section_id=16, m0=7.0, fignum=0, nyquist_factor=.5):
	# (mostly) diagnostic function to validate psa_forecast_1(). also a couple of pretty complex figures.
	CFF = combine_section_CFFs(section_id)
	recurrence_data = mean_recurrence(ary_in=CFF, m0=m0)
	# get opt. data:
	opt_data = numpy.load('dumps/gji_roc_lt_500/roc_sec_lt_%d_nits_1000_allprams.npy' % section_id)
	opt_data.sort(key=lambda x: (x['H']-x['F']))
	best_fits = opt_data[-1]
	b_0 = best_fits['b']
	nyquist_factor = best_fits['nyquist_factor']
	#
	nyquist_len = max(int(nyquist_factor*recurrence_data['mean_dN_fault']), 2)
	#nyquist_len = best_fits['nyquist_factor']
	
	fc=psa_forecast_1(ary_in=CFF, f_gt_lt=operator.lt, b_0=b_0, nyquist_factor=nyquist_factor)
	trends = get_trend_analysis(ary_in=CFF, nyquist_len=nyquist_len)
	#
	# see that CFF, trends line up. basically, see that the interval values match for the two data sets. match the data-sets up by event_date
	# and/or event_id. note that, in the future, using a dict. with {'event_number':[stuff...]} might be a bit more fool-proof for various
	# analyses.
	j_cff = 0
	j_trends = 0
	#
	while CFF[j_cff]['event_year']< trends[j_trends]['event_year']: j_cff+=1
	while trends[j_trends]['event_year']<CFF[j_cff]['event_year']: j_trends+=1
	print "CFF, trends match: %d, %d :: %f" % (j_cff, j_trends, CFF[j_cff]['event_year']-CFF[j_cff-1]['event_year'])
	print CFF[j_cff]
	print trends[j_trends]
	#
	# see that the intervals are the same:
	print "checking intervals..."
	for k_cff,k_trends in zip(xrange(j_cff, len(CFF)), xrange(j_trends, len(trends))):
		interval_cff = CFF[k_cff]['event_year']-CFF[k_cff-1]['event_year']
		if interval_cff!=trends['intervals'][k_trends]:
			print "intervals: %f, %f *** %s ***" % (interval_cff, trends['intervals'][k_trends], interval_cff==trends['intervals'][k_trends])
		#
	print "finished interval agreement check. now, check for fit agreement."
	
	for k_cff,k_trends in zip(xrange(j_cff, len(CFF)), xrange(j_trends, len(trends))):
		# fitting algorithm is like: b,a = numpy.linalg.lstsq([X,w], Y)
		b_trends = trends['lin_fit_b'][k_trends]
		a_trends = trends['lin_fit_a'][k_trends]
		#
		Y = [x for x in CFF['event_year'][k_cff+1-nyquist_len:k_cff+1]-CFF['event_year'][k_cff-nyquist_len:k_cff]]	# so. like T[1:10] - T[0:9]
		#																											 for alpha=10
		X = [[x,1.0] for x in CFF['event_year'][k_cff+1-nyquist_len:k_cff+1]]
		b,a = numpy.linalg.lstsq(X,Y)[0]
		#
		if b!=b_trends or a!=a_trends: print "b: %f/%f, a: %f/%f" % (b_trends, a_trends, b, a)
		#if k_cff>j_cff + 10: break
	print 'fit agreement finished. you should see any case where b!=b_trends (aka, local fit != (pre)scripted fit).'
	
	print '****'
	#
	fc_lens = [len(x) for x in fc]
	fc_lens.sort()
	#fc_delta_ts_2 = [max(x)-min(x) for x in fc]
	fc_delta_ts = [x[-1]-x[0] for x in fc]
	#if fc_delta_ts!=fc_delta_ts_2:
	#	print 'problem with delta ts'
	#	return fc_delta_ts, fc_delta_ts_2
	fc_delta_ts.sort()
	
	#
	print 'finished. now do some plotting...'
	
	#
	m0s = zip(*[[rw['event_year'], rw['event_magnitude']] for rw in CFF if rw['event_magnitude']>m0])
	y_b_max=1.5
	y_b_min = -1.5
	min_y_int=min(trends['intervals'])
	max_y_int=max(trends['intervals'])
	#
	f=plt.figure(fignum)
	plt.clf()
	plt.ion()
	#ax_b = plt.gca()
	ax_rate=f.add_axes([.1, .1, .8, .4])
	ax_b=f.add_axes([.1, .55, .8, .4], sharex=ax_rate)
	#plt.figure(fignum+1)
	#plt.clf()
	#ax_rate = plt.gca()
	ax_rate.set_yscale('log')
	ax_b.set_ylabel('Interval slope $b=d(\\Delta t)/dt$')
	ax_rate.set_ylabel('Interval $\\Delta t_i = t_i - t_{i-1}$')
	ax_rate.set_xlabel('Simulation time $t$ years')
	#
	ax_b.plot([CFF['event_year'][0], CFF['event_year'][-1]], [0., 0.], 'k-', zorder=0)
	ax_b.plot([CFF['event_year'][0], CFF['event_year'][-1]], [b_0, b_0], 'k--', zorder=0)
	#
	ax_b.plot(trends['event_year'], trends['lin_fit_b'], '-', label='Interval slope $b$')
	ax_b.vlines(m0s[0], y_b_min, y_b_max, color='r', lw=2, linestyle='solid')
	#
	#ax_rate.plot(trends['event_year'], [1./x for x in trends['intervals']], '-', label='rates')
	ax_rate.plot(trends['event_year'], trends['intervals'], '-', label='intervals')
	ax_rate.vlines(m0s[0], min_y_int, max_y_int, color='r', lw=2, linestyle='solid')
	#
	#lbl_str = 'Alert interval'
	for j,rw in enumerate(fc):
		#if j>0 and 1==1:
		#	lbl_str = None
		ax_b.fill_between([min(rw), max(rw)], [y_b_min, y_b_min], [y_b_max, y_b_max], color='m', alpha=.2, zorder=6)
		ax_rate.fill_between([min(rw), max(rw)], [min_y_int, min_y_int], [max_y_int, max_y_int], color='m', alpha=.2, zorder=6)
	#b_artist = plt.Line2D((0,0),(0,0), color='m', linestyle='-')
	b_artist = mpl.patches.Patch(facecolor='m', edgecolor='m', alpha=.2)
	#
	ax_rate.legend(loc=0, numpoints=1)
	#
	#ax_b.legend(loc=0, numpoints=1)
	b_handles, b_labels = ax_b.get_legend_handles_labels()
	ax_b.legend([h for h in b_handles]+[b_artist], [l for l in b_labels] + ['Alert intervals'], loc=0, numpoints=1)
	#
	plt.figure(fignum+2)
	plt.clf()
	ax_n = plt.gca()
	#ax_n.set_yscale('log')
	#ax_n.set_xscale('log')
	ax_n.plot(fc_lens, [x/float(len(fc_lens)) for x in xrange(1, len(fc_lens)+1)], 'b.-', label='alert sequence length')
	ax_n.legend(loc='best', numpoints=1)
	ax_n.set_xlabel('alert segment lengths $n$')
	#
	ax_dt = plt.twiny(plt.gca())
	#ax_dt.set_yscale('log')
	#ax_dt.set_xscale('log')
	ax_dt.plot(fc_delta_ts, [x/float(len(fc_delta_ts)) for x in xrange(1, len(fc_delta_ts)+1)], 'g.-', label='alert intervals $\\Delta t$')
	ax_dt.legend(loc='best', numpoints=1)
	ax_dt.set_xlabel('Alert durations $\\Delta t$')
	#
	#plt.ylabel('cumulative count $N$')
	ax_n.set_ylabel('Cumulative probability $P(N)$ that $len(alert)<n$')
	ax_dt.set_ylabel('Cumulative probability $P(N)$ that $\\Delta t(alert)<\\Delta t$')
	#
#
def plot_psa_metric_figure(CFF=None, section_id=None, m0=7.0, fignum=0, nyquist_factor=None, b_0=None, opt_data='dumps/gji_roc_lt_500/roc_sec_lt_%d_nits_1000_allprams.npy', lw=2., **kwargs):
	# create forecast metric figure(s). we'll work on generalizing this to use multiple forecasts later. for now, let's be specific
	# to psa_forecast_1()
	#
	# handle inputs:
	if (CFF==None or len(CFF)==0):
		#ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id		
		CFF=combine_section_CFFs(section_id)		# ... and this function knows to use the string format above...
		#											# ... and note we do CFF=ary_in in the if-else clause below...
	if isinstance(CFF, str):
		CFF = numpy.load(CFF)
	#
	# fetch some possible kwargs...
	fs_label  = kwargs.get('fs_label', 12)
	fs_title  = kwargs.get('fs_title', 12)
	fs_legend = kwargs.get('fs_legend', 12)
	#
	# get opt. data:
	# optimization data will probably be provided as a string pointing to an optimization run output, which will be list of dicts[ {},{}...]
	# from a monte-carlo algorithm. for now, what is the format?
	if isinstance(opt_data, str):
		# were we specific, or did we just pass a string format (this part needs to eventually be quite a bit smarter).
		if '%d' in opt_data:
			opt_data = opt_data % section_id		# ... and if section_id is null, there's not much to be done... eventually, we need 
													# to re-organize the derived data so that CFF, etc. contain more meta-data.
		opt_data = numpy.load(opt_data)
			#
	# ... and at this point, opt_data will be an array-o-dicts (or at least it should be)
	# now, determine forecast parameterization from optimized bits.
	opt_data.sort(key=lambda x: (x['H']-x['F']))
	best_fits = opt_data[-1]
	b_0 = best_fits['b']
	nyquist_factor = best_fits['nyquist_factor']
	#
	recurrence_data = mean_recurrence(ary_in=CFF, m0=m0)
	nyquist_len = max(int(nyquist_factor*recurrence_data['mean_dN_fault']), 2)
	#
	# alert segments:
	fc=psa_forecast_1(ary_in=CFF, f_gt_lt=operator.lt, b_0=b_0, nyquist_factor=nyquist_factor)
	trends = get_trend_analysis(ary_in=CFF, nyquist_len=nyquist_len)
	#
	fc_lens = [len(x) for x in fc]
	fc_lens.sort()
	fc_delta_ts = [x[-1]-x[0] for x in fc]
	fc_delta_ts.sort()
	#
	#
	m0s = zip(*[[rw['event_year'], rw['event_magnitude']] for rw in CFF if rw['event_magnitude']>m0])
	y_b_max=1.5
	y_b_min = -1.5
	min_y_int=min(trends['intervals'])
	max_y_int=max(trends['intervals'])
	#
	f=plt.figure(fignum)
	plt.clf()
	plt.ion()
	#ax_b = plt.gca()
	ax_rate=f.add_axes([.1, .1, .8, .4])
	ax_b=f.add_axes([.1, .55, .8, .4], sharex=ax_rate)
	ax_mags = ax_rate.twinx()
	#plt.figure(fignum+1)
	#plt.clf()
	#ax_rate = plt.gca()
	ax_rate.set_yscale('log')
	ax_b.set_ylabel('Interval slope $b=d(\\Delta t)/dt$', size=fs_label)
	ax_rate.set_ylabel('Interval $\\Delta t_i = t_i - t_{i-1}$', size=fs_label)
	ax_rate.set_xlabel('Simulation time $t$ years', size=fs_label)
	#
	ax_b.plot([CFF['event_year'][0], CFF['event_year'][-1]], [0., 0.], 'k-', zorder=0, lw=lw)
	ax_b.plot([CFF['event_year'][0], CFF['event_year'][-1]], [b_0, b_0], 'k--', zorder=0, lw=lw)
	#
	ax_b.plot(trends['event_year'], trends['lin_fit_b'], '.-', label='Interval slope $b$', lw=lw)
	ax_b.vlines(m0s[0], y_b_min, y_b_max, color='r', linestyle='solid', zorder=4, lw=lw)
	#
	#ax_rate.plot(trends['event_year'], [1./x for x in trends['intervals']], '-', label='rates')
	ax_rate.plot(trends['event_year'], trends['intervals'], '.-', label='intervals', lw=lw)
	#ax_rate.vlines(m0s[0], min_y_int, max_y_int, color='r', linestyle='solid', lw=lw)
	#
	min_mag = min(CFF['event_magnitude'])-.25
	ax_mags.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='g', alpha=.7, zorder=3,lw=lw)
	ax_mags.vlines(m0s[0], min_mag, m0s[1], color='r', lw=lw, linestyle='solid', zorder=4)
	#
	#lbl_str = 'Alert interval'
	for j,rw in enumerate(fc):
		#if j>0 and 1==1:
		#	lbl_str = None
		ax_b.fill_between([min(rw), max(rw)], [y_b_min, y_b_min], [y_b_max, y_b_max], color='m', alpha=.2, zorder=6)
		ax_rate.fill_between([min(rw), max(rw)], [min_y_int, min_y_int], [max_y_int, max_y_int], color='m', alpha=.2, zorder=6)
	#
	# legend?
	#b_artist = mpl.patches.Patch(facecolor='m', edgecolor='m', alpha=.2)
	#ax_rate.legend(loc=0, numpoints=1)
	#
	#b_handles, b_labels = ax_b.get_legend_handles_labels()
	#ax_b.legend([h for h in b_handles]+[b_artist], [l for l in b_labels] + ['Alert intervals'], loc=0, numpoints=1)
	#
	plt.figure(fignum+2)
	plt.clf()
	ax_n = plt.gca()
	#ax_n.set_yscale('log')
	#ax_n.set_xscale('log')
	ax_n.plot(fc_lens, [x/float(len(fc_lens)) for x in xrange(1, len(fc_lens)+1)], 'b.-', label='alert sequence length')
	ax_n.legend(loc='best', numpoints=1, prop={'size':fs_legend})
	ax_n.set_xlabel('alert segment lengths $n$', size=fs_label)
	#
	ax_dt = plt.twiny(plt.gca())
	#ax_dt.set_yscale('log')
	#ax_dt.set_xscale('log')
	ax_dt.plot(fc_delta_ts, [x/float(len(fc_delta_ts)) for x in xrange(1, len(fc_delta_ts)+1)], 'g.-', label='alert intervals $\\Delta t$')
	ax_dt.legend(loc='best', numpoints=1, prop={'size':fs_legend})
	ax_dt.set_xlabel('Alert durations $\\Delta t$', size=fs_label)
	#
	#plt.ylabel('cumulative count $N$')
	ax_n.set_ylabel('Cumulative probability $P(N)$ that $len(alert)<n$', size=fs_label)
	ax_dt.set_ylabel('Cumulative probability $P(N)$ that $\\Delta t(alert)<\\Delta t$', size=fs_label)
	#
	# note: evaluate_alert_segments() can also be used to create a figure.
	#z=evaluate_alert_segments(alert_segments=fc, CFF=CFF, section_id=section_id, m0=7.0, do_plot=True, fnum=11)
	
#
def forecast_metric_1(ary_in=None, m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False, fnum=0, set_name=None, f_gt_lt=operator.gt, section_id=16, over_ride=False):
	'''
	# depricated:
	# this can be tossed in exchange for the more generalized format:
	# x = psa_forecast_1(CFF)		# (and other forecast metrics can be substituted of course)
	# y = evaluate_alert_segments(x,CFF)
	# so, until we get rid of this not the over_ride=False default parameter. if True, then run the old code. otherwise, wrap the new bits.
	# 
	#
	# forecast based on seismic acceleration. specifically, for a set of data (a fault-section time-series),
	# find the slopes of the inter-occurrence intervals. slope<0 implies acceleration, implies hazard.
	# the evaluation metric will be basically a ROC type:
	# 1) total time "at alert" serves as "false alarm"
	# 2) fraction of "predicted" events. this will vary. for starters:
	#     - slope at time t_i of event or time t_i-1 of the previous event is below b_0 (typically, b_0=0.0).
	#     - the "_i-1" condition is fair because 1) we don't know the slope is changin until the even happens
	#       and 2) in a real scenario, we'd have aftershocks, which woluld turn the post-seismic slope positive.
	#
	# parameters:
	# m0: recurrence magnidue
	# b_0: critical slope; issue alert for f_gt_lt(b,b_0)==True, aka, b>b_0 (see below)
	# nyquist_factor: fraction of inter-event sequence length to sample
	# do_plot: do plots ?
	# fnum: plot on figure(fnum)...
	# set_name: a name to give the set if ary_in is an array, not a string, or if we want to specify it.
	# (we should use the plotting routine from vc_parser.plot_CFF_ary(), namely the top subfig. )
	# f_gt_lt: function for greater/lesser evaluation... use operator.{gt, ge, lt, le}, or any two parameter input.
	'''
	#
	#
	if (ary_in==None or len(ary_in)==0) and section_id!=None:
		ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		# otherwise, assume we've been passed a proper CFF object:
		CFF = ary_in
		# give the set a name  so we don't return the whole input data object...
		if set_name==None: set_name='unnamed_CFF'
	#
	if set_name==None:
		# still none? must be a string type ary_in...
		set_name=ary_in
	#
	########################
	# depricating....
	if not over_ride:
		# ary_in=None, m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False, fnum=0, set_name=None, f_gt_lt=operator.lt, section_id=16, detail=False
		x=psa_forecast(ary_in=ary_in, m0=m0, b_0=b_0, nyquist_factor=nyquist_factor, do_plot=do_plot, fnum=fnum, set_name=set_name, f_gt_lt=f_gt_lt, section_id=section_id, detail=False)	# though we might want detail=True...
		# evaluate_alert_segments(alert_segments=None, CFF=None, section_id=None, m0=7.0, do_plot=True, fnum=0)
		return evaluate_alert_segments(alert_segments=x,CFF=CFF, section_id=section_id, do_plot=do_plot, fnum=fnum)

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
		#if this_b >= b_0:
		if not f_gt_lt(this_b, b_0):		# abstracting the gt/lt option...
			# null case
			# if we've been collecting "alert" events, stop. if not, just troll along...
			if len(alert_segments[-1])>0:
				#
				alert_segments+= [[]]
		
		else:
			# f_gt_lt(this_b, b_0) is True
			# accelerating case (if f_gt_lt = operator.lt):
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
				#pass
				#
			if this_mag>=m0:
				# this is "the" earthquake. add this entry (it's probably already there) from the previous entry.
				#
				alert_segments[-1]+=[[trend_data[i]['event_year'], this_b]]
			#
		#
	#
	while len(alert_segments)>0 and len(alert_segments[-1])==0: alert_segments.pop()
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
	n_total = 0
	#
	# was an alert issued at the time of an m>m0 event?
	# we should clean this up, but for now, make and use dictionaries...
	#alert_dict = {x[0]:x[1] for x in 
	#
	# this approach is problematic and i think might be producing some artifact effects.
	#alert_dict = {}
	#for seg in alert_segments:
	#	#
	#	for rw in seg:
	#		alert_dict[rw[0]] = rw[1]
	#
	j_alert_start=0
	for i, rw in enumerate(CFF):
		if rw['event_magnitude']<m0: continue
		n_total+=1
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
		#
		## note: this is keying off a float index, so it has potential to be unreliable.
		#if alert_dict.has_key(rw['event_year']):						
		#	n_predicted += 1
		#else:
		#	n_missed += 1
		#
		try:
			while len(alert_segments)>j_alert_start and rw['event_year']>alert_segments[j_alert_start][-1][0] and j_alert_start<(len(alert_segments)-1): j_alert_start+=1
		except:
			print "alert indexing error. ", len(alert_segments)
		#
		for alert in alert_segments[j_alert_start:]:
			if rw['event_year'] > alert[0][0] and rw['event_year'] <= alert[-1][0]:
				n_predicted +=1
				#n_missed -=1
				break

	n_missed = n_total-n_predicted
	#
	if do_plot:
		# diagnostic plots of forecast metric:
		plt.figure(fnum)
		plt.clf()
		plt.title('VC PSA Hazard Metric')
		#
		ax_mags = plt.gca()
		ax_metric = ax_mags.twinx()		
		#
		ax_ints = ax_metric.twinx()
		ax_ints.set_yscale('log')
		ax_metric.set_yscale('linear')
		#
		min_mag = min(CFF['event_magnitude'])
		max_mag = max(CFF['event_magnitude'])
		intervals = [x-CFF['event_year'][i] for i,x in enumerate(CFF['event_year'][1:])]
		#
		#metric_pad_factor = min_mag
		#metric_pad_factor = 0.
		#min_metric = 0.
		min_metric = alert_segments[0][0][1]
		max_metric = alert_segments[0][0][1]
		#
		# do a quick spin to get min/max values and other useful stats:
		for segment in alert_segments:
			X,Y = zip(*segment)
			min_metric = min(min(Y), min_metric)
			max_matric = max(max(Y), max_metric)
		#
		for segment in alert_segments:
			X,Y = zip(*segment)
			min_metric = min(min(Y), min_metric)
			#
			# ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b']  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='m', zorder=1, alpha=.5)
			#ax_trend2.fill_between([x['event_year'] for x in trend_data], [1.  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0.0 for x in trend_data], color='m', zorder=1, alpha=.25)
			#
			# show the metric value:
			#ln_metric = ax_metric.fill_between(X,[y for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [y<0. for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			ln_metric = ax_metric.fill_between(X,[y for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [f_gt_lt(y,0.) for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			ax_metric.plot([x for i,x in enumerate(X) if Y[i]<0], [0.0 for y in Y if y<0], 'm-')
			ax_metric.plot([x for i,x in enumerate(X) if Y[i]<0], [y for y in Y if y<0], 'm-')
			#
			# and just an "alert!" box:
			#ln_metric = ax_metric.fill_between(X,[min_metric for y in Y],y2=[0.0 for y in Y], color='m', alpha=.15, where = [y<0. for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			ln_metric = ax_metric.fill_between(X,[min_metric for y in Y],y2=[0.0 for y in Y], color='m', alpha=.15, where = [f_gt_lt(y,0.) for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			#
			#ax_mags.fill_between(X, [min_mag for x in X], [m0 for x in X], zorder=5, alpha=.2, color='c')
		#
		ln_ints = ax_ints.plot(CFF['event_year'][1:], intervals, 'b.-', lw=2, ms=7, label='Intervals $\Delta t = t_i - t_{i-1}$')
		ax_metric.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		ln_mags = ax_mags.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='g', alpha=.7, lw=1.5, label='magnitudes')
		X_big_mags, Y_big_mags = zip(*[[x['event_year'], x['event_magnitude']] for x in CFF if x['event_magnitude']>m0])
		ax_mags.vlines(X_big_mags, [min_mag for x in Y_big_mags], Y_big_mags, color='r', alpha=.9, lw=2.75)
		#
		# cosmetics:
		ax_mags.set_xlabel('Event Year $t$')
		ax_ints.set_ylabel('Inter-event interval $\Delta t$')
		ax_metric.yaxis.set_ticks([])
		ax_mags.set_ylabel('Earthquake magnitude $m$')
		#
		ax_metric.set_ylim(ymin=1.15*min_metric, ymax = -min_metric/10.)
	#
	print "preliminary report:"
	#H = float(n_predicted)/(float(n_predicted)+n_missed)
	H=float(n_predicted)/float(n_total)
	F = total_alert_time/total_total_time
	#print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, total_alert_time/total_total_time)
	#print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, float(n_predicted)/(float(n_predicted)+n_missed))
	#print "score: H-F: %f " % (float(n_predicted)/(float(n_predicted)+n_missed) - total_alert_time/total_total_time)
	print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, F)
	print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, H)
	print "score: H-F, H/F: %f / %f " % (H-F, H/F)

	
	#
	#return alert_segments
	# , 'alert_segments':alert_segments
	#return {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'alert_segments':alert_segments, 'ary_in_name':ary_in, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
	return {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'alert_segments':alert_segments, 'ary_in_name':set_name, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
#
def simple_mpp_optimizer(sections=[], section_names=None, start_year=0., m0=7.0, b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=1000, keep_set=True, dump_file=None, n_cpus=None, f_gt_lt=operator.lt):
	# (production run function for the ROC diagram and fits).
	#
	# mpp implementation of simple_metric_optimizer() to optimize a set of forecasts. note that sections[] can include not only 
	# individual fault segments but groups of segments like: sections = [25, 23, 30, 31, 32, 33, [30,31,32,33], ...] in the event that
	# we want to do faultwise analysis. also, we'll use this trick to do the composite set, rather than use a "do_composite" flag or something.
	#
	# trap some default shortcuts:
	if sections in ('emc', 'EMC'):
		sections = emc_sections + [emc_sections]
		section_names = ['%d' % x for x in emc_sections] + ['emc']
	if sections == 'napa':
		sections = napa_sections + [napa_sections]
		section_names = ['%d' % x for x in napa_sections] + ['napa']
	#
	# we can use a Pool():
	if n_cpus==None: n_cpus = mpp.cpu_count()
	P = mpp.Pool()
	pool_handlers = []
	pool_results = []
	#
	if section_names == None:
		section_names = []
		for sec_id in sections:
			if hasattr(sec_id, '__len__'):
				section_names += ['data_set(%d)' % len(sec_id)]
			else:
				section_names += ['section_%d' % sec_id]
			#
		#
	#
	print "simple_mpp: assign jobs to pool"
	for i, sec_id in enumerate(sections):
		#
		pool_handlers += [P.apply_async(simple_metric_optimizer, kwds={'CFF':{'sections':sec_id, 'start_year':start_year}, 'm0':m0, 'b_min':b_min, 'b_max':b_max, 'd_b':d_b, 'nyquist_min':nyquist_min, 'nyquist_max':nyquist_max, 'd_nyquist':d_nyquist, 'nits':nits, 'keep_set':keep_set, 'set_name':section_names[i], 'dump_file':None, 'f_gt_lt':f_gt_lt})]
	P.close()
	#
	print "simple mpp: fetch results ([R.get() for R in pool_handlers] )"
	# this (list comprehension approach) fails silently and awkwardly if no data are returned.
	#pool_results = [R.get()[0] for R in pool_handlers]
	#	
	for R in pool_handlers:
		R_return = R.get()
		if hasattr(R_return, '__len__'):
			pool_results+=[R_return[0]]
		else:
			pool_results+=[None]
		#pool_results+=[R.get()[0]]
	
	print "results fetched."
	#
	# make a recarray:
	#
	# duh... this won't work because we have a list data type (there may be a way to include lists, etc in recarrays, but i don't know
	# it.) do we really want to keep the alert segments? they're easy enough to recover given the paramset. let's just kill them.
	#
	# screw this. there are strings... or one string anyway, so we have to fix its length and all of that. do this alter (or not at all)
	# ... or maybe we'll get the recarray; for string, the type is 's%d' % s_len or 'S%d' % s_len (i think).
	#print pool_results[0]
	#pool_results = numpy.rec.array([rw.values() for rw in pool_results], names=pool_results[0].keys(), formats = [type(x).__name__ for x in pool_results[0].itervalues()])
	#
	#
	if dump_file!=None:
		with open(dump_file,'w') as f:
			cPickle.dump(pool_results, f)
	#
	return pool_results
#
def simple_metric_optimizer(CFF=None, m0=7.0, b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01,  nits=1000, keep_set=True, set_name='data_set', dump_file=None, f_gt_lt=operator.lt, f_score=operator.sub, section_id=16, opt_func=psa_forecast_1):
	'''
	# this does an MC thing to optimize a forecast metric for a CFF (presumably a section catalog like 'data/VC_CFF_timeseries_section_16.npy'
	# or the data therein. it can be run independently or (it's primary purpose) wrapped by simple_mpp_optimizer()
	#
	#, opt_func=forecast_metric_1, opt_args=[], opt_kwargs={}):
	
	# a simple MC optimizer:
	# CFF: it's interval data. a standard CFF file (like the default) will suffice, otherwise
	# a rec_array basically with a ['event_year'], ['event_magnitude'] columns i think. if the arg is a string, load the file,
	# otherwise assume it's a ready-to-rock array.
	#
	# nits: mc nits
	# f_gt_lt_operator: greater/lesser operator. if b<b_0 --> if operator.lt(b,b_0) --> (in general) --> if f_gt_lt(b, b_0)
	# f_score: a two-parameter score function (which can be modified using a lambda function as necessary of course.
	#          variables will be passed to f_score like score = f_score(H,F).
	#          excellent choices inlcude: opeartor.sub(), opeartor.div()
	#
	#### save these for later...
	# opt_func: optimization function
	# opt_arts: args for opt_func
	# opt_kwargs: kwargs for opt_func
	# opt_funct will be called like f(CFF, *args, **kwargs)
	########
	#
	# later, we need to generalize the evaluation part -- namely, how we determine what is optimal. on the ROC,
	# Y=float(fc_data['n_predicted'])/(fc_data['n_predicted'] + float(fc_data['n_missed']))
	# X=fc_data['total_alert_time']/fc_data['total_time']
	# and then we use some sort of steepness metric...
	'''
	if f_score == None:
		#f_score = operator.div
		f_score = operator.sub
	#
	#opt_func = forecast_metric_1
	# forecast_metric_1(ary_in='data/VC_CFF_timeseries_section_16.npy', m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False, fnum=0)
	#
	if section_id!=None and CFF==None:
		CFF = combine_section_CFFs(section_id)
	#
	if isinstance(CFF, dict):
		# for now, assume have the proper call signature in the key-val pairs (this will contain section/catalog names and stuff like that).
		CFF = combine_section_CFFs(**CFF)
	#
	if isinstance(CFF, str):
		CFF = numpy.load(CFF)
	#
	R_b   = random.Random()
	R_nyq = random.Random()
	delta_b = b_max-b_min
	delta_nyq = nyquist_max - nyquist_min
	# and let's do this right and randomly sample...
	if nits==None: nits = 1 + int(abs((delta_b/d_b)*((delta_nyq)/d_nyquist)))	# safe-guard for n->0
	#
	all_prams=[]
	best_prams={}
	best_score=0.	# assume it will be positive (we can control this).
	#
	## we can plot this to watch the progression:
	#plt.figure(11)
	#plt.clf()
	#plt.ion()
	#plt.plot(range(2), range(2), '-')
	#bestXY = []
	for i in xrange(nits):
		print "optimizer n_iteration: %d" % i
		this_b   = b_min       + delta_b*R_b.random()
		this_nyq = nyquist_min + delta_nyq*R_nyq.random()
		#
		# this needs to be abstracted to take a set of prams...
		# so eventually we'd want something like opt_func(**opt_prams)
		#fit_data = opt_func(ary_in=CFF, m0=m0, b_0=this_b, nyquist_factor=this_nyq, do_plot=False, set_name=set_name, f_gt_lt=f_gt_lt)
		alert_segments = opt_func(ary_in=CFF, m0=m0, b_0=this_b, nyquist_factor=this_nyq, do_plot=False, set_name=set_name, f_gt_lt=f_gt_lt)
		fit_data = evaluate_alert_segments(alert_segments=alert_segments, CFF=CFF, do_plot=False)
		fit_data['b']=this_b
		fit_data['nyquist_factor']=this_nyq
		#
		# get hits/falsies:
		if fit_data.has_key('H'):
			hit_rate = fit_data['H']
		else:
			hit_rate = float(fit_data['n_predicted'])/(float(fit_data['n_predicted'])+float(fit_data['n_missed']))
		#
		if fit_data.has_key('F'):
			false_alarm_rate = fit_data['F']
		else:
			false_alarm_rate = float(fit_data['total_alert_time'])/float(fit_data['total_time'])		# 
		#
		# and there are cases where we get zero alert time, which results in an error using geometric score algorithms
		# (aka, s=H/F). so we'll handle this as best we can. what about a None value? or, just move on? None might throw an error
		# when we cast into recarrays. numpy.nan? we'll see...
		#
		try:
			#this_score = hit_rate/alert_rate	# could also be hit_rate-alert_rate (see below).
			this_score = f_score(hit_rate, false_alarm_rate)
		except:
			print "***************\nscore exception\n****************"
			#this_score = None
			continue
		#
		# note: the more contemporary area-skill = ingegral(H(f)-f), aka the area between H(F) and H=F,
		# can be emulated by:
		# this _score = hit_rate-alert_rate		#, (which i think was our first metric) so we "integrate" under each point
		#										# these should produce similar, but maybe not identical results.
		#
		fit_data['score']=this_score
		print this_score, hit_rate, false_alarm_rate
		if this_score>best_score:
			best_score=this_score
			#best_prams = fit_data
			# for best_prams, let's forego the time-series of alert segments:
			#best_prams = {key:val for key,val in fit_data.iteritems() if key!='alert_segments'}
			best_prams = {key:val for key,val in fit_data.iteritems() if not hasattr(val, '__len__')}
			#
			# some diagnostic bits:
			#print "best prams(%f): %s" % (best_score, str([[key,val] for key, val in fit_data.iteritems() if not key in ('alert_segments', 'ary_in_name')]))
			#print "hit=%f, alert=%f, score=%f\n*****\n" % (hit_rate, alert_rate, this_score)
			#plt.plot([alert_rate], [hit_rate], 'o')
			#bestXY += [[alert_rate, hit_rate]]
			#plt.draw()
			
		#
		if keep_set:
			#all_prams+=[fit_data]
			all_prams += [{key:val for key,val in fit_data.iteritems() if key!='alert_segments'}]
		#
	#
	#plt.plot(zip(*bestXY)[0],zip(*bestXY)[1], '-')
	#
	if dump_file!=None:
		# we'll want (optionally) to save all the data, not just the best_prams. cPickle seems to be a bit twitchy
		# (seems to not recognize how to encode an embeded list-object and runs away, probably writing all the
		# memory in the world to the output file -- and so, the line below .dump([best, all]) probably
		# effectively fails (??). SO, if keep_set==True: dump up to three lists separately:
		# best_prams, all_prams{minus alerts list}, alerts list
		with open(dump_file, 'w') as f:
			#cPickle.dump([best_prams, all_prams], f)
			cPickle.dump(best_prams, f)
		if keep_set:
			all_prams_dump_file = dump_file.replace('.npy', '_allprams.npy')	# for a cheap fix...
			with open(all_prams_dump_file, 'w') as f:
				# alerts key will be: 'alert_segments'
				#cPickle.dump(
				cPickle.dump([{key:val for key,val in rw.iteritems() if key!='alert_segments'} for rw in all_prams], f)
			#
			# ... and for now, let's forego the alert_segments lists. it's probably just as fast to reproduce them...unless we want to
			# save archives for some reason.
	#
	print "[from simple_metric_optimizer()]: finished optimizing (%s)" % (best_prams)
	return best_prams, all_prams

	

def plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=0.0, nyquist_factor=.5, do_spp=False, do_plot=False, do_clf=True, n_cpus=None):
	'''
	# depricated? this script seems to be showing its age as well...
	#
	# simple ROC diagram. see the optimized one...
	# also note that the parameter list for forecast_metric_1() has changed, so generally speaking this 
	# function needs some maintenance or to be retired.
	#
	# scatter plot of hit_ratio vs alert_time_ratio for as many data as we throw at it.
	# note, this uses a single value(s) for (b_0, nyquist_factor). see optimized versions as well.
	#
	# ... and this works ok, but it is a bit limited in exactly what sets we can optimise. see simple_mpp_optimizer()
	# and simple_metric_optimzer() above.
	'''		
	# 
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
	#if resultses[0].has_key('alert_segments'): [x.pop('alert_segments') for x in resultses]
	resultses[0] = {key:val for key,val in resultses[0].iteritems() if key!='alert_segments'}	# a little bit faster and more elegant...
	return resultses
#
def plot_section_ROC_curve(roc_data=None, section_id=None, fignum=0, num_points=100, do_clf=True, label_str=None, markers='.-'):
	# construct an ROC curve for a section (or more generally, a collection of ROC optimizer data collected by some unknown parsing).
	# roc_data input is the 'full_pram' set from simple_metric_optimizer()[1] or the output file (which can be a string-filename
	# or a roc_data = numpy.load(roc_data); see for example numpy.load('dumps/roc_detail_data/roc_sec_16_allprams.npy')
	#
	# ... and later on, code up something for section_id so we can build a roc_data set on demand, if necessary.
	#
	if isinstance(roc_data, str): roc_data = numpy.load(roc_data)
	if isinstance(roc_data[0],dict):
		# list of dicts. let's turn it into a rec-array.
		col_names = [key for key,val in roc_data[0].iteritems() if not isinstance(val,str)]		#roc_data[0].keys()
		lst_data = [[rw[key] for key in col_names] for rw in roc_data]
		for i,rw in enumerate(roc_data):
			# H,F are pre-calculated now...
			#H = float(rw['n_predicted'])/(rw['n_predicted'] + rw['n_missed'])
			#if rw['total_time']!=0.:
			#	F = float(rw['total_alert_time'])/rw['total_time']
			#else:
			#	F=numpy.nan
			H=rw['H']
			F=rw['F']
			#
			score_lin = H-F
			score_geom = H/F
			#
			#lst_data[i] += [H,F, score_lin, score_geom]
			lst_data[i] += [score_lin, score_geom]
			#for j in xrange(lst_data[i]):
			#	if lst_data[i][j]==None: lst_data[i][j]=numpy.nan
			lst_data[i] = [x if x!=None else numpy.nan for x in lst_data[i]]
			#
		col_formats = [type(roc_data[0][key]).__name__ for key in col_names]
		#col_names += ['H', 'F', 'score_lin', 'score_geom']
		col_names += ['score_lin', 'score_geom']
		while len(col_formats)<len(col_names): col_formats+=[type(1.0).__name__]
		#
		roc_data = numpy.core.records.fromarrays(zip(*lst_data), names=col_names, formats = col_formats)
		#
	# now we have a rec-array object like [[total_time, b, alpha, n_missed, alert_time, n_predicted, score, H, F, score_lin, score_geom] ]
	# we want to bin these by F and plot the max value for each bin (or the distribution or whatever). H,F are always 0<{F,H}<1
	#
	bin_size=1.0/num_points
	#
	roc_curve_data = [[] for x in xrange(num_points)]
	#roc_curve_data = {x:{'F':bin_size*x, 'rws':[]}}	# let's just start with it crayon simple.
	
	for i,rw in enumerate(roc_data):
		bin_num = int(rw['F']/bin_size)
		roc_curve_data[bin_num]+=[i]
	#
	roc_X = numpy.arange(0., 1.+bin_size, bin_size)
	#roc_Y = [None] + [len(x) for j,x in enumerate(roc_curve_data)]
	roc_Y = [None] + [max([roc_data['H'][i] for i in rw]) if len(rw)>0 else None for j,rw in enumerate(roc_curve_data)]
	#
	#roc_X, roc_Y = zip(*[[x,y] for x,y in zip(roc_X, roc_Y) if not (numpy.isnan(x) or numpy.isnan(y))])
	#
	plt.figure(fignum)
	if do_clf: plt.clf()
	if label_str!=None:
		curve_1 = plt.plot(roc_X, roc_Y, markers, zorder=1)
	else:
		curve_1 = plt.plot(roc_X, roc_Y, markers, zorder=1, label=label_str)
	# 
	plt.plot([0., 1.], [0., 1.], 'r-', lw=2, zorder=0)
	plt.xlabel('False alarm rate $F$')
	plt.ylabel('Hit rate $H$')
	
	plt.figure(fignum+1)
	if do_clf: plt.clf()
	if label_str!=None:
		plt.plot(roc_data['F'], roc_data['H'], '.')
	else:
		plt.plot(roc_data['F'], roc_data['H'], '.', label=label_str)		# ... and there's a better way to do this using plt.legend()...
	plt.plot([0., 1.], [0., 1.], 'r-', lw=2, zorder=0)
	plt.xlabel('False alarm rate $F$')
	plt.ylabel('Hit rate $H$')
		
	#return roc_data	
	return roc_curve_data
#
def plot_ROC(roc_array=None, section_id=16, fignum=0, num_points=100, do_clf=True, file_format='dumps/roc_detail_data/roc_sec_%d_nits_2500_allprams.npy', marker_str='.-', label_str='data'):
	#
	# roc_array should be from get_ROC_from_optimizer(), or equivalent. we'll want to thoroughly trap the input eventually.
	# see below, we'll want a recarray (or you could get away with a dict of lists like {'F':[], 'H':[]}
	#
	if roc_array==None and section_id!=None: roc_array=get_ROC_from_optimizer(section_id=section_id, file_format=file_format)
	if isinstance(roc_array[0], dict) or isinstance(roc_array,str): roc_array=get_ROC_from_optimizer(opt_data=roc_array)
	#
	# for now, just plot the ROC. the idea is to allow custom plotting -- plot multiple ROCs (and other stuff) on the same figure.
	#
	plt.figure(fignum)
	if do_clf: plt.clf()
	#
	plt.plot(roc_array['F'], roc_array['H'], marker_str, label=label_str)
	#
	return roc_array
#
def randomize_CFF(CFF_in=None, section_id=None):
	R=random.Random()
	if CFF_in==None and section_id!=None:
		CFF_in = vc_parser.combine_section_CFFs(sec_id)
	#
	CFF_rand=CFF_in.copy()
	#
	intervals = (CFF_in['event_year'][1:] - CFF_in['event_year'][0:-1]).tolist()
	intervals.sort(key=lambda x:R.random())
	random.shuffle(CFF_rand)
	CFF_rand['event_year'][0]=0.0
	for j,rw in enumerate(CFF_rand[1:]):
		CFF_rand[j+1]['event_year'] = CFF_rand[j]['event_year'] + intervals[j]
	#
	return CFF_rand
#
def random_like_CFF(CFF_in=None, section_id=None):
	# a different randomizer. instead of just mixing up a given sequence, create a random seqeunce 'like' the input.
	# in this case, it should have the same length and span approximately the same time-frame.
	R=random.Random()
	if CFF_in==None and section_id!=None:
		CFF_in = combine_section_CFFs(section_id)
	#
	mean_interval = (max(CFF_in['event_year'])-min(CFF_in['event_year']))/float(len(CFF_in))
	#
	CFF_rand=CFF_in.copy()
	CFF_rand[0]['event_year']=min(CFF_in['event_year'])
	CFF_rand[0]['event_number']=0
	#
	for j,rw in enumerate(CFF_rand[1:]):
		CFF_rand['event_year'][j+1]=CFF_rand[j]['event_year'] + mean_interval*R.random()
		CFF_rand['event_number'][j+1]=j+1
	#
	return CFF_rand
	
	#
def random_like_CFF_2(CFF_in=None, section_id=None):
	# a different randomizer. instead of just mixing up a given sequence, create a random seqeunce 'like' the input.
	# in this case, it should have the same length and span approximately the same time-frame.
	R=random.Random()
	if CFF_in==None and section_id!=None:
		CFF_in = combine_section_CFFs(section_id)
	#
	mean_interval = (max(CFF_in['event_year'])-min(CFF_in['event_year']))/float(len(CFF_in))
	mags = sorted(CFF_in['event_magnitude'].tolist(), key=lambda x: R.random())
	#return mags, CFF_in['event_magnitude'].tolist()
	times =[min(CFF_in['event_year'])]
	#
	for j, m in enumerate(mags[1:]):
		times += [times[-1]+R.random()*mean_interval]
	#
	print len(times), len(mags)
	#
	ary_out = numpy.core.records.fromarrays([range(len(mags)), times, mags ], names=['event_number', 'event_year', 'event_magnitude'], formats = ['int', 'float', 'float'])
	#
	return ary_out
	
	#
#
def get_ROC_from_optimizer(opt_data=None, section_id=16, file_format='dumps/roc_detail_data/roc_sec_%d_nits_2500_allprams.npy'):
	'''
	# modular ROC plotting:
	# input optimization output, return a rec_array ready for ROC plotting.
	#
	# provide opt_data as a string --> file name of optimization data set, list-o-dictionaries, each dict. is an optimization instance
	# (aka, from a monte-carlo method). (aka, the full-detail return from simple_metric_optimizer() OR,
	# opt_data=None; use section_id + file_format to load an optimization output file. see first couple-o-lines below.
	'''
	if isinstance(opt_data, str):
		opt_data=numpy.load(opt_data)
	#
	if opt_data==None and section_id!=None:
		opt_data = numpy.load(file_format % section_id)
	#
	type_dict={type('abc').__name__:'S128'}
	# let's go ahead and compile everything here. we might use this to analyize the b_0 and nyquist_factor dependencies as well.
	# first, make a list, then a recarray... (eventually, we want to reorganize this so we can allow a list of dicts. or recarray seamlessly).
	col_names = opt_data[0].keys() + [ky for ky in ['H', 'F', 'score_lin', 'score_geom'] if ky not in opt_data[0].keys()]
	#
	#roc_datas = [[val for key,val in rw.iteritems()] for rw in opt_data]
	roc_datas = []
	for i,rw in enumerate(opt_data):
		new_row = [val for key,val in rw.iteritems()]
		H = float(rw['n_predicted'])/(float(rw['n_predicted']) + rw['n_missed'])
		if rw['total_time']!=0.:
			F = rw['total_alert_time']/rw['total_time']
		else:
			F=None
		new_row += [H,F,H-F, H/F]
		#
		roc_datas += [new_row]

	#
	formats = [type(x).__name__ if not isinstance(x,str) else 'S128' for x in roc_datas[0]]
	
	# numpy.core.records.fromarrays(zip(*total_scores), names=['b_0', 'nyquist_factor', 'score', 'stdev'], formats = [type(x).__name__ for x in total_scores[0]])
	roc_datas = numpy.core.records.fromarrays(zip(*roc_datas), names=col_names, formats=formats)
	#
	
	
	return roc_datas
	
#
def plot_ROC_optimized_prams(roc_data=None, fignum=0, section_id=None, nits=1000):
	'''
	# make some plots for ROC optimization. roc_data should be the return value from something like:
	# opts_124=vc_parser.simple_metric_optimizer(nits=100, b_min=-.15, b_max=.25, keep_set=True, section_id=124)[1]
	# where the [1] gives the results from the full space explored by the MC simulation (and it's pretty huge since it includes
	# all the "alerts" sections).
	# simple_metric_optimizer() returns a list: [{dict of best prams}, {all prams}]
	# simple_metric_optimizer()[1] returns a list of dict objects like: 
	{'F': 0.26949814476171824,
 'H': 0.8216216216216217,
 'b': 0.19840002124380476,
 'n_missed': 66,
 'n_predicted': 304,
 'nyquist_factor': 0.5959081323587745,
 'score': 0.55212347685990348,
 'total_alert_time': 13372.468020446951,
 'total_time': 49619.888969070511}
	#
	#
	# roc_data: optimization data, namely the full parameter space output from simple_metric_optimizer(keep_set=True) {plus a bunch of
	# other inputs).
	'''
	#
	#if isinstance(roc_data, int): roc_data = 'interval_metric_output/roc_interval_metric_%d_LT.npy' % roc_data
	if isinstance(roc_data, str): roc_data=numpy.load(roc_data)
	if roc_data==None and section_id!=None:
		roc_data = simple_metric_optimizer(section_id=section_id, keep_set=True, m0=7.0, b_min=-.1, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, nits=nits, f_gt_lt=operator.lt)[1]
	#
	b_vals = [rw['b'] for rw in roc_data]
	alpha_vals = [rw['nyquist_factor'] for rw in roc_data]
	scores = [rw['score'] for rw in roc_data]		# ... or we might calc our own score, just to be sure.
	H,F = [], []
	scores_lin  = []
	scores_geom = []
	#
	for rw in roc_data:
		H += [float(rw['n_predicted'])/(rw['n_predicted'] + rw['n_missed'])]
		F += [rw['total_alert_time']/rw['total_time']]
		scores_lin += [H[-1]-F[-1]]
		scores_geom += [float(H[-1])/F[-1]]
	#
	f = plt.figure(fignum)
	plt.clf()
	plt.plot(F,H, '.')
	plt.plot(numpy.arange(2.), numpy.arange(2.), 'r-', zorder=0)
	plt.xlabel('$F$')
	plt.ylabel('$H$')
	#
	f=plt.figure(fignum+1)
	plt.clf()
	ax3d = f.add_subplot(111, projection='3d')
	ax3d.plot(b_vals, alpha_vals, scores_lin, '.')
	#for j,x in enumerate(b_vals):
	#	ax3d.plot([b_vals[j], b_vals[j]], [alpha_vals[j], alpha_vals[j]], [0., scores_lin[j]], '-o', lw=2)
	ax3d.set_xlabel('threshold slope, $b_0$')
	ax3d.set_ylabel('nyquist_factor')
	#ax3d.set_zlabel('ROC metric, (Percent Predicted) - (false alarm)')
	ax3d.set_zlabel('ROC metric, (H-F)')
	plt.title('Linear ROC metric')
	#
	# can i get a PCA analysis from this?
	#from  matplotlib.mlab import PCA
	#my_pca = PCA(numpy.array(zip(b_vals, alpha_vals)))
	#print "pca: ", my_pca
	my_pca = pca_tools.PCA_transform(zip(b_vals, alpha_vals))
	#
	print "pca axes (eigen vector matrix): ", my_pca.eig_vals, my_pca.eig_vecs
	# ... and it looks like this approach is not super conclusive. typically, you'd track along the largest axis. in this case, i think we'd track along the scaled diagonal
	# between them. so maybe try a [x,x] projection...
	# now, the idea will be to map a vector along the primary axis to nyquist_factor - b_0 space and show that we can turn this into a more traditional ROC analysis.
	min_xy = [min(b_vals), min(alpha_vals)]
	vx = numpy.array([my_pca.mus, my_pca.to_PCA([.5,0.])])
	vy = numpy.array([my_pca.mus, my_pca.to_PCA([0.,.5])])
	min_z = min(scores_lin)
	#ax3d.plot(zip(*vx)[0], zip(*vx)[1], [min_z for z in vx], color='r', lw=1.5, ls='-', marker='s')
	#ax3d.plot(zip(*vy)[0], zip(*vy)[1], [min_z for z in vx], color='g', lw=1.5, ls='-', marker='s')
	trace = numpy.arange(-.6, .6, .1)
	trace = zip(trace, trace)
	trace_pca = [list(my_pca.to_PCA(tr))+[min_z] for tr in trace]
	ax3d.plot(*zip(*trace_pca), color='r', lw=1.5, ls='-', marker='s', ms=5, alpha=.5)
	#
	# now, find (approximately) where this would end up in the ROC metric (H-F) point-cloud. get (b_0, alpha) from (x', y'), then find NN from MC set.
	# these data sets should be fairly small, so we could be sloppy about this (aka, just do a loop-loop). but instead, we'll use scipy.spatial.KDTree
	# or cKDTree, which i assume is the c implementation, or core class or something. more research and/or some benchmarking is needed.
	#
	K = cKDTree(numpy.array(zip(b_vals,alpha_vals)))
	NNpts = K.query(numpy.array(zip(*zip(*trace_pca)[0:2])), k=1, eps=0., p=2.)
	ax3d.plot(*zip(*[[b_vals[j], alpha_vals[j], scores_lin[j]] for j in NNpts[1]]), marker='.', color='r')
	ax3d.plot(*zip(*[[b_vals[j], alpha_vals[j], scores_lin[j]] for j in NNpts[1]]), marker='.', color='r')
	ax3d.plot(*zip(*[[b_vals[j], alpha_vals[j], min_z] for j in NNpts[1]]), marker='.', color='r')
	#
	###############################
	#
	f2=plt.figure(fignum+2)
	plt.clf()
	ax3d2 = f2.add_subplot(111, projection='3d')
	ax3d2.plot(b_vals, alpha_vals, scores_geom, '.')
	ax3d2.set_xlabel('threshold slope, $b_0$')
	ax3d2.set_ylabel('nyquist_factor')
	#ax3d.set_zlabel('ROC metric, (Percent Predicted) - (false alarm)')
	ax3d2.set_zlabel('ROC metric, (H/F)')
	plt.title('Geom. ROC metric')
	#
	f3=plt.figure(fignum+3)
	plt.clf()
	ax3d2 = f3.add_subplot(111, projection='3d')
	ax3d2.plot(F, b_vals, H, '.')
	ax3d2.plot([0.,1.,1.], [0., 1., 1.], [0., 1., 1.], 'r-', lw=2)
	ax3d2.set_xlabel('$F$')
	ax3d2.set_ylabel('$b$')
	ax3d2.set_zlabel('$H$')
	#
	f4=plt.figure(fignum+4)
	plt.clf()
	ax3d2 = f4.add_subplot(111, projection='3d')
	ax3d2.plot(F, alpha_vals, H, '.')
	ax3d2.plot([0.,1.,1.], [0., 1., 1.], [0., 1., 1.], 'r-', lw=2)
	ax3d2.set_xlabel('$F$')
	ax3d2.set_ylabel('nyquist factor')
	ax3d2.set_zlabel('$H$')
	#
	return my_pca

#
def plot_aggregate_metric(scores_in, n_top=None):
	# make a pretty (3d) plot of the aggregate type optimizer solutions.
	# scores_in are a recarray or dict like {'b_0':[], 'nyquist_factor':[], 'score':[]}
	if isinstance(scores_in, str):
		scores_in = numpy.load(scores_in)
	#scores_in.sort(order=('b_0', 'nyquist_factor'))
	#
	#
	#
	f2 = plt.figure(3)
	f2.clf()
	ax3d = f2.add_subplot(111, projection='3d')
	ax3d.plot(scores_in['b_0'], scores_in['nyquist_factor'], scores_in['score'], '.')
	ax3d.set_xlabel('threshold slope, $b_0$')
	ax3d.set_ylabel('nyquist_factor')
	#ax3d.set_zlabel('ROC metric, (Percent Predicted) - (false alarm)')
	ax3d.set_zlabel('ROC metric, (H-F)')
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
def plot_raw_ROC(scores_in=None, fignum=0, do_clf=True, plot_f_out=None):
	# a (simplified?) ROC plotter for optimized (segment partitioned) sets. this should be pretty generalized so we pass
	# a list of dicts (might trap other inputs as well later) and plot them. we'll sort out legends and stuff later.
	# ... and maybe remove the 'best', since we only plot here; there is no fitting... and so maybe we should move this to the paper_figs
	# script...
	# inputs: scores_in: numpy array or list-o-dictionaries with ROC fit data, aka, the results from an optimization. can also be a string
	# file-name to a .npy dump of said list-o-dictionaries.
	# produces a figure of H(F) for each optimization row in the input.
	#
	# here's the basic cascade of scores_in inputs:
	if scores_in==None:
		#scores_in = 'dumps/optimize_faultwise_best_scores.npy' # (i think).
		scores_in = 'dumps/gji_roc_lt_500/roc_sec_lt_18_nits_1000_allprams.npy'	
		
	if isinstance(scores_in, str):
		scores_in = numpy.load(scores_in)
	#
	# did we get a list of dicts instead of a recarray?
	if isinstance(scores_in[0], dict):
		# this should work, but it has not yet been tested... but it probably won't work because of string types...
		scores_in = numpy.rec.array([[val for val in rw.itervalues() if not hasattr(val, '__len__')] for rw in scores_in], names=[key for key,val in scores_in[0].iteritems() if not hasattr(val, '__len__')], formats = [type(x).__name__ for x in scores_in[0].itervalues() if not hasattr(x, '__len__')])
		#
	#
	F = [scores_in['total_alert_time'][i]/t for i, t in enumerate(scores_in['total_time'])]
	H = [float(N)/(float(N)+scores_in['n_missed'][i]) for i, N in enumerate(scores_in['n_predicted'])]
	#print "key:%s, <F>=%f, <H>=%f" % (key, numpy.mean(X), numpy.mean(Y))
	#print numpy.mean(b_col), numpy.mean(nq_col)
	#
	# calculate the total/mean score. we'll just calc the score since we can't necesssarily depend on the scoring
	# metric to remain constant.
	#mean_score = numpy.mean(numpy.array(H)/numpy.array(F))		# this might not be right; note we allow for different scoring functions.
	mean_score = numpy.mean(scores_in['score'])				# is 'score' part of the data set?
	#
	plt.figure(fignum)
	if do_clf:
		plt.clf()
		plt.plot([0., 1.], [0., 1.], 'r-', lw=2, alpha=.8, zorder=0)
	#plt.plot(F,H, label='%s: $<H/F>=%.3f$' % (key, mean_score), **datas['plot_kwargs'])
	#plt.plot(numpy.mean(F), numpy.mean(H), '*', ms=18, alpha=.75, color=datas['plot_kwargs']['color'], zorder=datas['plot_kwargs']['zorder'])
	plt.plot(F,H, 'o', label='%s: $<H/F>=%.3f$' % ('F', mean_score) )
	plt.plot(numpy.mean(F), numpy.mean(H), 'r*', ms=18, alpha=.75)
	#

def plot_best_opt_prams(scores_in=None, plot_f_out=None):
	'''
	# (this can probably be depricated and removed...)
	# plots from faultwise optimization (aka, from: optimize_metric_faultwise() )
	# in particular, plots an ROC diagram for  several populations of faultwise forecast fits.
	#
	# note that scores_in is a recarray with keys ['b_0'], ['nyquist_factor'], 'total_alert_time', 'total_time'
	# 'n_missed', 'n_predicted'. let's also add a script to handle a list of dicts...
	#
	'''
	# a sloppy tweak: this maps one col-name to another col name, aka b-->b_0
	col_name_subs = {}
	#
	if scores_in==None:
		#scores_in = 'dumps/optimize_faultwise_best_scores.npy' # (i think).
		scores_in = 'dumps/optimize_faultwise_trig_105_best_scores.npy'	# which i think is the same as the above, but more nits.
																		# maybe we want to compare a few sets? also get a "default" set (see below).
	if isinstance(scores_in, str):
		scores_in = numpy.load(scores_in)
	#
	# did we get a list of dicts instead of a recarray?
	if isinstance(scores_in[0], dict):
		# this should work, but it has not yet been tested... but it probably won't work because of string types...
		scores_in = numpy.rec.array([[val for val in rw.itervalues() if not hasattr(val, '__len__')] for rw in scores_in], names=[key for key,val in scores_in[0].iteritems() if not hasattr(val, '__len__')], formats = [type(x).__name__ for x in scores_in[0].itervalues() if not hasattr(x, '__len__')])
		#
		# sloppy tweak:
		if 'b_0' not in (scores_in.dtype.names): col_name_subs['b_0']='b'
		# 
	#
	# get mean b_0, nyquist_factor values:
	#mean_best_scores = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=numpy.mean(scores_in['b_0']),
	mean_best_scores = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=numpy.mean(scores_in[col_name_subs.get('b_0', 'b_0')]), nyquist_factor=numpy.mean(scores_in['nyquist_factor']), do_spp=False, do_plot=False, do_clf=True, n_cpus=None)
	col_names   = [key for key,val in mean_best_scores[0].iteritems() if not hasattr(val, '__len__')]
	col_formats = [type(val).__name__ for val in mean_best_scores[0].itervalues() if not hasattr(val, '__len__')]
	mean_best_scores = numpy.rec.array([ [val for val in rw.itervalues() if not hasattr(val, '__len__')] for rw in mean_best_scores], names=col_names, formats=col_formats)
	#
	# just comment this out so we can run Napa stuffp; we'll need to put it back/modify later.
	mean_best_regional_scores = None
	'''
	#####################
	# mean aggregate. 
	#####################
	#note this was run using a MPP process and not cpu-wise aggregated. just take the mean values. for now, hardcode the file.
	# so note this is not actually aggregating like the other sets. the other sets aggregate over fits to individual fault-catalogs.
	# this aggregates over 4 fits to the same regional catalog.
	mean_aggregate_emc_fits = numpy.load('dumps/optimize_fullemc_forecast_4_x_2500.npy')
	a_maf = numpy.mean([float(rw['n_predicted']) for rw in mean_aggregate_emc_fits])
	b_maf = numpy.mean([float(rw['total_alert_time']) for rw in mean_aggregate_emc_fits])
	c_maf = numpy.mean([float(rw['n_missed']) for rw in mean_aggregate_emc_fits])
	d_maf = numpy.mean([float(rw['total_time']) for rw in mean_aggregate_emc_fits])
	H_maf = a_maf/(a_maf + c_maf)
	F_maf = b_maf/d_maf
	score_maf = H_maf/F_maf
	b_0_maf = numpy.mean([rw['b'] for rw in mean_aggregate_emc_fits])
	nyquist_fact_maf = numpy.mean([rw['nyquist_factor'] for rw in mean_aggregate_emc_fits])
	#
	#print "mean_aggregate: b_0=%f, n_f=%f, H=%f, F=%f" % (b_maf, nyquist_fact_maf, H_maf, F_maf)
	
	#
	mean_best_regional_scores = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=b_0_maf, nyquist_factor=nyquist_fact_maf, do_spp=False, do_plot=False, do_clf=True, n_cpus=None)
	col_names   = [key for key,val in mean_best_regional_scores[0].iteritems() if not hasattr(val, '__len__')]
	col_formats = [type(val).__name__ for val in mean_best_regional_scores[0].itervalues() if not hasattr(val, '__len__')]
	mean_best_regional_scores = numpy.rec.array([ [val for val in rw.itervalues() if not hasattr(val, '__len__')] for rw in mean_best_regional_scores], names=col_names, formats=col_formats)
	######
	'''
	#
	default_scores = plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=0.0, nyquist_factor=.5, do_spp=False, do_plot=False, do_clf=True, n_cpus=None)
	# and we don't really care about fault name, etc. (we're plotting the points unnamed). make it a recarray (and copy this up above later...):
	col_names   = [key for key,val in default_scores[0].iteritems() if not hasattr(val, '__len__')]
	col_formats = [type(val).__name__ for val in default_scores[0].itervalues() if not hasattr(val, '__len__')]
	default_scores = numpy.rec.array([ [val for val in rw.itervalues() if not hasattr(val, '__len__')] for rw in default_scores], names=col_names, formats=col_formats)
	
	#
	#b_col = scores_in['b_0']
	b_col = scores_in[col_name_subs.get('b_0', 'b_0')]
	nq_col = scores_in['nyquist_factor']
	#
	# information gain:
	plt.figure(0)
	plt.clf()
	plt.plot([0., 1.], [0.,1.], '-', lw=2.5, label='Random: $<y/x>=1.0$')
	#
	data_sets =  {'Default':{'data':default_scores, 'plot_kwargs':{'zorder':5, 'marker':'o', 'color':'g', 'ms':8, 'ls':'', 'alpha':.7}},
	'Regional-mean optimized':{'data':mean_best_regional_scores, 'plot_kwargs':{'zorder':3, 'marker':(3,1,0.), 'color':'b', 'ms':10, 'ls':'', 'alpha':.7}},
	'Faultwise-mean optimized':{'data':mean_best_scores, 'plot_kwargs':{'zorder':4, 'marker':(2,1,0.), 'color':'r', 'ms':10, 'ls':'', 'alpha':.7}}, 
	'Faultwise optimized':{'data':scores_in, 'plot_kwargs':{'zorder':6, 'marker':(4,1,45.), 'color':'m', 'ms':14, 'ls':''}}
	}
	#
	for key, datas in data_sets.iteritems():
		# note the marker tuple: (num_sides/points, style, rotation_angle). styles: {0:polygon, 1:star-like, 2:asterisk/plus, 3:circle}
		scores = datas['data']
		if datas==None: continue
		#
		X = [scores['total_alert_time'][i]/t for i, t in enumerate(scores['total_time'])]
		Y = [float(N)/(float(N)+scores['n_missed'][i]) for i, N in enumerate(scores['n_predicted'])]
		#print "key:%s, <F>=%f, <H>=%f" % (key, numpy.mean(X), numpy.mean(Y))
		#print numpy.mean(b_col), numpy.mean(nq_col)
		#
		# calculate the total/mean score. we'll just calc the score since we can't necesssarily depend on the scoring
		# metric to remain constant.
		mean_score = numpy.mean(numpy.array(Y)/numpy.array(X))
		plt.plot(X,Y, label='%s: $<H/F>=%.3f$' % (key, mean_score), **datas['plot_kwargs'])
		plt.plot(numpy.mean(X), numpy.mean(Y), '*', ms=18, alpha=.75, color=datas['plot_kwargs']['color'], zorder=datas['plot_kwargs']['zorder'])
	# and the regional aggregate score:
	#plt.plot(F_maf, H_maf, 'r*', ms=14, zorder=7, alpha=.8, label='EMC regional')
	#plt.plot(F_maf, H_maf, 'k*', ms=17, zorder=6, alpha=.8)
	plt.legend(loc='lower right', numpoints=1)
	
	plt.xlabel('percent alert time, "false alarm" rate $F$')
	plt.ylabel('percent predicted,  "hit" rate $H$')
	#
	if plot_f_out != None:
		plt.savefig(plot_f_out)
	#
	# best-fit parameters:
	plt.figure(1)
	plt.clf()
	ax=plt.gca()
	X=scores_in[col_name_subs.get('b_0', 'b_0')].copy()
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
	ax2.set_ylabel('score')
	#
	f2 = plt.figure(3)
	f2.clf()
	ax3d = f2.add_subplot(111, projection='3d')
	ax3d.plot(scores_in[col_name_subs.get('b_0', 'b_0')], scores_in['nyquist_factor'], scores_in['score'], '.')
	ax3d.set_xlabel('b_0')
	ax3d.set_ylabel('nyquist_factor')
	ax3d.set_zlabel('score')
	#
	mean_b = numpy.mean(scores_in[col_name_subs.get('b_0', 'b_0')])
	std_b  = numpy.std(scores_in['score'])
	#
	mean_nf = numpy.mean(scores_in['nyquist_factor'])
	std_nf = numpy.std(scores_in['nyquist_factor'])
	#
	print "mean_b: %f +/- %f" % (mean_b, std_b)
	print "mean_nf:    %f +/- %f" % (mean_nf, std_nf)


#
##############################################
# end ROC and predictability bits
##################################
#
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
def recurrence_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, keep_figs=False, output_dir='CDF_figs'):
	'''
	# (these should be the non-conditional weibull plots, equivalent to cond_wiebul(t0=0).
	# for each section_id, fetch the mean_recurrence data.
	# 1) plot each cumulative probability
	# 2) including a weibull fit.
	# 3) plot each cum. prob onto a single figure
	# 4) a cumulative figure
	'''
	#
	if section_ids in (None, [], ()): section_ids = list(emc_section_filter['filter'])
	section_ids=list(section_ids)
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
		# faultwise mean/composite:
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
		if keep_figs==False and sec_id!=-2:
			plt.close(i)
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
##################################################
# development and diagnostic bits:
#
def non_converging_weibul_example(fnum0=7, nits=100000, n_cpus=None):
	'''
	# some... many of the faultwise Weibull fits (aka, fits to "waiting time", or hazard function, for a single fault segment
	# do not converge. an example is Section 25, t0=513.073 (t0 = <mean interval> * 2.0).
	# basically, this nominally appears to be converging (in the MC method), but at tau and beta --> 0, which is sort of
	# nonsensical.
	'''
	#
	# first, get the data:
	A=hazard_function_fitting_test(section_id=25)
	# times:
	t0=513.07251464973058
	T=A[t0]	# we should index these differently, but for now it works well enough.
	T.sort()
	Y=[float(k)/len(T) for k in xrange(len(T))]
	#
	plt.figure(fnum0)
	plt.clf()
	plt.plot(T, Y, '.-')
	#
	# and if you like, you can adjust the starting points, etc., but this probably blows up.
	try:
		fit, cov = spo.curve_fit(lambda x, beta, chi: f_weibull(x=x, beta=beta, chi=chi, x0=513.073), xdata=T, ydata=Y, p0=[1.5, 500.])
		print "converging fit worked! ", fit, cov
	except:
		print "fit broke, as expected..."
	#
	# now, try an MC method:
	print "try MC method... (%s)" % str(n_cpus)
	prams_dict = {'chi':[0., 650.], 'beta':[0.,6.], 'x0':[t0,t0]}
	#
	Zmin, Z=mcfitter(func=f_weibull, X=T, Y=Y, prams_dict=prams_dict, nits=nits, n_cpus=n_cpus)
	#
	'''
	# depricated: MPP functionality is now incorporated into mcfitter().
	if n_cpus==1:
		# simple SPP method:
		print "doing spp"
		Zmin, Z=mcfitter(func=f_weibull, X=T, Y=Y, prams_dict=prams_dict, nits=nits)
	else:
		print "doing mpp"
		if n_cpus==None: n_cpus = mpp.cpu_count()
		#
		# just use processes:
		results=[]
		pool = mpp.Pool(n_cpus)
		for i in xrange(n_cpus):
			results+=[pool.apply_async(mcfitter, kwds={'func':f_weibull, 'X':T, 'Y':Y, 'prams_dict':prams_dict, 'nits':nits/n_cpus, 'n_cpus':1})]
		pool.close()
		pool.join()
		#
		#Zmin, Z=[], []
		Zmin, Z = results[0].get()
		 
		for j, result in enumerate(results):
			if j==0: continue
			#
			a,b = result.get()
			Zmin = numpy.append(Zmin, a)
			Z    = numpy.append(Z, b)
			
			#Amin, A = result.get()
			#Zmin += Amin.tolist()
			#Z += A.tolist()
		#
		#Zmin = numpy.core.records.fromarrays(zip(*Zmin), dtype=A.dtype)
		#Z = numpy.core.records.fromarrays(zip(*Z), dtype=A.dtype)
		Zmin.sort(order='chi_sqr')
	'''
		
		
	print "mins: ", Zmin[-1]
	#
	f=plt.figure(fnum0+1)
	plt.clf()
	f.add_axes([.1, .1, .8, .8], projection='3d')
	#
	# some short-hand:
	plot_len = min(10000, len(Z))
	x = Z['beta'][0:plot_len]
	y = Z['chi'][0:plot_len]
	z = Z['chi_sqr'][0:plot_len]
	#
	ax3d = plt.gca()
	ax3d.scatter(x,y,z, c=z, cmap=cm.jet, alpha=.4)
	plt.plot(Zmin['beta'], Zmin['chi'], Zmin['chi_sqr'], 'r.-', alpha=.7)
	#ax3d.plot_trisurf(x,y,z, cmap=cm.jet, linewidth=.1)
	Zmin.sort(order=('beta', 'chi'))
	plt.plot(Zmin['beta'], Zmin['chi'], 'r.--', alpha=.7)
	#
	# contour plot?
	#triang = tri.Triangulation(x, y)
	triang = tri.Triangulation(Z['beta'], Z['chi'])
	little_tri = tri.Triangulation(x,y)\
	# ... and this does not appear to work in 3d; is tricontourf() supported by Axes3D?
	# it looks like we'd need to grid these data to plot the contours "behind" the 3D plot (on the axes planar surfaces).
	#ax3d.tricontourf(x, y, z, zdir='z', offset = 0., cmap=cm.jet)
	#ax3d.tricontourf(x, y, z, zdir='x', offset = -2., cmap=cm.jet)
	#ax3d.tricontourf(x, y, z, zdir='y', offset = -100., cmap=cm.jet)
	
	#plt.gca().set_aspect('equal')
	plt.figure(fnum0+4)
	plt.clf()
	#plt.tricontourf(triang, z)
	plt.tricontourf(triang, Z['chi_sqr'], 32)
	plt.plot(Zmin['beta'], Zmin['chi'], 'r.--', alpha=.9)
	plt.xlabel('$\\beta$')
	plt.ylabel('nyquist factor')
	#
	f=plt.figure(fnum0+2)
	plt.clf()
	Xbeta=Zmin['beta'].tolist()
	Xbeta.sort()
	plt.plot(Xbeta, Zmin['chi_sqr'], 'bo-')
	ax2=plt.gca().twiny()
	Xchi = Zmin['chi'].tolist()
	Xchi.sort()
	ax2.plot(Xchi, Zmin['chi_sqr'], 'rs-')
	#
	plt.figure(fnum0+3)
	plt.clf()
	plt.plot(T, Y, '.-')
	Xfit = [min(T) + x*(max(T)-min(T))/500.0 for x in xrange(500)]
	Yfit = [f_weibull(x=x, beta=Zmin[-1]['beta'], chi=Zmin[-1]['chi'], x0=t0) for x in Xfit]
	plt.plot(Xfit, Yfit, '--')
	#
	return Zmin, Z
#
def hazard_function_fitting_test(section_id=16, file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5]):
	'''
	# we're getting some screwy behavior trying to fit the hazard functions (waiting times). let's spell it all out here
	# and possibly introduce a MC fitting method.
	'''
	#
	lw=2.5
	ms=5.
	max_x = 0.
	min_x = 0.
	dT_composite_faultwise = []
	full_catalog = []
	# control color cycling:
	colors_ =  mpl.rcParams['axes.color_cycle']
	#
	ary_in = numpy.load(file_path_pattern % section_id)
	ary_m0 = [rw for rw in ary_in if rw['event_magnitude']>m0]
	ary_m0 = numpy.core.records.fromarrays(zip(*ary_m0), dtype=ary_in.dtype)
	#
	# first, get the recurrence intervals:
	delta_Ts = ary_m0['event_year'][1:] - ary_m0['event_year'][0:-1]
	mean_dt = numpy.mean(delta_Ts)
	t0s = [mean_dt*x for x in [0., .5, 1.0, 1.5, 2.0, 2.5]]
	print "t0s: ", t0s
	#
	time_serieses = {}
	#
	fnum=0
	for t0 in t0s:
		plt.figure(fnum)
		plt.clf()
		#
		X = [t for t in delta_Ts if t>=t0]
		X.sort()
		Y = [(1.0+k)/float(len(X)) for k in xrange(len(X))]
		#
		print "t0=%f, len(X): %d" % (t0, len(X))
		plt.plot(X, Y, '.-', label='$t_0=%.3f' % t0, ms=8, zorder=1)
		#
		try:
			fit_prams, fit_cov = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=t0), xdata=numpy.array(X), ydata=numpy.array(Y), p0=numpy.array([mean_dt, 1.5]))
			#
			dx = (max(X)-min(X))/500.
			fit_X = [min(X) + i*dx for i in xrange(500)]
			fit_Y = [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=t0) for x in fit_X]
			#	
			plt.plot(fit_X, fit_Y, '--', label='$\\beta=%.3f, \\chi=%.3f$' % (fit_prams[1], fit_prams[0]), lw=2, zorder=4, alpha=.7)
		except:
			print "fit for t0 = %.3f failed." % t0
		plt.legend(loc=0, numpoints=1)
		#
		fnum+=1
		#
		time_serieses[t0]=X
	return time_serieses
#
def tau_t0_fig(section_ids=None, glob_pattern='VC_CDF_WT_figs/VC_CDF_WT_*.npy'):
	# plot tau vs t0 and t0_index (which would be a fractional t0... how' bout vs t0 and ts t0/tau.
	# basically, we expect a break where the weibull distribution stops fitting. we will specifically
	# be seeing the difference beteewn converging (n_max -> 600) and MC fits (in this case, 100000 iterations).
	#
	# general take-away talking points: we see a change (fits break, beta --> big, etc. for \tau_r>1 t_0r>1).
	#
	ary_files = glob.glob(glob_pattern)
	#
	# keep track of the chi_0 and beta_0 fits in a dict.
	chi_beta_0s = {}
	#
	#full_array = numpy.load(ary_files[0])
	for j, fl in enumerate(ary_files):
		this_array = numpy.load(fl)
		#
		# it's possible we won't have a zero time point. just take the first one (but don't assume it's sorted):
		this_array.sort(order='t0')
		chi_beta_0s[this_array[0]['section_id']] = {'chi':this_array[0]['chi'], 'beta':this_array[0]['beta']}
		#
		# note "reduce()" syntax below.
		#if j==0:
		#	full_array = numpy.load(ary_files[0])
		#	continue
		#full_array = numpy.append(full_array, this_array)
	full_array = reduce(numpy.append, [numpy.load(x) for x in ary_files])		# apply numpy.append recursively through the list.
																				# short-hand for loop logic above.
	#
	#print len(full_array)
	t0s = full_array['t0']
	chis = full_array['chi']
	r_chis = [rw['chi']/chi_beta_0s[rw['section_id']]['chi'] for rw in full_array]
	r_t0s   = [rw['t0']/chi_beta_0s[rw['section_id']]['chi'] for rw in full_array]	# and these won't be exactly t0/2, t0, etc.
																					# bc the original t0 values were determined from 
																					# mean recurrence intervals, not chi/tau.
	betas = full_array['beta']
	#
	# now, we want to separate the MC from the convergent fits:
	MC_data = numpy.array([rw for rw in full_array if 'MC' in rw['fit_type']], dtype=full_array.dtype)
	fit_data = numpy.array([rw for rw in full_array if 'spo.' in rw['fit_type']], dtype=full_array.dtype)
	#
	t0s_mc = MC_data['t0']
	chis_mc = MC_data['chi']
	r_chis_mc = [rw['chi']/chi_beta_0s[rw['section_id']]['chi'] for rw in MC_data]
	r_t0s_mc   = [rw['t0']/chi_beta_0s[rw['section_id']]['chi'] for rw in MC_data]
	betas_mc = MC_data['beta']
	#
	t0s_fit = fit_data['t0']
	chis_fit = fit_data['chi']
	r_chis_fit = [rw['chi']/chi_beta_0s[rw['section_id']]['chi'] for rw in fit_data]
	r_t0s_fit   = [rw['t0']/chi_beta_0s[rw['section_id']]['chi'] for rw in fit_data]
	betas_fit = fit_data['beta']
	#
	#print "lens: %d, %d" % (len(MC_data), len(fit_data))
	#
	#colors = [(1 if 'MC' in x else 0) for x in full_array['fit_type']]
	
	#
	plt.figure(0)
	plt.clf()
	#
	ax=plt.gca()
	plt.plot(r_t0s_fit, chis_fit, 's', alpha=.5, zorder=1, color='b', label='spo.curve_fit')
	plt.plot(r_t0s_mc, chis_mc, 's', alpha=.5, zorder=1, color='g', label='MC')
	ax.legend(loc=0, numpoints=1)
	#
	plt.ylabel('$\\tau$')
	plt.xlabel('reduced $t_0$, $t_{0r} = t_0/\\tau_0$')
	ax = plt.gca().twinx()
	#
	ax.plot(r_t0s_fit, r_chis_fit, 'o', zorder=2, color='r', label='spo.curve_fit')
	ax.plot(r_t0s_mc, r_chis_mc, 'o', zorder=2, color='y', label='MC')
	ax.legend(loc=0, numpoints=1)
	#
	ax.set_ylabel('reduced $\\tau$, $\\tau_r = \\tau/\\tau_0$')
	ax.set_xlabel('reduced $t_0$, $t_{0r} = t_0/\\tau_0$')
	ax.legend(loc=0, numpoints=1)
	#
	#
	plt.figure(1)
	plt.clf()
	plt.plot(r_t0s_fit, betas_fit, 'bo', label='spo.curve_fit')
	plt.plot(r_t0s_mc, betas_mc, 'go', label='MC')
	ax=plt.gca()
	ax.set_xlabel('reduced $t_0$, $t_{0r} = t_0/\\tau_0$')
	ax.set_ylabel('$\\beta$')
	ax.legend(loc=0, numpoints=1)
	#
	plt.figure(2)
	plt.clf()
	plt.plot(r_chis_fit, betas_fit, 'bo', label='spo.curve_fit')
	plt.plot(r_chis_mc, betas_mc, 'go', label='MC')
	ax=plt.gca()
	ax.set_xlabel('reduced $\\tau$, $\\tau_r = \\tau/\\tau_0$')
	ax.set_ylabel('$\\beta$')
	ax.legend(loc=0, numpoints=1)
	#

	#
	return full_array
###################################
# end development and diagnostic bits
###############################
#
#######################################
# statistics figure scripts:
#
def expected_waiting_time_t0(section_ids=None, catalog=None, m0=7.0, fits_data_file_CDF='CDF_EMC_figs/VC_CDF_Weibull_fits_dump.npy', WT_catalog_format='data/VC_CFF_timeseries_section_%d.npy', sim_file=default_sim_file, n_t0=10000, fnum=0, do_local_fit=True):
	# (these really are EWT, post-nomenclature correction).
	# do_local_fit: fit the weibull distribution along the way (for each time-step)? obviously, we get a better fit, but what
	# does it mean? could be used for tabulated uncertainties. unfortunately, it seems to break down on the right hand side of the dist.
	# ... and work out the details later...
	# if we have a catalog, pass it in (preferably a recarray). otherwise, provide section_ids).
	# TODO:
	# n_t0: number of data points (number ot t0 values) ... and note, this is an error (in optimization anyway). we should limit this to just
	#   the set of {t0} for t0=t_{new_earthquake}, since there's no new information between earthquakes. so fix this later.
	#
	# note: for now,use a pre-compiled catalog. assembling event-section_id catalogs requires some parsing. these catalogs
	# have been queried and stored in VC/vc_parser/data as *.npy files. there is a script here somewhere to compile an event-section-wise
	# catalog; also checkout the quakelib.Event() bits (which i think are not included (yet) in the new quakelib).
	#
	# the various parameters are (probably) already calculated. we can look into recalculating them later; note that
	# the WT fits can get a bit ugly because some fits don't naturally converge, so we use a MC method... which is, of
	# course questionable because it does not seem to converge, but it does seem to give the best fit.
	#
	# fits_data_file_CDF: looks like we don't use this any longer.
	#
	if section_ids=='emc': section_ids = list(emc_sections) 	# globally defined list of emc_sections
	if isinstance(section_ids, int) or isinstance(section_ids,float): section_ids = [section_ids]
	section_ids = list(set(section_ids))
	#
	# load the pre-calced fits? if they don't exist, run them.
	# (but then, we never ended up using these...)
	'''
	try:
		# but it looks like we don't use this any longer...
		cdf_fits = numpy.load(fits_data_file_CDF)	# recarray with dtype: dtype=[('t0', '<f8'), ('section_id', '<i8'), ('chi', '<f8'), ('beta', '<f8'), ('sigma_chi', '<f8'), ('sigma_beta', '<f8'), ('chi_sqr', '<f8'), ('fit_type', 'S16')])
	except:
		# ... ok, this is the right idea, but it's complicated, so let's handle it externally for now.
		# anyway, we might not want to do it this way because it might run for a long, long time...
		# fits failed. run a new set:
		# waiting_time_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], keep_figs=False, output_dir='VC_CDF_WT_figs', mc_nits=100000, n_cpus=None)
		print "fits data not found. run waiting_time_figs()"
		return None
		
		#wt_temp = waiting_time_figs(section_ids=section_ids, m0=m0, output_dir='temp_weibul_fits', mcnits=2000)
	'''	
	#
	# we want the expected \Delta t to the next "big one" as a function of t_0 (aka, <Delta t> (t_0) ), and in this case
	# t_0 is basically current ellapsed time. nominally, we shoul do this 1) directly from data, 2) using t0=0 fits, 3)
	# custom fits for each value of t0 (or at least use prams[t_0 closest to t_0].
	#
	# we can pull waiting time data from the sim-file, or we can use the pre-calculated bits...
	# (see: /home/myoder/Documents/Research/Yoder/VC/vc_parser/VC_CDF_WT_figs/*.npy).
	# for now, let's pull the data from ths sim-file, but we'll use pre-calculated fit data (maybe).
	#
	#
	if catalog==None:
		# (providing a catalog will be the exception (maybe during development processes, etc.). generally, load
		# fron file.npy.
		#
		#catalog=[]
		#with h5py.File(sim_file, 'r') as vc_data:
		#	# ... and here it is... this is a bad idea; the whole point of those pre-calculated events is
		#	# that it is difficult and expensive to grab section -- event combination data (aka, event_table does not
		#	# know if the event occurred on a valid section, etc. use the pre-calculated numpy.dump() files.
		#	vc_events = vc_data['event_table']
		#	#catalog = [rw for rw in vc_events if (rw['event_magnitude']>=m0 and rw['section_id
		#	pass
		#
		
		#str_mag = str(m0).replace('.','')
		
		# this combine_section_CFFs() bit needs to be checked out... also, a start_year should be added to dump the early-sim junk.
		#
		#catalog = numpy.load(WT_catalog_format % (section_ids[0]))
		#for j, section_id in enumerate(section_ids[1:]):
		#	catalog=numpy.append(catalog, numpy.load(WT_catalog_format % (section_id)))
		#	#print "(%d) appending catalog: %d" % (len(catalog), section_id)
		catalog = combine_section_CFFs(section_ids, start_year=None, end_year=None)
		#
		#return catalog
		catalog.sort(order='event_year')
		#
		# nominally, we should perform a set() operation to be sure that we don't have any dupes (from compiling these catalogs),
		# but preliminary analysis suggests that we don't.
		#
	#
	event_years = [rw['event_year'] for rw in catalog if rw['event_magnitude']>=m0]
	recurrence_intervals = [event_years[i]-event_years[i-1] for i in xrange(1,len(event_years))]
	#
	recurrence_intervals.sort()
	delta_t_total = 2.0*max(recurrence_intervals)
	#
	expected_delta_ts = []	# will be like [ [{.25}, {.5}, {.75}], ... ]
	model_delta_ts_0 = []	# from base weibull fit (t0=0)
	model_delta_ts_1 = []	# from the weibull fitting as we go...
	#
	dt0 = delta_t_total/float(n_t0)
	this_t0 = 0.0
	#
	# preliminary fit:
	#
	X=recurrence_intervals
	N_ri = len(X)
	Y=numpy.array([j/float(N_ri) for j in numpy.arange(1., N_ri+1, 1.)])
	#
	
	#return [X,Y, numpy.array([numpy.mean(X), 1.5])]
	
	fit_prams = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=0.), xdata=numpy.array(X), ydata=Y, p0=numpy.array([numpy.mean(X), 1.5]))
	
	#fit_prams_0, fit_cov_0 = spo.curve_fit(f_weibull, xdata=numpy.array(X), ydata=numpy.array([j/float(len(X)) for j in numpy.arange(1., len(X)+1)]), p0=numpy.array([mean_dT, 1.5]))
	fit_tau_0 = fit_prams[0][0]
	fit_beta_0 = fit_prams[0][1]
	#
	#print "fit_prams: ", fit_prams
	fit_prams_dyn = fit_prams[0].copy()
	fit_tau_1 = fit_tau_0
	fit_beta_1 = fit_beta_0
	#prev_fit_prams = fit_prams_dyn.copy()
	#
	# this would be a lot faster if we advanced by event rather than d_t0...
	while this_t0<delta_t_total:
		these_delta_ts = [dt for dt in recurrence_intervals if dt>=this_t0]
		these_delta_ts.sort()	# i think this is now redundant...
		#
		# now, get medianses:
		# (there's a smart recursive way to do this, but we're just goint to code it):
		N = len(these_delta_ts)
		med_25 = numpy.median(these_delta_ts[:int(math.floor(N/2.))])
		med_5 = numpy.median(these_delta_ts)
		med_75 = numpy.median(these_delta_ts[int(math.ceil(N/2.)):])	# provide default values in the event
		#if numpy.isnan(med_75): med_75 = expected_delta_ts[-1][3]+this_t0							# a median can't be found.
		#if numpy.isnan(med_25): med_25=this_t0
		#
		#if this_t0>120.: print med_25, med_5, med_75
		expected_delta_ts+= [[this_t0, med_25-this_t0, med_5-this_t0, med_75-this_t0]]
		model_delta_ts_0 += [[this_t0] + [f_inv_weibull(P=p, tau=fit_tau_0, beta=fit_beta_0, t0=this_t0)-this_t0 for p in [.25, .5, .75]]]
		#
		# dynamic fit:
		if do_local_fit:		# eventually maybe get rid of this, but there are a few other "if"s that will have to be done.
			try:
			#if True:
				# try to get new fit parameters. if the fit fails, continue to use the previous values:
				N_local = len(these_delta_ts)
				Y=[j/float(N_local) for j in xrange(1, N_local+1)]
				#
				new_fit_prams_dyn = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=this_t0), xdata=numpy.array(these_delta_ts), ydata=numpy.array(Y), p0=numpy.array([numpy.mean(these_delta_ts), 1.5]))[0]
				#
				# if the fit fails, we'll kick over to exception handling at this point.
				#
				#prev_fit_prams = fit_prams_dyn.copy()
				fit_prams_dyn = new_fit_prams_dyn.copy()					# this will only execute if the fit succeeds. otherwise, 
				fit_tau_1 = new_fit_prams_dyn[0]
				fit_beta_1 = new_fit_prams_dyn[1]
				new_fit_prams_dyn = None							# we keep the old values (do we need to copy?)
			
			except:
				#print new_fit_prams_dyn
				#print "excepting: "
				#print prev_fit_prams
				#print fit_prams_dyn
				#
				#fit_tau_1  = fit_tau_0
				#fit_beta_1 = fit_beta_0
				#
				#fit_prams_dyn = prev_fit_prams.copy()
				#return None
			 	#fit failed. we could have a go at MC, but let's just punt and use the previous set -- aka, do nothing.
				# eventually, we might want to know that the fit failed, but for now just pass...
				pass
				#
				#fit_tau_1 =  (fit_tau_1 + fit_tau_0)/2.	# note: this will be a creaping average...
				#fit_beta_1 = (fit_beta_1 + fit_beta_0)/2.
		
			#model_delta_ts_1 += [[this_t0] + [f_inv_weibull(P=p, tau=fit_prams_dyn[0], beta=fit_prams_dyn[1], t0=this_t0)-this_t0 for p in [.25, .5, .75]]]
			model_delta_ts_1 += [[this_t0] + [f_inv_weibull(P=p, tau=fit_tau_1, beta=fit_beta_1, t0=this_t0)-this_t0 for p in [.25, .5, .75]]]
			#
		#
		#
		this_t0+=dt0
		#if this_t0>175.:break
	#
	expected_delta_ts = numpy.core.records.fromarrays(zip(*expected_delta_ts), names=['t0', 'med25', 'med50', 'med75'], formats = [type(x).__name__ for x in expected_delta_ts[0]])
	#return model_delta_ts_0
	model_delta_ts_0 = numpy.core.records.fromarrays(zip(*model_delta_ts_0), names=['t0', 'med25', 'med50', 'med75'], formats = [type(x).__name__ for x in model_delta_ts_0[0]])
	
	#
	plt.figure(fnum)
	plt.clf()
	xlim = max([t0 for j,t0 in enumerate(expected_delta_ts['t0']) if expected_delta_ts['med50'][j]>0])
	#plt.gca().set_xlim(right=max(expected_delta_ts['t0'])*1.25)
	plt.gca().set_xlim(right=xlim*1.25)
	plt.plot(expected_delta_ts['t0'], expected_delta_ts['med25'], 'b-', lw=1.5)
	plt.plot(expected_delta_ts['t0'], expected_delta_ts['med50'], 'k-', lw=2.5, alpha=.8)
	plt.plot(expected_delta_ts['t0'], expected_delta_ts['med75'], 'b-', lw=1.5)
	plt.fill_between(expected_delta_ts['t0'], expected_delta_ts['med25'], expected_delta_ts['med75'], color='y', alpha=.7)
	#
	# model_0:
	#plt.plot(model_delta_ts_0['t0'], model_delta_ts_0['med25'], 'b--', lw=2.5)
	#plt.plot(model_delta_ts_0['t0'], model_delta_ts_0['med50'], 'm--', lw=2.5)
	#plt.plot(model_delta_ts_0['t0'], model_delta_ts_0['med75'], 'b--', lw=2.5)
	#
	if do_local_fit:
		model_delta_ts_1 = numpy.core.records.fromarrays(zip(*model_delta_ts_1), names=['t0', 'med25', 'med50', 'med75'], formats = [type(x).__name__ for x in model_delta_ts_1[0]])
		# model_1:
		plt.plot(model_delta_ts_1['t0'], model_delta_ts_1['med25'], 'c-', lw=1.5, alpha=.6)
		plt.plot(model_delta_ts_1['t0'], model_delta_ts_1['med50'], 'm-', lw=1.5, alpha=.6)
		plt.plot(model_delta_ts_1['t0'], model_delta_ts_1['med75'], 'b-', lw=1.5, alpha=.6)
	
	return catalog
#
#def waiting_time_single_curve(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0 = 5.0, mc_nits=100000, n_cpus=None, fignum=None, sim_file=default_sim_file, output_type='dict'):
def conditional_RI_single_curve(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0 = 5.0, mc_nits=100000, n_cpus=None, fignum=None, sim_file=default_sim_file, output_type='dict'):
	# calculate (and plot if fignum!=None) a single waiting_time distribution. illustrate that, since single (or few)-fault sets are
	# not weibull-random, we can calculate specific probabilities for a set of faults.
	# this is P(t;t_0), right -- waiting-time probability distribution?
	#
	# some preliminary bits:
	max_x = 0.
	min_x = 0.
	#
	if section_ids==None or len(section_ids)==0:
		with h5py.File(sim_file) as vc_data:
			section_ids = set(vc_data['block_info_table']['section_id'])
		#
	#
	section_ids=list(section_ids)
	#
	# get a catalog of unique rows (section-wise catalogs can contain duplicate entries).
	catalog = combine_section_CFFs(sections=section_ids, ary_in_format=file_path_pattern, start_year=10000., end_year=None)
	catalog.sort(order='event_year')
	mean_rec_data = mean_recurrence(ary_in=catalog, m0=m0)
	X = mean_rec_data['dists']['dT'].tolist()
	#
	N = float(len(X))	
	X.sort()
	this_X = [x for x in X if (x-t0)>=0.]
	this_X.sort()		# just in case...
	N = float(len(this_X))
	Y = [float(j)/N for j in range(1, int(N)+1)]
	#
	# keep track of min/max X values for later global plot... or we could use them here as well.
	max_x = max(max_x, max(X))
	min_x = min(min_x, min(X))
	#
	mean_dT = numpy.mean(X)
	#
	#print "preliminary fit for section_id=%d" % sec_id
	fit_prams_0, fit_cov_0 = spo.curve_fit(f_weibull, xdata=numpy.array(X), ydata=numpy.array([j/float(len(X)) for j in numpy.arange(1., len(X)+1)]), p0=numpy.array([mean_dT, 1.5]))
	#
	try:
		# f_weibull(x=None, chi=1.0, beta=1.0, x0=None)
		fit_type = 'spo.curve_fit'	# try a converging fit; if it fails, we'll use a fit_type='MC' (monte carlo)
		fit_prams, fit_cov = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=t0), xdata=numpy.array(this_X), ydata=numpy.array(Y), p0=numpy.array([max(1.0, mean_dT-t0), 1.5]))
		#			
		print "fitted (converging)...", fit_cov
		pram_sigmas = numpy.sqrt(numpy.diag(fit_cov))		# the uncertainties (variances) of the fit parameters are the diagonals of the covariance matrix... right?
		mean_chi_sqr = numpy.mean([(f_weibull(x=this_X[k], chi=fit_prams[0], beta=fit_prams[1], x0=t0)-Y[k])**2. for k, xx in enumerate(this_X)]) # in xrange(len(X))])
		stdErr = mean_chi_sqr/math.sqrt(N)
		print 'cov: ', fit_cov
		print "fit_prams(%d/%s)[t0 = %f]: %s / %s" % (0, str(section_ids), t0, str(fit_prams), str(fit_prams_0))
		#print "fit_prams(%d/%s)[t0 = %f]: %s / s" % (i, str(section_ids), t0, str(fit_prams))
		#
		#best_fit_dict[sec_id][t0] = [sec_id, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr]
		print "assign bits..."
		fit_vals = [t0, section_ids, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr, fit_type]
	except:
		return None
		fit_type = 'MC_%d' % mc_nits
		# converging fit failed. do an MC method:
		print "converging fit failed. try an MC approach:"
		prams_dict = {'chi':[0., 2.5*mean_dT], 'beta':[0.,6.], 'x0':[t0,t0]}	# can we automatedly guess at these prams?
																			# maybe at some point, we should devise a quasi-
																			# converging MC method.
		Zmin, Z=mcfitter(func=f_weibull, X=this_X, Y=Y, prams_dict=prams_dict, nits=mc_nits, n_cpus=n_cpus)
		Zmin.sort(order='chi_sqr')
		best_fit = Zmin[0]
		#
		fit_prams = [best_fit['chi'], best_fit['beta']]
		sigma_chi  = (prams_dict['chi'][1]-prams_dict['chi'][0])/mc_nits
		sigma_beta = (prams_dict['beta'][1]-prams_dict['beta'][0])/mc_nits
		#
		#best_fit_dict[sec_id][t0] = [sec_id, best_fit['chi'], best_fit['beta'], sigma_chi, sigma_beta, best_fit['chi_sqr']]
		fit_vals = [t0, section_ids, best_fit['chi'], best_fit['beta'], sigma_chi, sigma_beta, best_fit['chi_sqr'], fit_type]
	#
	if output_type.lower()=='dict':
		fit_vals = {key:val for key,val in zip(*[['t0', 'section_ids', 'tau', 'beta', 'sigma_tau', 'sigma_beta', 'chi_sqr', 'fit_type'],fit_vals])}
	
	# note: "sigmas" might really be variances. curve_fit() return variances; if we want stdev, take the sqrt(var)...
	return fit_vals
#
def waiting_time_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], keep_figs=False, output_dir='VC_CDF_WT_figs', mc_nits=100000, n_cpus=None):
	'''
	# this function will be removed and replaced by its more modular cousin: conditional_RI_figs()
	# though note, the faultwise-aggregate functionality will not be automatic... but it will be easy to script separately if
	# necessary.
	# also note: aggregate processing (aka, processing a set of figures rather than a single figure) will be moved offline to 
	# paper/figure specific scripts.
	#
	# deprication note: this script will probably be depricated and replaced with something a bit simpler. this script
	# fits all the individual section_ids and compiles a collective catalog as it goes. 1) this complexity is not necessary,
	# and 2) we don't really need the "faultwise mean" analysis we've been doing (i don't think). the newer version will
	# take a more generalized list/list of dicts (or something) input so that catalogs can be better customized, like
	# section_ids = [1,2,3, [1,2], 4,5,6, [1,3,5,6]], etc.
	#
	# (waiting time probabilities, not 'Expected Waiting Time')
	# waiting time (or equivalently, "hazard function"?, plots. aka, probability of an earthquake given that one
	# has not occured for t0. equivalently, probability of interval Delta t > t0.
	#
	# prams:
	# (most are pretty self explanatory)
	#
	# for each section_id, fetch the mean_recurrence data.
	# 1) plot each cumulative probability
	# 2) including a weibull fit.
	# 3) do this for a set of waiting times:0, <t>/2, <t>, 1.5<t>, 2<t>, 2.5<t>
	# 4) a cumulative figure
	'''
	#
	if section_ids in (None, [], ()):
		section_ids='emc'
		#
	if isinstance(section_ids, str) and section_ids.lower()=='emc':
		section_ids = list(emc_section_filter['filter'])
	#
	section_ids = list(set(section_ids))
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
	fit_columns = ['t0', 'section_id', 'chi', 'beta', 'sigma_chi', 'sigma_beta', 'chi_sqr', 'fit_type']	# fit variables; map to dict.
	fit_type_str_len = 32
	fit_columns_types = ['float', 'int', 'float', 'float', 'float', 'float', 'float', 'S%d' % fit_type_str_len]
	for j, sec_id in enumerate(section_ids):
		#this_color = colors_[j%len(colors_)]
		#i=j+1
		i=j
		sections[sec_id] = {'fig':i}
		plt.figure(i, figsize=(12,10))
		plt.clf()
		#
		these_t0_factors = t0_factors		# default t0 factors, but we might change them for composite figures...
		#
		if sec_id == -1:
			# note: the full catalog gets compiled from the various catalog_id files (which are pre-calculated)
			#        (see the non-aggregate bits below).
			#
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
			these_t0_factors = [0., 50., 100., 150., 200., 250.]		# these numbers are fine for large catalogs (CA, so-CA, etc.,
																		# but might be problematic for some figures.
		#
		elif sec_id == -2:
			# "faultwise" composite:
			X = dT_composite_faultwise
			lw=5.
			ms=10.
			this_lbl = '(faultwise mean)'
			this_lbl_composite = 'all faults (faultwise)'
			sec_name = 'faultwise mean'
			these_t0_factors = [0., 50., 100., 150., 200., 250.]		# these numbers are fine for large catalogs (CA, so-CA, etc.,
																		# but might be problematic for some figures.
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
		best_fit_dict[sec_id]=[]
		#
		# do a preliminary fit for the whole set (should be equivalent to t0=0.0
		print "preliminary fit for section_id=%d" % sec_id
		fit_prams_0, fit_cov_0 = spo.curve_fit(f_weibull, xdata=numpy.array(X), ydata=numpy.array([j/float(len(X)) for j in numpy.arange(1., len(X)+1)]), p0=numpy.array([mean_dT, 1.5]))
		print fit_prams_0
		this_chi_0 = fit_prams_0[0]
		this_beta_0 = fit_prams_0[1]
		#
		print "--------------"
		for i_t, t0 in enumerate(this_t0s):
			this_color = colors_[i_t%len(colors_)]
			max_x = 0.
			min_x = 0.
			#this_X = [x-t0 for x in X if (x-t0)>=0.]
			this_X = [x for x in X if (x-t0)>=0.]
			if len(this_X)<5: continue		# ... because it won't fit...
			this_X.sort()
			#
			# skip it if there aren't very many data:
			#if len(this_X)<5: continue
			#
			N = float(len(this_X))
			Y = [float(j)/N for j in range(1, int(N)+1)]
			max_x = max(max_x, max(X))
			min_x = min(min_x, min(X))
			#
			plt.plot([x for x in this_X], Y, '.-', color = this_color, label='data, $t_0=%.3f$' % t0)
			# curve_fit() tends to break -- aka, not converge. in that event, do an MC fit
			print "these lens: ", len(this_X), len(Y)
			try:
				# f_weibull(x=None, chi=1.0, beta=1.0, x0=None)
				fit_type = 'spo.curve_fit'	# try a converging fit; if it fails, we'll use a fit_type='MC' (monte carlo)
				fit_prams, fit_cov = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=t0), xdata=numpy.array(this_X), ydata=numpy.array(Y), p0=numpy.array([max(1.0, mean_dT-t0), 1.5]))
				#			
				print "fitted (converging)...", fit_cov
				pram_sigmas = numpy.sqrt(numpy.diag(fit_cov))		# the uncertainties of the fit parameters are the diagonals of the covariance matrix... right?
				mean_chi_sqr = numpy.mean([(f_weibull(x=this_X[k], chi=fit_prams[0], beta=fit_prams[1], x0=t0)-Y[k])**2. for k, xx in enumerate(this_X)]) # in xrange(len(X))])
				stdErr = mean_chi_sqr/math.sqrt(N)
				print fit_cov
				print "fit_prams(%d/%d)[t0 = %f]: %s / %s" % (i, sec_id, t0, str(fit_prams), str(fit_prams_0))
				#
				#best_fit_dict[sec_id][t0] = [sec_id, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr]
				fit_vals = [t0, sec_id, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr, fit_type]
			except:
				fit_type = 'MC_%d' % mc_nits
				# converging fit failed. do an MC method:
				print "converging fit failed. try an MC approach:"
				prams_dict = {'chi':[0., 2.5*mean_dT], 'beta':[0.,6.], 'x0':[t0,t0]}	# can we automatedly guess at these prams?
																					# maybe at some point, we should devise a quasi-
																					# converging MC method.
				Zmin, Z=mcfitter(func=f_weibull, X=this_X, Y=Y, prams_dict=prams_dict, nits=mc_nits, n_cpus=n_cpus)
				Zmin.sort(order='chi_sqr')
				best_fit = Zmin[0]
				#
				fit_prams = [best_fit['chi'], best_fit['beta']]
				sigma_chi  = (prams_dict['chi'][1]-prams_dict['chi'][0])/mc_nits
				sigma_beta = (prams_dict['beta'][1]-prams_dict['beta'][0])/mc_nits
				#
				#best_fit_dict[sec_id][t0] = [sec_id, best_fit['chi'], best_fit['beta'], sigma_chi, sigma_beta, best_fit['chi_sqr']]
				fit_vals = [t0, sec_id, best_fit['chi'], best_fit['beta'], sigma_chi, sigma_beta, best_fit['chi_sqr'], fit_type]
			#
			fit_data = {col:val for col, val in zip(*[fit_columns, fit_vals])}
			#best_fit_dict[sec_id] += {col:val for col, val in zip(*[fit_columns, fit_vals])}
			best_fit_dict[sec_id] += [fit_vals]	# and we'll make a structured array at the end of it all.
			#
			X_fit = numpy.arange(min(this_X), max(this_X)*1.5, (max(this_X)-min(this_X))/500.)
			#
			# plot local t0 fit:
			plt.plot([x for x in X_fit], [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=t0) for x in X_fit], '--', color=this_color, lw=lw, ms=ms, label='$t_0=%.2f, \\beta=%.3f, \\tau=%.3f$' % (t0, fit_prams[1], fit_prams[0]))
			#
			# plot using t0=0 fit:
			plt.plot([x for x in X_fit], [f_weibull(x=x, chi=fit_prams_0[0], beta=fit_prams_0[1], x0=t0) for x in X_fit], '-.', color=this_color, lw=lw, ms=ms, label=None)
			#
			#
		print best_fit_dict[sec_id]
		#best_fit_dict[sec_id] = numpy.core.records.fromarrays(zip(*best_fit_dict[sec_id]), names=fit_columns, formats = [(type(x).__name__ for x in best_fit_dict[sec_id][0]])
		best_fit_dict[sec_id] = numpy.core.records.fromarrays(zip(*best_fit_dict[sec_id]), names=fit_columns, formats = fit_columns_types)
		#
		# save best fits:
		best_fit_dict[sec_id].dump('%s/VC_CDF_WT_fits_m%s_section_%d.npy' % (output_dir, str(m0).replace('.', ''), sec_id))
		#
		plt.legend(loc=0, numpoints=1)
		plt.gca().set_ylim([0., 1.1])
		plt.title('CDF for m>7 on fault section %s' % sec_name)
		plt.legend(loc=0, numpoints=1)
		plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t$ (years)' % m0)
		plt.ylabel('Probability $P(t)$')
		plt.savefig('%s/VC_CDF_WT_m%s_section_%d.png' % (output_dir, str(m0).replace('.', ''), sec_id))
		
		if keep_figs==False and not sec_id<0: plt.close(i)	# ?? sec_id in [list-o-sec_ids]...
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
#def waiting_time_figs_2(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir='VC_CDF_WT_figs', mc_nits=100000, n_cpus=None, start_year=10000, end_year=None):
def conditional_RI_figs(section_ids=[], file_path_pattern='data/VC_CFF_timeseries_section_%d.npy', m0=7.0, t0_factors = [0., .5, 1.0, 1.5, 2.0, 2.5], output_dir='VC_CDF_WT_figs', mc_nits=100000, n_cpus=None, start_year=10000, end_year=None, keep_figs=False, **kwargs):
	'''
	# ... note as soon as we work out of these early yoder et al. 2015ab papers, we can dump def waiting_time_figs() entirely... (see above)
	# ... and phase 2 of this retro-fit will be to separate this into a loop (which we'll move elsewhere) and a caller for a single cRI
	#  figure. groups of figurues will be handled by a callign script.
	#
	# and note that we screwed up the notation/nomenclature. these are "Conditional Recurrence Probabilities", P(t,t_0),
	# where waiting times are \Delta t = t_0 + t .
	# newer version of waiting_time_figs. section_ids [] is generalized to take integers or lists. the automatic aggregate
	# bits are removed. a simple way to get an aggregate of all elements is like:
	# my_list = [1,2,3,4]
	# input_list = my_list + [mylist] --> [1,2,3,4, [1,2,3,4]], which will produce figs for each element and the combined catalog [1,2,3,4].
	#
	# (waiting time probabilities, not 'Expected Waiting Time')
	# waiting time (or equivalently, "hazard function"?, plots. aka, probability of an earthquake given that one
	# has not occured for t0. equivalently, probability of interval Delta t > t0.
	#
	# prams:
	# (most are pretty self explanatory)
	, start_year=start_year=10000, end_year=end_year=None
	# start_year, end_year: start/end years of simulation. it's typically a good idea to sluff off the leading 10,000 years or so
	# (i assume this depends on the simulation size) if you've got a short catalog, 0 will be fine, but there tend to be
	# artifacts in the early parts of the 
	# fit_file_pattern: use an existing fit (for re-making figures), or calculate new fits (which can be time consuming). to always re-fit,
	# se fit_file_pattern=None
	#
	# for each section_id, fetch the mean_recurrence data.
	# 1) plot each cumulative probability
	# 2) including a weibull fit.
	# 3) do this for a set of waiting times:0, <t>/2, <t>, 1.5<t>, 2<t>, 2.5<t>
	# 4) a cumulative figure
	# TODO: ... and we should save all of these data (including the fits), so we can re-run these much faster.
	'''
	#
	# fetch some possible kwargs:
	# font sizes:
	fs_title  = kwargs.get('fs_title', 12)
	fs_legend = kwargs.get('fs_legend', 12)
	fs_label  = kwargs.get('fs_label', 12)
	#
	if section_ids in (None, [], ()):
		section_ids='emc'
		#
	if isinstance(section_ids, str) and section_ids.lower()=='emc':
		section_ids = list(emc_section_filter['filter'])
	#
	# set() breaks when trying to hash list objects, so do this manually:
	#section_ids = list(set(section_ids))
	for i, x in enumerate(section_ids):
		while section_ids.count(x)>1: section_ids.pop(i)
	#
	plt.ion()
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
	lw=2.5
	ms=5.
	max_x = 0.
	min_x = 0.
	full_catalog = []
	# control color cycling:
	colors_ =  mpl.rcParams['axes.color_cycle']
	sections ={'all_cdf':{'fig':0}}
	#best_fit_array = []		# ...and we'll cast this as a recarray later.
	best_fit_dict = {}
	fit_columns = ['t0', 'section_id', 'chi', 'beta', 'sigma_chi', 'sigma_beta', 'chi_sqr', 'fit_type']	# fit variables; map to dict.
	fit_type_str_len = 32
	#fit_columns_types = ['float', 'int', 'float', 'float', 'float', 'float', 'float', 'S%d' % fit_type_str_len]
	fit_columns_types = ['float', 'S256', 'float', 'float', 'float', 'float', 'float', 'S%d' % fit_type_str_len]
	#	
	plt.figure(i, figsize=(12,10))
	plt.clf()
	#
	for j, sec_id in enumerate(section_ids):
		#this_color = colors_[j%len(colors_)]
		#i=j+1
		i=j
		#sections[sec_id] = {'fig':i}
		if keep_figs: plt.figure(i, figsize=(12,10))
		plt.clf()
		#
		#these_t0_factors = t0_factors		# default t0 factors, but we might change them for composite figures...
		#
		# handle some section_id input types:
		# (but the better approach here might be to always pass these as tuples, not lists). tuples/integers can be
		# passed to dict objects as indices
		#if not hasattr(sec_id, 'append'):
		#	if isinstance(sec_id, tuple):
		#		sec_id=list(sec_id)		# can we always make these tuples?
		#	else:
		#		sec_id=[sec_id]
		#	print "wrapping sec_id: ", sec_id
		if not hasattr(sec_id, '__len__'): sec_id=tuple([sec_id])		# ignoring, for now, the possibility of strings, etc.
		if hasattr(sec_id, '__len__') and hasattr(sec_id, 'append'): sec_id=tuple(sec_id)
		#
		# we'll need an index for dict objects. they can take integers, etc. or tuples but not lists...
		# later, look into using integer or tuple types. if our sec_id lists --> tuples, we can pass sec_id
		#sec_id_index = sec_id if hasattr(sec_id, 'append') else tuple(sec_id)
		#
		# this bit about combinining catalogs needs to be reviewed and possibly re-written. how, specifically, do we want to 
		# facilitate custom catalogs?
		ary_in = combine_section_CFFs(sec_id, start_year=start_year, end_year=end_year)
		#ary_in = numpy.load(file_path_pattern % sec_id)
		#
		if full_catalog in (None, []) or len(full_catalog)==0:
			full_catalog = ary_in
		else:
			full_catalog = numpy.append(full_catalog, ary_in)
		#
		# get some mean interval stats and dN and dT (number of events and time between m>m0 events):
		mean_rec_data = mean_recurrence(ary_in=ary_in, m0=m0)
		X = mean_rec_data['dists']['dT'].tolist()			# time intervals; dNs are ['dists']['dN']
		#
		print "sec_id: ", sec_id
		this_lbl = 'section(s) %s' % ', '.join(map(str, sec_id))
		#sec_name = this_lbl
		#
		N = float(len(X))	
		X.sort()
		#
		# keep track of min/max X values for later global plot... or we could use them here as well.
		max_x = max(max_x, max(X))
		min_x = min(min_x, min(X))
		#
		# get t_0 values as factors of <t>.
		mean_dT = numpy.mean(X)
		this_t0s = [mean_dT * x for x in t0_factors]
		#
		best_fit_dict[sec_id]=[]
		#
		# do a preliminary fit for the whole set (should be equivalent to t0=0.0
		print "preliminary fit for section_id(s)=%s" % ', '.join([str(x) for x in sec_id])
		fit_prams_0, fit_cov_0 = spo.curve_fit(f_weibull, xdata=numpy.array(X), ydata=numpy.array([j/float(len(X)) for j in numpy.arange(1., len(X)+1)]), p0=numpy.array([mean_dT, 1.5]))
		print fit_prams_0
		this_chi_0 = fit_prams_0[0]
		this_beta_0 = fit_prams_0[1]
		#
		print "--------------"
		for i_t, t0 in enumerate(this_t0s):
			this_color = colors_[i_t%len(colors_)]
			max_x = 0.
			min_x = 0.
			#
			#this_X = [x for x in X if (x-t0)>=0.]		# why didn't we just say if x>t0?
			this_X = [x for x in X if x >= t0 ]
			#print "tyring... %s/%f because it's too short (%d)" % (str(sec_id), t0, len(this_X))
			if len(this_X)<=5: 
				print "skipping %s/%f because it's too short (%d)" % (str(sec_id), t0, len(this_X))
				continue					# ... because we won't be able to fit it...
			this_X.sort()					# ... though it should already be sorted...
			#
			N = float(len(this_X))
			Y = [float(j)/N for j in range(1, int(N)+1)]
			max_x = max(max_x, max(X))
			min_x = min(min_x, min(X))
			#
			#plt.plot([x for x in this_X], Y, '.-', color = this_color, label='data, $t_0=%.3f$' % t0)
			plt.plot(this_X, Y, '.-', color = this_color, label='data, $t_0=%.3f$' % t0)
			# curve_fit() tends to break -- aka, not converge. in that event, do an MC fit
			print "these lens: ", len(this_X), len(Y)
			try:
				# f_weibull(x=None, chi=1.0, beta=1.0, x0=None)
				fit_type = 'spo.curve_fit'	# try a converging fit; if it fails, we'll use a fit_type='MC' (monte carlo)
				fit_prams, fit_cov = spo.curve_fit(lambda x, chi, beta: f_weibull(x=x, chi=chi, beta=beta, x0=t0), xdata=numpy.array(this_X), ydata=numpy.array(Y), p0=numpy.array([max(1.0, mean_dT-t0), 1.5]))
				#			
				print "fitted (converging)...", fit_cov
				pram_sigmas = numpy.sqrt(numpy.diag(fit_cov))		# the uncertainties of the fit parameters are the diagonals of the covariance matrix... right?
				mean_chi_sqr = numpy.mean([(f_weibull(x=this_X[k], chi=fit_prams[0], beta=fit_prams[1], x0=t0)-Y[k])**2. for k, xx in enumerate(this_X)]) # in xrange(len(X))])
				stdErr = mean_chi_sqr/math.sqrt(N)
				print fit_cov
				print "fit_prams(%d/%d)[t0 = %f]: %s / %s" % (i, sec_id, t0, str(fit_prams), str(fit_prams_0))
				#
				#
				#fit_vals = [t0, sec_id, fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr, fit_type]
				fit_vals = [t0, str(sec_id), fit_prams[0], fit_prams[1], pram_sigmas[0], pram_sigmas[1], mean_chi_sqr, fit_type]
			except:
				fit_type = 'MC_%d' % mc_nits
				# converging fit failed. do an MC method:
				print "converging fit failed. try an MC approach:"
				prams_dict = {'chi':[0., 2.5*mean_dT], 'beta':[0.,6.], 'x0':[t0,t0]}	# can we automatedly guess at these prams?
																					# maybe at some point, we should devise a quasi-
																					# converging MC method.
				Zmin, Z=mcfitter(func=f_weibull, X=this_X, Y=Y, prams_dict=prams_dict, nits=mc_nits, n_cpus=n_cpus)
				Zmin.sort(order='chi_sqr')
				best_fit = Zmin[0]
				#
				fit_prams = [best_fit['chi'], best_fit['beta']]
				sigma_chi  = (prams_dict['chi'][1]-prams_dict['chi'][0])/mc_nits
				sigma_beta = (prams_dict['beta'][1]-prams_dict['beta'][0])/mc_nits
				#
				#				
				#fit_vals = [t0, sec_id, best_fit['chi'], best_fit['beta'], sigma_chi, sigma_beta, best_fit['chi_sqr'], fit_type]
				fit_vals = [t0, str(sec_id), best_fit['chi'], best_fit['beta'], sigma_chi, sigma_beta, best_fit['chi_sqr'], fit_type]
			#
			fit_data = {col:val for col, val in zip(*[fit_columns, fit_vals])}
			#best_fit_dict[sec_id] += {col:val for col, val in zip(*[fit_columns, fit_vals])}
			best_fit_dict[sec_id] += [fit_vals]	# and we'll make a structured array at the end of it all.
			#
			X_fit = numpy.arange(min(this_X), max(this_X)*1.5, (max(this_X)-min(this_X))/500.)
			#
			# plot local t0 fit:
			plt.plot([x for x in X_fit], [f_weibull(x=x, chi=fit_prams[0], beta=fit_prams[1], x0=t0) for x in X_fit], '--', color=this_color, lw=lw, ms=ms, label='$t_0=%.2f, \\beta=%.3f, \\tau=%.3f$' % (t0, fit_prams[1], fit_prams[0]))
			#
			# plot using t0=0 fit:
			plt.plot([x for x in X_fit], [f_weibull(x=x, chi=fit_prams_0[0], beta=fit_prams_0[1], x0=t0) for x in X_fit], '-.', color=this_color, lw=lw, ms=ms, label=None)
			#
			#
		print "best fit dict (final):", best_fit_dict[sec_id]
		#best_fit_dict[sec_id] = numpy.core.records.fromarrays(zip(*best_fit_dict[sec_id]), names=fit_columns, formats = [(type(x).__name__ for x in best_fit_dict[sec_id][0]])
		if len(best_fit_dict[sec_id])>0:
			best_fit_dict[sec_id] = numpy.core.records.fromarrays(zip(*best_fit_dict[sec_id]), names=fit_columns, formats = fit_columns_types)
			best_fit_dict[sec_id].dump('%s/VC_CDF_WT_fits_m%s_section_%s.npy' % (output_dir, str(m0).replace('.', ''), '_'.join([str(x) for x in sec_id])))
		#
		plt.legend(loc=0, numpoints=1, prop={'size':fs_legend})
		plt.gca().set_ylim([0., 1.1])
		plt.title('CDF for m>%s on fault section %s' % (str(m0), ', '.join([str(x) for x in sec_id])), size=fs_title)
		plt.legend(loc=0, numpoints=1, prop={'size':fs_legend})
		plt.xlabel('$m=%.2f$ Recurrence interval $\\Delta t_r$ (years)' % m0, size=fs_label)
		plt.ylabel('Probability $P(\\Delta t_r)$', size=fs_label)
		plt.savefig('%s/RI_conditional_CDF_m%s_section_%s.png' % (output_dir, str(m0).replace('.', ''), '_'.join([str(x) for x in sec_id])))
		
		#if keep_figs==False and not sec_id<0: plt.close(i)	# ?? sec_id in [list-o-sec_ids]...
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
##############################################################
# VQ data mining, etc. type helper functions and scripts:
#
# map related stuff (that might get moved to vc_geodetic? maybe. this will all be reorganized in PyVQ
#

def get_fault_model_extents(section_ids=None, sim_file=allcal_full_mks, n_cpus=None):
	# note: this, presently, is fully guessing at the (x,y) <--> (lon,lat) conversion.
	# note also: this could probably be sped up considerably by shifting the min(), max() to the MPP process, the basic strategy
	# to reduce the amount of data pickled back to the parent process.
	if n_cpus==None: n_cpus = mpp.cpu_count()
	#
	if section_ids in ('emc', 'EMC'): section_ids = emc_section_filter['filter']
	#section_ids = (section_ids or emc_section_filter['filter'])
	#print "section_ids: ", section_ids
	#
	with h5py.File(sim_file) as vc_data:
		section_ids = (section_ids or set(vc_data['block_info_table']['section_id'].tolist()))
		if n_cpus==1:
			block_info = numpy.array([rw for rw in vc_data['block_info_table'] if rw['section_id'] in section_ids], dtype=vc_data['block_info_table'].dtype)
		else:
			pool = mpp.Pool(n_cpus)
			results = []
			tbl = vc_data['block_info_table']
			N_len = len(tbl)
			dN = int(N_len/n_cpus)
			#
			for i in xrange(n_cpus):
				#
				N_term = (1+i)*dN
				if i==N_len-1: N_term = N_len
				results+=[pool.apply_async(in_rec_array, args=(), kwds={'rec_array':tbl[(i)*dN:N_term], 'col_name':'section_id', 'in_list':section_ids})]
				#
			pool.close()
			pool.join()
			#for res in results: print "type: ", type(res)
			#
			if len(results)>1:
				#
				#block_info = numpy.array(reduce(numpy.append, [numpy.array(x.get()) for x in results]), dtype=tbl.dtype)
				block_info = numpy.array([], dtype=tbl.dtype) # numpy.array(results[0].get(), dtype=tbl.dtype)
				for res in results:
					try:
						x = res.get()
						if len(x)>0:
							#
							block_info = numpy.append(block_info, x)
							print "appending... ", len(x)
						else:
							print "len 0 return..."
					except:
						print "failing to append res: ", res, res.get()
				#a=[numpy.append(block_info, x.get()) for x in results[1:]]
			else:
				# reduce() will throw an error if you give it only one value.
				block_info = results[0].get()
				
		#	
		# names=output_names, formats = [type(x).__name__ for x in outputs[0]])
		#t0=time.time()
		#print "do reduce: %s", str(t0)
		min_x = min([min(x) for x in [block_info['m_x_pt1'],  block_info['m_x_pt2'],block_info['m_x_pt3'],  block_info['m_x_pt4']]])
		max_x = max([max(x) for x in [block_info['m_x_pt1'],  block_info['m_x_pt2'],block_info['m_x_pt3'],  block_info['m_x_pt4']]])
		#print "finished reduce X: ", time.time(), time.time()-t0
		
		#
		#t0=time.time()
		#print "do reduce Y: %s", str(t0)
		min_y = min([min(x) for x in [block_info['m_y_pt1'],  block_info['m_y_pt2'],block_info['m_y_pt3'],  block_info['m_y_pt4']]])
		max_y = max([max(x) for x in [block_info['m_y_pt1'],  block_info['m_y_pt2'],block_info['m_y_pt3'],  block_info['m_y_pt4']]])
		#print "finished reduce Y: ", time.time(), time.time()-t0
		#
	#
	min_lat, min_lon = xy_to_lat_lon(x=min_x, y=min_y, sim_file=sim_file, return_format='tuple')
	max_lat, max_lon = xy_to_lat_lon(x=max_x, y=max_y, sim_file=sim_file, return_format='tuple')
	#
	return {'lat_min': min_lat, 'lat_max':max_lat, 'lon_min':min_lon, 'lon_max':max_lon}
#
def vc_basemap(projection='cyl', resolution='i', **kwargs):
	#
	titlefont1 = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
	titlefont2 = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14, weight='bold')
	sectionkeyfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=7)
	ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
	framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
	legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
	smtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9, weight='bold')
	cbticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
	ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
	#
	water_color             = '#bed5ff'
	land_color              = '#ffffff'
	seq_land_color          = '#ffffff'
	boundary_color          = '#000000'
	coastline_color         = '#9a9a9a'
	country_color           = '#9a9a9a'
	state_color             = '#9a9a9a'
	fault_color             = '#000000'
	alt_fault_color         = '#737373'
	selected_fault_color    = '#FFFFFF'
	map_tick_color          = '#000000'
	map_frame_color         = '#000000'
	grid_color              = '#000000'
	cb_fontcolor            = '#000000'
	#
	boundary_width          = 2.0
	coastline_width         = 2.0
	country_width           = 2.0
	state_width             = 2.0
	fault_width             = 0.5
	forecast_fault_width    = 6.0
	seq_fault_width_max     = 6.0
	seq_fault_width_min     = 3.0
	map_frame_width         = 1.0
	grid_width              = 0.5
	num_grid_lines          = 5
	#
	sp_line_color           = '#000000'
	#sp_line_colormap        = sequence_cmap
	sp_line_width           = 2.0
	#
	t0_dt_main_line_color   = '#000000'
	t0_dt_sub_line_color    = '#737373'
	t0_dt_main_line_width   = 2.0
	t0_dt_sub_line_width    = 1.0
	#t0_dt_range_color       = sequence_cmap(0.99)
	#
	# catch some defaults:
	try:
		projection
	except:
		projection = 'cyl'
	llcrnrlon=kwargs['llcrnrlon']
	llcrnrlat=kwargs['llcrnrlat']
	urcrnrlon=kwargs['urcrnrlon']
	urcrnrlat=kwargs['urcrnrlat']
	#
	lat_0 = llcrnrlat + (urcrnrlat - llcrnrlat)/2.
	lon_0 = llcrnrlon + (urcrnrlon - llcrnrlon)/2.
	#
	bm = bmp.Basemap(projection=projection, resolution=resolution, **kwargs)
	aspect = bm.aspect
	#
	bm.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=water_color)
	
	# Fill the continents and lakes.
	bm.fillcontinents(color=land_color, lake_color=water_color)
	
	# draw coastlines, edge of map.
	bm.drawcoastlines(color=coastline_color, linewidth=coastline_width)

	# draw countries
	bm.drawcountries(linewidth=country_width, color=country_color)

	# draw states
	bm.drawstates(linewidth=state_width, color=state_color)

	# draw parallels.
	parallels = numpy.linspace(llcrnrlat, urcrnrlat, num_grid_lines+1)
	mm_parallels = bm.drawparallels(
		parallels,
		labels=[1,0,0,0],
		color=grid_color,
		fontproperties=ticklabelfont,
		fmt='%.2f',
		linewidth=grid_width,
		dashes=[1, 10]
	)
	bm.drawrivers(color='b')

	# draw meridians
	meridians = numpy.linspace(llcrnrlon, urcrnrlon, num_grid_lines+1)
	mm_meridians = bm.drawmeridians(
		meridians,
		labels=[0,0,1,0],
		color=grid_color,
		fontproperties=ticklabelfont,
		fmt='%.2f',
		linewidth=grid_width,
		dashes=[1, 10]
	)
	
	return bm
#
def get_vc_fault_polygon(fault_block_vertices):
	# a sort of sloppy script to make a polygon out of a set of block vertices (see get_fault_traces() ).
	# assume the blocks are an array with fixed z levels...
	# but in the end, i guess we still have to make vectors out of the blocks
	print "this is not finished yet..."
	#
	Z = zip(*fault_block_vertices)[2]
	#
	min_z = min(Z)
	max_z = max(Z)
	#
	#reduced_blocks = [block in fault_block_vertices if block[2] in (min_z, max_z)]
	#
	# ... and then finish it later...
#
def simple_fault_trace(fault_block_vertices=None, vert_cols=[2,3,4,5]):
	# falut_block_vertices from get_fault_traces output (just the vertices).
	# like: [section_id, block_id, [verts: xyz_UL, xyz_LL, xyz_LR, xyz_UR]] (or something like this)
	# each vert is a lists-like: [x,y,z]
	#
	vecs = []
	# just mash all the block vertices into a set of points, make a set() of them, and assume they're in
	# the correct sequence?
	# for now, just keep [x,y]
	#
	#traces = []	# list of traces. each entry will be [[x,y],[x,y]...] in order of tip-to-tail.
	#
	for rw in fault_block_vertices:
		#verts = rw[2:6]
		verts = [rw[i] for i in vert_cols]	# note, these can be named or numbered depending on the object type...
		#print verts
		#points+= [(x[0],x[1]) for x in verts if (x[0],x[1]) not in points]
		#
		vecs += [[verts[0][0:2], verts[-1][0:2]] ]
		if len(verts)>3: vecs += [[verts[1][0:2], verts[2][0:2]]]
	#
	return vecs
#
def simple_fault_trace2(fault_block_vertices=None, vert_cols=[2,3,4,5]):
	# slightly more compact (and will probably replace) simple_fault_trace()
	# falut_block_vertices from get_fault_traces output (just the vertices).
	# like: [section_id, block_id, [verts: xyz_UL, xyz_LL, xyz_LR, xyz_UR]] (or something like this)
	# each vert is a lists-like: [x,y,z]
	#
	print "this is not working yet..."
	return "this is not working yet..."
	#
	vecs = []
	# just mash all the block vertices into a set of points, make a set() of them, and assume they're in
	# the correct sequence?
	# for now, just keep [x,y]
	#
	#traces = []	# list of traces. each entry will be [[x,y],[x,y]...] in order of tip-to-tail.
	#
	for rw in fault_block_vertices:
		#verts = rw[2:6]
		verts = [rw[i] for i in vert_cols]	# note, these can be named or numbered depending on the object type...
		#print verts
		#points+= [(x[0],x[1]) for x in verts if (x[0],x[1]) not in points]
		#
		#vecs += [zip(*[verts[0][0:2], verts[-1][0:2]])]
		#if len(verts)>3: vecs += [zip(*[verts[1][0:2], verts[2][0:2]])]
		vecs += [ [verts[0][0:2], verts[-1][0:2]] ]
		if len(verts)>3: vecs += [ [verts[1][0:2], verts[2][0:2]] ]
	#return vecs
	#
	# now, paste these verts together into connected fault traces.
	# note that each vector is a mini-trace, so let's put them together:
	#
	traces = [vecs.pop()]
	#
	print "len(vecs): %d" % len(vecs)
	#
	# ... but i'm not sure this actually works, because fault blocks don't always (???) actually touch?
	while len(vecs)>0:
		for j, trace in enumerate(traces):
			found_one=False
			for i,vec in enumerate(vecs):
				#print vec, trace
				if vec[0]==trace[-1]:
					traces[j]+=vec
					found_one=True
				if vec[-1]==trace[0]:
					traces[j] = vec + trace
					found_one=True
				#
		
			if found_one==False and len(vecs)>0:
				print "adding new trace..."
				traces += [vecs.pop()]
	
	#
	print "lens: %d, %d" % (len(traces), len(vecs))
	return traces
#
def get_fault_blocks(section_ids=None, sim_file=allcal_full_mks):
	'''
	# this is a little bit mis-named. returns a dict. of faults, like:
	# {fault_id:[list-o-blocks]...}
	#
	# eventually, we'll want different types of fault traces: show all the blocks, show sections, show just a trace, maybe
	# top/bottom trace, top/bottom/edges polygon (aka, hollow it out). we want, ultimately, to return a set of [[X,Y, Z?] ]
	# pairs/triplets for plotting.
	#
	fault_col_names = ('section_id', 'block_id', 'UL', 'LL', 'LR', 'UR')
	'''
	faults = {}
	#
	# first, get all block_ids and separate into faults (sections?):
	with h5py.File(sim_file) as vc_data:
		block_data = vc_data['block_info_table']
		if section_ids == None: section_ids = set(block_data['section_id'].tolist())
		#
		print "get_fault_blocks() using section_ids=", section_ids
		block_data_iterator = block_data.__iter__()
		#for rw in block_data:
		i_max = len(block_data)-2		# aka, the [-2] position...
		i_blk=0
		for rw in block_data_iterator:
			if rw['section_id'] not in section_ids:
				continue
				#
				# this should be a fast way to skip extraneous rows, but the iterator does not seem to work properly...
				'''
				print "skipping section_id: %d" % rw['section_id']
				# do a fast-skip thorought these elements (without searching the list again).
				#continue
				bogus_section_id=rw['section_id']
				while rw['section_id']==bogus_section_id and i_blk<i_max:
					# manually increment j and rw:
					try:
						rw=block_data_iterator.next()
						i_blk+=1
					except:
						break
				bogus_section_id=None
				print "new section_id: %d" % rw['section_id']
				#
				'''	
			#
			# now we have a valid section:
			if faults.has_key(rw['fault_id'])==False:
				faults[rw['fault_id']] = []
				print "adding fault_id: %d/%d" % (rw['fault_id'], rw['section_id'])
			#if not hasattr(faults[rw['fault_id']], '__len__'): faults[rw['fault_id']]=[]
			#
			# and for now, let's not separate the sections; just note the section_id:
			UL, LL, UR = [[rw['m_%s_pt%d' % (xyz, j)] for xyz in ['x', 'y', 'z']] for j in [1,2,4]]	# note: this may need to be
																									# corrected for newer vc where there are 3, not 4 vertices.
			LR = numpy.array(LL) + numpy.array(UR)-numpy.array(UL)
			#
			faults[rw['fault_id']] += [[rw['section_id'], rw['block_id'], UL, LL, LR.tolist(), UR]]
			i_blk+=1
			#
		#for key in faults.iterkeys():
		#	# ... can't encode list-in-list, only simpler structures... for now, just return the list.
		#	print "key-rw: ", faults[key][0]
		#	faults[key] = numpy.core.records.fromarrays(zip(*faults[key]), names=fault_col_names, formats = [type(x).__name__ for x in faults[key][0]])
			
	#
	# returns a dict like {fault_id:[list-o-blocks]...}
	# where each list-o-blocks is like [section_id, block_id, vert0, vert1, vert2, vert3
	return faults
#
def get_nearest_section_ids(lat_0=emc_event['lat'], lon_0=emc_event['lon'], n_sections=5, section_ids=None, sim_file=default_sim_file, dist_mode=0, verbose=False, fignum=None):
	'''
	# return a list of the section_ids nearest to (lon, lat).
	# lat, lon: closest to this position
	# n_sections: number of sections to return
	# section_ids: an initial section_id filter, either for speed or some other purpose.
	# sim_file: the sim_file (data) to use
	# dist_mode: 0: section center (mean position), 1: closest any point/vertex.
	# ... but note that as it stands, this is the center/closest point of the blocks, not the section itsef (that will be harder
	# and we'll save it for later), so these modes are all but identical.
	#
	# note: there may be a sloppy handling of the 3 vs 4 verts data format that will (in some cases) not properly resolve depth.
	# we'll probably, in cases where there are only 3 verts, simply repeat the final vert which should give (relatively) accurate
	# x,y positions, but not depth.
	'''
	#
	# convert input lat/lon to vc x,y. this will change using the new VC, since the native components for new_vc will be lat, lon.
	# this will work and has a nice, general syntax:
	#xy = lat_lon_to_xy(lat, lon, sim_file=sim_file, lat0=None, lon0=None, chi=111.1, return_format='dict')
	#x,y = [xy[s] for s in ['x', 'y']
	#
	x0,y0 = lat_lon_to_xy(lat_0, lon_0, sim_file=sim_file, lat0=None, lon0=None, chi=111.1, return_format='tuple')
	#
	my_blocks = []	# we'll make it a recarray. note we're not going to return the block/section data, just the section_ids and their dists.
	#				# for mode=1 (any point), we'll return the nearest dist.
	#
	with h5py.File(sim_file) as vc_data:
		block_data = vc_data['block_info_table']
		#if section_ids == None: section_ids = set(block_data['section_id'].tolist())
		#XYZ = ['x', 'y', 'z']
		# do we have 3 or 4 data points?
		j_verts = range(1,4) + [3]
		if 'm_x_pt_4' in block_data.dtype.names: j_verts=range(1,5)		# ... or _z_, _y_... . indeed, we have a 4th point, use it.
		#
		for rw in block_data:
			if section_ids!=None and rw['section_id'] not in section_ids:
				continue
			#
			#X = [ [rw['m_%s_pt%d' % (xyz, j)] for xyz in ['x', 'y', 'z']]  for j in j_verts]		# rows: [pos1, pos2, pos3...]
			if dist_mode==0:
				X = [ [rw['m_%s_pt%d' % (xyz, j)] for j in j_verts]  for xyz in ['x', 'y', 'z']]		# columns: [X, Y, Z]
				mean_pos = [numpy.mean(x) for x in X]
				dist_sqr = (mean_pos[0]-x0)**2. + (mean_pos[1]-y0)**2.		# we're just rank-ordering, so save the time of sqrt().
				#dist = numpy.linalg.norm(mean_pos[0:2])
				#
				my_blocks += [[rw['block_id'], rw['section_id'], dist_sqr, mean_pos[0], mean_pos[1]]]
			else:
				# dist_mode==1 or anything else:
				X = [ [rw['m_%s_pt%d' % (xyz, j)] for xyz in ['x', 'y', 'z']]  for j in j_verts]		# ... rows of unique positions.
				my_blocks +=  [ [rw['block_id'], rw['section_id'], (x[0]-x0)**2. + (x[1]-y0)**2., x[0], x[1]] for x in X] 
			#
		#
	#
	# we might make a recarray, but since this is all internal (to this function), it's not really necessary...
	#my_blocks = numpy.core.records.fromarrays(zip(*my_blocks), names=['block_id', 'section_id', 'dist', 'lat', 'lon'], formats=[type(x).__name__ for x in my_blocks[0]])
	#my_blocks.sort(order='dist')
	my_blocks.sort(key=lambda x: x[2])
	#
	# now, we want the nearest n_sections.
	return_section_ids = {}
	for rw in my_blocks:
		# my_blocks is sorted by distance. now get the closest unique section_id.
		if not return_section_ids.has_key(rw[1]): return_section_ids[rw[1]]=rw
		# ... and just skip all the rest...
		if len(return_section_ids)>=n_sections: break
	#
	section_rows = [val for val in return_section_ids.itervalues()]
	#section_rows.sort(key = lambda x: x[2])	# (though it should already be sorted).	
	#
	# just the section_ids:
	return_section_ids = [rw[1] for rw in section_rows[0:n_sections]]
	#return_section_ids = [rw[1] for rw in unique_rows[0:n_sections]]
	#
	if fignum!=None:
		plt.figure(fignum)
		plt.clf()
		print "sections: ", return_section_ids
		#ft = get_fault_traces(section_ids = return_section_ids, fignum=fignum, lat_lon=True)
		ft = get_block_traces(section_ids=return_section_ids, fignum=fignum, lat_lon=True)
		plt.plot(lon_0, lat_0, 'r*', ms=18)
	#
	#return my_blocks
	verbose = False
	if verbose:
		return section_rows
	else:
		return return_section_ids
#
def get_block_traces(fault_blocks=None, section_ids=None, sim_file=allcal_full_mks, fignum=None, lat_lon=True ):
	# like get_fault_traces(), but don't bother with the faultwise partitioning...
	#
	# some error handling:
	if isinstance(section_ids, int): section_ids=[section_ids]
	#
	if fault_blocks==None:
		with h5py.File(sim_file) as vc_data:
			blks = vc_data['block_info_table']
			fault_blocks = [rw.tolist() for rw in blks if rw['section_id'] in section_ids]
			# numpy.core.records.fromarrays(zip(*field_data_prime), names=['x', 'y', 'z', 'dx', 'dy', 'dz', 'dxyz', 'dxy'], formats=[type(x).__name__ for x in field_data_prime[0]])
			fault_blocks = numpy.core.records.fromarrays(zip(*fault_blocks), dtype=blks.dtype)
	#
	traces = []
	#
	for rw in fault_blocks:
		UL, LL, UR = [[rw['m_%s_pt%d' % (xyz, j)] for xyz in ['x', 'y', 'z']] for j in [1,2,4]]	# note: this may need to be
																								# corrected for newer vc where there are
																								# 3, not 4 vertices.
		LR = numpy.array(LL) + numpy.array(UR)-numpy.array(UL)
		#
		if lat_lon==True: 
			UL = xy_to_lat_lon(UL[0], UL[1], sim_file=sim_file, return_format='lonlat_list') + [UL[2]]
			UR = xy_to_lat_lon(UR[0], UR[1], sim_file=sim_file, return_format='lonlat_list') + [UR[2]]
			LR = xy_to_lat_lon(LR[0], LR[1], sim_file=sim_file, return_format='lonlat_list') + [LR[2]]
			LL = xy_to_lat_lon(LL[0], LL[1], sim_file=sim_file, return_format='lonlat_list') + [LL[2]]
					#
		traces += [[UL, UR]]
		traces += [[LL, list(LR)]]
		#
	#
	if fignum!=None:
		plt.figure(fignum)
		plt.clf()	
		#
		for rw in traces:
			X,Y,Z = zip(*rw)
			plt.plot(X,Y, '-', color='r')
	#
	return traces		
#
def get_fault_traces(fault_blocks=None, section_ids=None, sim_file=allcal_full_mks, fignum=None, lat_lon=True ):
	# get some simple fault traces for plotting (or whatever):
	if isinstance(section_ids, int): section_ids=[section_ids]
	#
	# get faultwise collections of blocks (by section_id)
	if fault_blocks == None:
		fault_blocks = get_fault_blocks(section_ids=section_ids, sim_file=sim_file)
	#
	fault_traces = {}
	#
	for fault_id,fault in fault_blocks.iteritems():
		fault_traces[fault_id] = simple_fault_trace(fault_block_vertices=fault, vert_cols=[2,3,4,5])
	#
	# if lat_lon, then convert from xy to lat_lon:
	# use: xy_to_lat_lon(x, y, sim_file=allcal_full_mks, lat0=None, lon0=None, chi=111.1, return_format='dict')
	if lat_lon:
		for fault_id, trace in fault_traces.iteritems():
			for i, pair in enumerate(trace):
				# each entry in the trace is like [ [x1,x2], [y1, y2] ]
				#print "pair: ", pair
				fault_traces[fault_id][i] = [xy_to_lat_lon(rw[0], rw[1], sim_file=sim_file, return_format='lonlat_list') for rw in pair]
				
	#
	# ... and for now, plot it here. we'll move this out of funct. later...
	if fignum!=None:
		plt.figure(fignum)
		plt.clf()
		for fault_id, trace in fault_traces.iteritems():
			for pair in trace:
				#X,Y = zip(*pair)
				#plt.plot(X,Y, '-')
				plt.plot([x[0] for x in pair], [y[1] for y in pair], '-')
		#
	#
	return fault_traces
#
def plot_fault_traces(traces=None, section_ids=None, sim_file=allcal_full_mks, fignum=0, lat_lon=True, do_clf=False, plot_color=None, zorder=5, line_style = '-'):
	if traces==None:
		traces = get_block_traces(section_ids=section_ids, sim_file=sim_file, fignum=None, lat_lon=lat_lon)
	#
	plt.figure(fignum)
	if do_clf: plt.clf()	
	#
	if plot_color==None: plot_color='r'
	for i, rw in enumerate(traces):
		X,Y,Z = zip(*rw)
		lbl=None
		if i==0: lbl='sec_id=%s' % ', '.join(map(str, section_ids))
		plt.plot(X,Y, '%s%s' % (plot_color, line_style), label=lbl, zorder=zorder)
	plt.legend(loc=0, numpoints=1)
#
def get_anss_seismicity(section_ids=None, sim_file=allcal_full_mks, start_date=None, end_date=None, m_c=3.0, n_max=999999, n_cpus=None):
	# make a map of real seismicity around our model area. use the fault model to determine extents.
	#
	# emc section_ids: emc_section_filter['filter']
	#
	# handle some default values and book-keeping:
	# ...
	if end_date==None: end_date=dtm.datetime.now(pytz.timezone('UTC'))
	if start_date==None: start_date = end_date-dtm.timedelta(days=500)
	#
	#
	ll_range = get_fault_model_extents(section_ids=section_ids, sim_file=sim_file, n_cpus=n_cpus)
	#
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	#
	# ... but we'll use an ANSStools function.
	earthquake_catalog = ANSStools.catfromANSS(lon=[ll_range['lon_min'], ll_range['lon_max']], lat=[ll_range['lat_min'], ll_range['lat_max']], minMag=m_c, dates0=[start_date, end_date], Nmax=n_max, fout=None)
	#
	# now, let's cast the catalog as a recarray (which will soon be done in ANSStools):
	if isinstance(earthquake_catalog, numpy.recarray)==False:
		# might need to use: numpy.core.records.recarray
		# (can also use "names=[]", "formats=[]" syntax, but note that type(datetime).__name__ does
		# not produce a viable name-type for numpy.rec.array().
		# pass
		earthquake_catalog = numpy.rec.array(earthquake_catalog, dtype=[('event_date', 'M8[us]'), ('lat','f'), ('lon','f'), ('mag','f'), ('depth','f')])
	#	
	return earthquake_catalog	
#
def seismicity_map(section_ids=None, sim_file=allcal_full_mks, start_date=None, end_date=None, n_cpus=None, fignum=0, map_size=[10,8], etas_gridsize=.25, etas_mc=3.0, etas_contour_intervals=24, etas_cat_len=500, etas_catalog=None, p_map=0.0, map_resolution = 'i'):
	# make a map of real seismicity around our model area. use the fault model to determine extents.
	#
	# emc section_ids: emc_section_filter['filter']
	#
	# handle some default values and book-keeping:
	#
	# allow knowing users to pass a sweep index with the etas_catalog, so we can speed up the vc_catalog fetch.
	sweep_index=None
	if isinstance(etas_catalog, dict):
		sweep_index = etas_catalog['sweep_index']
		etas_catalog = etas_catalog['catalog']
		
	# ...
	if end_date==None: 
		if isinstance(start_date, float) or isinstance(start_date, int):
			with h5py.File(sim_file) as vc_data:
				end_date = 365.25*(max(vc_data['event_table']['event_year']))
			#
		#
		else:
			end_date=dtm.datetime.now(pytz.timezone('UTC'))
	if start_date==None:
		start_date = end_date-dtm.timedelta(days=500)
	#
	#
	ll_range = get_fault_model_extents(section_ids=section_ids, sim_file=sim_file, n_cpus=n_cpus)
	#
	lon_0 = ll_range['lon_min'] + (ll_range['lon_max']-ll_range['lon_min'])/2.
	lat_0 = ll_range['lat_min'] + (ll_range['lat_max']-ll_range['lat_min'])/2.
	#
	plt.figure(fignum, figsize=map_size)
	plt.clf()
	#bm = vc_basemap(llcrnrlon=ll_range['lon_min'], llcrnrlat=ll_range['lat_min'], urcrnrlon=ll_range['lon_max'], urcrnrlat=ll_range['lat_max'], lon_0=lon_0, lat_0=lat_0, resolution='i', projection='cyl')
	bm = vc_basemap( projection='cyl', llcrnrlon=ll_range['lon_min'], llcrnrlat=ll_range['lat_min'], urcrnrlon=ll_range['lon_max'], urcrnrlat=ll_range['lat_max'], lon_0=lon_0, lat_0=lat_0, resolution=map_resolution)
	
	# note: we could also specify lon_0, lat_0, width, height {in meters}.
	#bm.drawcoastlines()
	#bm.drawmapboundary()
	#bm.drawlakes()
	#bm.drawrivers()
	plt.title('VC fault model map\n\n')
	plt.show()
	#
	# looks like there might be some GIT synching problems to be handled here. maybe a change from Buller didn't manage
	# to push up properly?
	print "start_date, end_date: ", start_date, end_date
	#
	if etas_catalog in (None, 'etas'):
		etas_catalog = BASScast.getMFETAScatFromANSS(lons=[ll_range['lon_min'], ll_range['lon_max']], lats=[ll_range['lat_min'], ll_range['lat_max']], dates=[dtm.datetime.now(pytz.timezone('UTC'))-dtm.timedelta(days=etas_cat_len), dtm.datetime.now(pytz.timezone('UTC'))], mc=etas_mc)
	elif etas_catalog == 'vc':
		# make a vc type catalog:
		etas_catalog = vc_ETAS_catalog(section_ids=section_ids, sim_file=sim_file, start_year=start_date, end_year=end_date, n_cpus=n_cpus, fignum=fignum, map_size=map_size, etas_mc=etas_mc, sweep_index=sweep_index)	# and a bunch of those prams can actually be removed.
	else:
		# use the provided catalog...
		pass
	#
	print "catalog selected. calculate etas."
	#
	# instantiate a BASScast object() (from the BASScast module).
	# basscast needs a list object, otherwise its sorting routines get screwed up (we might fix this on the basscast side as well...)
	#this_etas_catalog = etas_catalog
	#if hasattr(etas_catalog, 'tolist'): this_etas_catalog=etas_catalog.tolist()
	#
	#return etas_catalog
	etas = BASScast.BASScast(incat=etas_catalog, fcdate=end_date, gridsize=etas_gridsize, contres=etas_contour_intervals, mc=etas_mc, eqeps=None, eqtheta=None, fitfactor=5., contour_intervals=etas_contour_intervals, lons=[ll_range['lon_min'], ll_range['lon_max']], lats=[ll_range['lat_min'], ll_range['lat_max']], rtype='ssim', p_quakes=1.05, p_map=p_map)
	#
	#conts = etas.getContourSet(X_i=None, Y_i=None, Z_ij=None, contres=etas_contour_intervals, zorder=7, alpha=.15)
	#conts2 = etas.BASScastContourMap(fignum=3, maxNquakes=10, alpha=.75)
	#
	# and we can plot these contours anywhere like:
	# X,Y,Z = a.X_i, a.Y_i, a.Z2d
	# plt.contourf(X,Y,Z, alpha=.25) (noting that as X,Y will be in lat/lon coords.)
	#
	# now, draw some fault traces:
	# fault_blocks=None --> will fetch fault blocks from get_fault_blocks()
	fault_traces = get_fault_traces(fault_blocks=None, section_ids=section_ids, sim_file=sim_file, fignum=None, lat_lon=True)
	# format is like:
	#plt.figure(fignum)
	#plt.clf()
	#for fault_id, trace in fault_traces.iteritems():
	#	for pair in trace:
	#		#X,Y = zip(*pair)
	#		#plt.plot(X,Y, '-')
	#		plt.plot([x[0] for x in pair], [y[1] for y in pair], '-')
	# 
	# convert to map coordis using etas.cm()
	plt.figure(fignum)
	colors_ =  mpl.rcParams['axes.color_cycle']
	fault_color = 'b'
	for fault_index, (fault_id, trace) in enumerate(fault_traces.iteritems()):
		for pair in trace:
			#X,Y = zip(*pair)
			#plt.plot(X,Y, '-')
			#fault_color = colors_[fault_index%len(colors_)]
			#
			new_pair = [etas.cm(x,y) for x,y in pair]
			plt.plot([x[0] for x in new_pair], [y[1] for y in new_pair], '-', color=fault_color, lw=1.5)
		#	
	#
	return etas
#
def vc_ETAS_catalog(section_ids=None, sim_file=allcal_full_mks, start_year=None, end_year=None, n_cpus=None, fignum=0, map_size=[10,8], etas_mc=3.0, sweeps_index=None):
	'''
	# make an ETAS map from VC catalog events.
	# note: p_map gives the mapping representation of the Omori scaling exponent p. p_map -> 0.0 gives the time independend
	# RI type, historical rates, as opposed to the local time dependent rates.
	# we'll also need to work out bits about plotting large events, etc.
	#
	# ... and then we'll feed this into the existing BASS map thing above...
	'''
	#
	# first, get a catalog:
	vc_catalog = combine_section_CFFs(section_ids, start_year=start_year, end_year=end_year)
	# now, convert this catalog to an ETAS-like format. our present catalog looks like:
	# ('event_number',  'event_year',  'event_magnitude', 'cff_initial', 'cff_final', 'event_area')
	# sampel BASScast catalog call:
	# cat1=BASScast.getMFETAScatFromANSS(lons=[-117., -114.], lats=[32., 35.], dates=[dtm.datetime(1990,1,1, tzinfo=BASScast.pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], mc=3.0)
	# rows (basically from ANSS) are like:
	# [726468.958140625, 32.561, -115.805, 3.31, 4.19] --> [year, lat, lon, mag, depth]	(eventually, these too will be returned as recarrays).
	# so to do ETAS, we have to find the event's epicenter... and this is going to be messy and slow... but let's just code it:
	# (note this can also be done by just treating the individual blocks as little earthquakes, but it will give a somewhat different result.
	# also, we may still need to reconstruct the catalog, though this might be basically done in vc_geodetic. we could take the time series of
	# block displacements and treat each one like a little earthquake).
	vc_etas_catalog=[]
	with h5py.File(sim_file) as vc_data:
		sweeps = vc_data['event_sweep_table']
		blocks = vc_data['block_info_table']
		#
		# make a sweeps index (this may alredy exist in some form, but this will be quick). 
		# we want the sweeps index for the first sweep of each event.
		if sweeps_index==None: sweeps_index = {int(rw['event_number']):int(rw['block_id']) for rw in sweeps if rw['sweep_num']==0}	# ... and this takes a little while.
		#
		for rw in vc_catalog:
			if rw['event_magnitude']<etas_mc: continue
			# get initial block:
			blk = blocks[sweeps_index[rw['event_number']]]
			#print "blk: ", blk
			#print 
			x = numpy.mean([blk['m_x_pt%d' % k] for k in range(1,5)])	# these will be a little different with new VC since there are only 3 verts.
			y = numpy.mean([blk['m_y_pt%d' % k] for k in range(1,5)])
			z = - numpy.mean([blk['m_z_pt%d' % k] for k in range(1,5)])
			#
			z/=1000.		# in km...
			#
			# now lat/lon conversion:
			# def xy_to_lat_lon(x, y, sim_file=allcal_full_mks, lat0=None, lon0=None, chi=111.1, return_format='dict')
			lat, lon = xy_to_lat_lon(x,y, sim_file=sim_file, lat0=None, lon0=None, chi=111.1, return_format='tuple')
			#
			# submit year in days. ETAS will internally convert these values to seconds...
			vc_etas_catalog += [[float(rw['event_year']*365.25), lat, lon, float(rw['event_magnitude']), z]]
			#   #note: these float() conversions look stupid, but they do something to decouple the h5py object from the row (make a copy,
			#   # not a reference probably), which makes the recarray outputs not totally screwed up.
		#
		# okada_disps = numpy.rec.array(okada_disps, names=['x', 'y', 'z', 'dx', 'dy', 'dz'], formats=[type(x).__name__ for x in okada_disps[0]])
		# field_data = numpy.core.records.fromarrays(zip(*field_data), names=['x', 'y', 'z', 'dx', 'dy', 'dz', 'dxyz', 'dxy'], formats=[type(x).__name__ for x in field_data[0]])
		#
	#
	#vc_etas_catalog_ary=None
	try:
		# note: gotta convert the source data to an array or you get, for each field, something like (x,0., 0., 0., 0.) instead of just x.
		# or we can use the numpy.core.records.fromarrays() method...
		# both of these work:
		#vc_etas_catalog = numpy.rec.array(numpy.array(vc_etas_catalog), names = ['event_year', 'lat', 'lon', 'mag', 'depth'], formats = [type(x).__name__ for x in vc_etas_catalog[0]])
		vc_etas_catalog = numpy.core.records.fromarrays(zip(*vc_etas_catalog), names = ['event_year', 'lat', 'lon', 'mag', 'depth'], formats = [type(x).__name__ for x in vc_etas_catalog[0]])
	except:
		print "unable to convert to rec_array..."
		pass
	#
	return vc_etas_catalog, sweeps_index
	#


# CFF and catalog management helper functions.
#
# def mean_recurrence(ary_in='data/VC_CFF_timeseries_section_123.npy', m0=7.0, do_plots=False, do_clf=True):
def mean_recurrence(ary_in='data/VC_CFF_timeseries_section_123.npy', m0=7.0, do_plots=False, do_clf=True, use_precalc=True, section_ids=[], h5data=default_sim_file):
	# find mean, stdev Delta_t, N between m>=m0 events in ary_in.
	# for now, assume ary_in is a structured array, h5, or name of a structured array file.
	#
	# TODO: we want to test our pre-calcs against optimized array loading/searches, aka using [()] syntax.
	#if use_precalc==False and section_ids!=None and (hasattr(section_ids, '__len__') and len(section_ids)>0):
	if use_precalc==False and section_ids!=None:
		#if not hasattr(section_ids, '__len__'): section_ids=[section_ids]
		#if isinstance(section_ids, str): list(eval(section_ids))		# give it a shot anyway...
		#
		print "not-precalcing"
		ary_in = get_event_time_series_on_section(sim_file=h5data, section_id=section_ids, n_cpus=None, is_sorted=False, recarray=True)
	#
	if isinstance(ary_in, str)==True:
		ary_in = numpy.load(ary_in)
	# ... and for good measure, let's always sort these (just in case)... or at least try to.
	try:
		ary_in.sort(order='event_year')
	except:
		print "ary_in sorting failed, but continuing anyway, hoping for not nonsense..."
	#
	if not hasattr(ary_in, 'dtype') and isinstance(ary_in[0][0], str):
		# probably a cols + data list format. convert to recarray:
		ary_in = numpy.core.records.fromarrays(zip(*ary_in[1:]), names=ary_in[0], formats=[type(x).__name__ for x in ary_in[0]])
	#
	Ns, Js, Ts, Ms = zip(*[[x['event_number'], j, x['event_year'], x['event_magnitude']] for j, x in enumerate(ary_in) if float(x['event_magnitude'])>=m0])
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
	stress_drop = [(x['cff_initial']-x['cff_final'])**2. for x in ary_in[2:] if float(x['event_magnitude'])>=m0]	# large m stress-drop
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
	'''
	# get trend data for ary_in.
	# (re-write and clean-up from get_trend_analysis() so we don't make mistakes...
	#
	# main event: fit sub-sequences of length nyquist_len to a linear fit; keep the a,b fits. not that numpy.linalg.lstsq()
	# return the fit as (slope, intercept), which has created a bit of a mess.
	#
	# maybe also return record-breaking stats., though i don't think we're using them just yet.
	#
	# note: we're working with three tiers of the input data set: the raw set, the intervals set, and the fit-set. the indices line up like:
	# j_intervals --> j_raw-1
	# j_fit --> j_raw + nyquist len, aka fit[0][col] = raw[nyquist_len][col]
	# basically, we lose the first entry from the raw set when we take the intervals dt = t_i - t_{i-1}
	# then, we lose nyquist_len-1 when we do the line-fit over nyquist_len
	'''
	if isinstance(ary_in, str)==True:
		# ... and these are pre-calculated segment-time-series type recarrays, like:
		# ary_in = 'data/VC_CFF_timeseries_section_125.npy'
		ary_in = numpy.load(ary_in)
	#
	#intervals = numpy.array([ary_in['event_year'][i]-ary_in['event_year'][i-1] for i in xrange(1, len(ary_in))])
	# one array like: [ [N, t,dt] ]
	intervals = zip(ary_in['event_number'][1:], ary_in['event_year'][1:], ary_in['event_year'][1:]-ary_in['event_year'][:-1] )
	#
	#mean_interval = (max(ary_in['event_year'])-min(ary_in['event_year']))/float(len(ary_in['event_year'])-1)	# "-1" --> interval, not event.
	#
	#output_names = ['event_number', 'event_year', 'lin_fit_a', 'lin_fit_b', 'rb_ratio', 'interval_rate_ratio', 'intervals']
	output_names = ['event_number', 'event_year', 'intervals', 'lin_fit_b', 'lin_fit_a']
	# for convenience, wrap it up:
	intervals = numpy.core.records.fromarrays(zip(*intervals), names=output_names[0:3], formats=[type(x).__name__ for x in intervals[0]])
	#
	#fitses = [numpy.linalg.lstsq(X_numpy[i-nyquist_len:i], intervals[i-nyquist_len:i])[0] for i in xrange(nyquist_len, len(ary_in))]
	#
	#fitses = [list([Ns[i-1], Xs[i-1]]) + list(numpy.linalg.lstsq(X_numpy[i-nyquist_len:i], intervals[i-nyquist_len:i])[0]) for i in xrange(nyquist_len, len(intervals))]
	#
	#
	# create dictionary fits. this is a little bit inefficient, but it makes it easier to keep test and keep track of during diagnostics.
	#
	'''
	# expanded for-loop format, and see list-comprehension format below.
	fits = []
	for j,rw in enumerate(intervals[nyquist_len-1:]):
		# so we have some mixed indices. j is the index for this counter. our desired sequence range is [j:j+nyquist_len], so if nyquist_len=5,
		# we'll get indices 0-4. rw[col] = intervals[j+nyquist_len-1][col], so we defind k = j + nyquist_len -1
		k=j+nyquist_len-1
		this_fit = numpy.linalg.lstsq([[x,1.0] for x in intervals['event_year'][j:j+nyquist_len]], intervals['intervals'][j:j+nyquist_len])[0]
		#fits += [ [rw['event_number'], intervals['event_number'][k], rw['event_year'], intervals['event_year'][k], intervals['intervals'][k], rw['intervals'], this_fit[1], this_fit[0]] ]
		#
		fits += [{key:val for key,val in zip(output_names, [rw['event_number'], rw['event_year'], rw['intervals']] +list(this_fit))}]
	'''
	# so notes on calculating fitses (when we check this for the upteenth time...:
	# fitses is a list of dicts; each dict. contains keys from output_names above; looks like
	# [{'event_number': xx, 'event_year': xx , 'intervals': (which is just one interval value), 'lin_fit_b': xx, 'lin_fit_a': xx} ]
	# note the whole set spins through intervals[nyquist_len-1:], so the first value we can fit and forward.
	# then, the index j starts at 0, so an interval is like rw['event_year']-intervals[j]['event_year'] (aka, j starts at 0; but
	# the intervals[] index (say k) starts at nyquist_len-1.
	#
	fitses = [{key:val for key,val in zip(output_names, [rw['event_number'], rw['event_year'], rw['intervals']] + list(numpy.linalg.lstsq([[x,1.0] for x in intervals['event_year'][j:j+nyquist_len]], intervals['intervals'][j:j+nyquist_len])[0]))} for j,rw in enumerate(intervals[nyquist_len-1:])]
	#
	
	#
	# get record-breaking intervals
	# (are we using these anywhere? for now, comment them out. we can add them back if we need them.
	'''
	nrbs=[]
	#output_names += ['rb_ratio']
	for i in xrange(nyquist_len, len(intervals)):
		#
		nrbs_up=[intervals['intervals'][i-nyquist_len]]
		nrbs_dn=[intervals['intervals'][i-nyquist_len]]
		#
		#for j, interval in enumerate(intervals[i-nyquist_len:i]):
		for j, interval in enumerate(intervals['intervals'][i-nyquist_len+1:i]):
			# note on "+1": we've already assigned/evaluated i-nyquist_len.
			if interval > nrbs_up[-1]: nrbs_up+=[interval]
			if interval < nrbs_dn[-1]: nrbs_dn+=[interval]
		#
		nrbs += [math.log10(float(len(nrbs_up))/float(len(nrbs_dn)))]
		#rb_ratio = float(len(nrbs_up))/float(len(nrbs_dn))
		#outputs[i]+=[rb_ratio]
		#
	#nrbs = [None for i in xrange(nyquist_len)] + nrbs
	'''
	#
	print "lens: ", len(fitses), len(intervals)
	#
	outputs = numpy.core.records.fromarrays(zip(*[ [fit[key] for key in output_names] for j,fit in enumerate(fitses) ]), names=output_names, formats = [type(fitses[0][key]).__name__ for key in output_names])
	#

	return outputs	
#
def plot_trend(CFF_in=None, section_ids=[16], nyquist_len=10, fignum=0, do_clf=True, m0=7.0):
	# diagnostic tool for CFF interval fitting.
	# CFF in: a CFF array. for simpler operation (and reference), provide a section_id; the script will pull a CFF from
	# the standard place.
	# nyquist_len: length of subsequence to sample.
	#
	if CFF_in==None:
		if not hasattr(section_ids, '__len__'): section_ids=[section_ids]
		#
		CFF_in=combine_section_CFFs(sections=section_ids)
	#
	intervals = CFF_in['event_year'][1:]-CFF_in['event_year'][:-1]
	trends = get_trend_analysis(ary_in=CFF_in, nyquist_len=nyquist_len)
	#
	f=plt.figure(fignum)
	if do_clf: plt.clf()
	ax_main = f.add_axes([.1, .45, .85, .5])
	ax_main.set_xlabel('time $t$ (years)')
	ax_main.set_ylabel('Interval $\\Delta t$')
	ax_b    = f.add_axes([.1, .05, .85, .3], sharex=ax_main)
	ax_b.set_ylabel('Interval slope $b$')
	#
	ax_main.plot(CFF_in['event_year'][1:], intervals, '.-')
	XX_t, YY_t = [], []
	for i,t in enumerate(trends):
		# lines over which data are fit:
		X_t = [CFF_in['event_year'][i+1], CFF_in['event_year'][i + nyquist_len]]
		Y_t = [t['lin_fit_a'] + t['lin_fit_b']*X_t[-2], t['lin_fit_a'] + t['lin_fit_b']*X_t[-1]]
		#
		#ax_main.plot(X_t, Y_t, '.-')
		#XX_t += [CFF_in['event_year'][i + nyquist_len -1], CFF_in['event_year'][i + nyquist_len]]
		#YY_t += [t['lin_fit_a'] + t['lin_fit_b']*X_t[-2], t['lin_fit_a'] + t['lin_fit_b']*X_t[-1]]
		XX_t += [CFF_in['event_year'][i + nyquist_len]]
		YY_t += [t['lin_fit_a'] + t['lin_fit_b']*X_t[-1]]
		#
	ax_main.plot(XX_t, YY_t, 'g-')
	ax_b.plot(trends['event_year'], trends['lin_fit_b'], '.-', zorder=3)
	ax_b.plot([trends['event_year'][0],trends['event_year'][-1]], [0., 0.], 'k', zorder=0)
	ax_main.plot([trends['event_year'][0],trends['event_year'][-1]], [0., 0.], 'k', zorder=0)
	#
	# plot large events:
	max_dt = max(intervals)
	min_dt = min(intervals)
	max_b  = max(trends['lin_fit_b'])
	min_b = min(trends['lin_fit_b'])
	#
	for rw in CFF_in:
		if rw['event_magnitude']<m0: continue
		ax_main.plot([rw['event_year'], rw['event_year']], [min_dt, max_dt], 'r-', lw=2, alpha=.7, zorder=4)
		ax_b.plot([rw['event_year'], rw['event_year']], [min_b, max_b], 'r-', lw=2, alpha=.7, zorder=4)
#
def plot_CFF_ary(ary_in=None, fnum=0, nyquist_factor=.5, gt_lt_eval=operator.lt, section_id=16):
	'''
	# (depricate? this bit can probably be removed. it was originally intended to investigate the usefulness of plotting CFF
	# time series... which appear to be mostly useless. subsequent figure-scripts plotting  interval \Delta t and d(\Delta t)/dt
	# apper more useful, and the alert system appears to be a little bit more complicated with nyquist_factor and b_0, so this script might 
	# be more clutter than useful at this stage.
	#
	# this script is for some of the earlier CFF numpy array types. newer arrays may require different scripting. (but it seems to work ok).
	# basics: plots the triple time-series figure with [intervals, magnitudes, alert], [CFF, Delta CFF (?)], and [mags,slopes] ].
	# this was a development figure where be basically decided not to use CFF. and the alert figure got moved to ???.
	# for older data sets (most likely):
	#	# 2 cols: event_number, CFF_initial
	#   # 3 cols: event_number, event_year, CFF_initial
	#   # 4 cols: event_number, event_year, CFF_initial, CFF_final
	#  later versions will include column names.
	#
	# ary_in: can be a numpy.recarray or a string/filename like: 'data/VC_CFF_timeseries_section_125.npy'
	#   numpy array would be like numpy.load(ary_in)
	#
	# gt_lt_eval: how do we evaluate greter/less than? submit gt/lt functions: operator.{lt, le, gt, ge} or
	# strings: {'le', 'leq', 'ge', 'geq'}, and we trap a few others as well.
	#
	# (this is primarily a development script and may need to be depricated).
	#
	'''
	#
	# are we going to be a greater to less than camp?
	#gl_lt_eval = operator.lt		# "less than", "less that or equal", "greater than", "greater than or equal"
	#gl_lt_eval = operator.le
	#gl_lt_eval = operator.gt
	#gt_lt_eval = operator.ge		# greater than or equal to; ge(4,5)=False, ge(4,4)=True, ge(4,3)=True
	# and allow for string inputs:
	if isinstance(gt_lt_eval, str):
		gt_lt_eval = gt_lt_eval.lower()
		if gt_lt_eval in ('lt', 'less', 'lessthan', 'less_than'): gl_lt_eval = operator.lt
		if gt_lt_eval in ('lte', 'leq', 'lessthanorequal', 'less_than_or_equal'): gl_lt_eval = operator.le
		if gt_lt_eval in ('ge', 'greater', 'greaterthan', 'greater_than'): gl_lt_eval = operator.gt
		if gt_lt_eval in ('gte', 'geq', 'greaterthanorequal', 'greater_than_or_equal'): gl_lt_eval = operator.ge
	#
	if section_id!=None and (ary_in==None or len(ary_in)==0):
		# guess the ary_in file name from section_id (and later on, trap for multiple section_id values):
		ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id
		
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		CFF = ary_in
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
		# bottom:
		ax_mag = f.add_axes([.1, .05, .85, .25])
		#
		# CFF plot. (center):
		ax_CFF = f.add_axes([.1, .35, .85, .25], sharex=ax_mag)
		ax_dCFF = ax_CFF.twinx()	# over-plot stress (CFF) drop...
		#
		# top (intervals, etc.)
		ax_ints  = f.add_axes([.1, .65, .85, .25], sharex=ax_mag)
		ax_mag2  = ax_ints.twinx()
		ax_trend2 = ax_ints.twinx()
		#
		ax_trend = ax_mag.twinx()
		#
		X_init = CFF['event_year']
		X_finals = [x+.01 for x in X_init]	# this to show the stress drops after a mainshock (X_final is the time of the data point after the event).
		#
		Y0 = -1.*CFF['cff_initial']
		Y_final = -1.*(CFF['cff_final'])
		
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
		#ax_trend.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='r', zorder=2, alpha=.8)
		ax_trend.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], y2=[0.0 for x in trend_data], where=[gt_lt_eval(x['lin_fit_b'],0.) for x in trend_data], color='r', zorder=2, alpha=.8)
		ax_trend.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		ax_trend.set_ylabel('(log) interval slope $b$')
		#
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], 'r-', zorder=5, alpha=.8)
		#
		ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b']  for x in trend_data], y2=[0.0 for x in trend_data], where=[gt_lt_eval(x['lin_fit_b'],0.) for x in trend_data], color='m', zorder=1, alpha=.5)
		ax_trend2.fill_between([x['event_year'] for x in trend_data], [1.  for x in trend_data], y2=[0.0 for x in trend_data], where=[gt_lt_eval(x['lin_fit_b'],0.) for x in trend_data], color='m', zorder=1, alpha=.25)
		#		
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
def get_EMC_CFF(sections=None, file_out_root='data/VC_CFF_timeseries_section'):
	# make pre-calculated CFF,event time-series for EMC sections.
	# (wrapper script for get_CFF_on_section()  , so this maybe should move to vc_paper_emc_figs ?)
	#if sections==None: sections = [16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 69, 70, 73, 83, 84, 92, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 123, 124, 125, 126, 149]
	if sections==None: sections = emc_sections
	#
	event_sweeps_dict = None
	#with h5py.File(sim_file) as vc_data:
	#	# pre-calc a dictionary of event-sweeps. each event_number will include a list of sweep events.
	#	event_sweeps_dict = make_h5_indexed_dict_spp(h5_in=vc_data['event_sweep_table'], index_col='event_number')
	#
	for section in sections:
		X=get_CFF_on_section(section_id=section, event_sweeps_dict=event_sweeps_dict)
		X=numpy.array(X)
		#X.dump('data/VC_CFF_timeseries_section_%d.npy' % section)
		X.dump('%s_%d.npy' % (file_out_root, section))
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
def get_CFF_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, event_sweeps_dict=None, verbose=False):
	# CFF = shear_stress - mu*normal_stress
	#  blocks: dictionary of blocks indexed by event number: {<event_number>:[block_data, block_data, etc.], <>:[]}
	# for runs of multiple instances, this can be pre-calculated and shared.
	#
	# 1) get events for the section (use "events_by_section".. (?) get_event_time_series_on_section() (?) )
	# 2) for each event, gather all the blocks (use the event_sweep_table).
	# 3) for each block-sweep in the event, add the local CFF (aka, add the shear_init + mu*normal_init for each entry in the event.
	# 4) time (year) of the event can be taken from the event_table.
	#
	time_start = time.time()
	print "starting CFF(t) for section: %d (%s)" % (section_id, time.ctime())
	if n_cpus==None: n_cpus = mpp.cpu_count()
	CFF=[]
	# return as a list since we'll be adding columns.
	events = get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=section_id, n_cpus=n_cpus, recarray=False)
	t0=time.time()
	#
	# so, events is a time series of events (earthquakes) on a particular fault section.
	# for each event, the normal/shear stresses are calculated.
	#
	# now, fetch block/sweep level data. for now, let's return a simpler data set than we initially read.
	#
	with h5py.File(sim_file) as vc_data:
		vq_events = vc_data['event_table'][()]		# we might want a smart way of handling different-labeled data (aka, i think new vq is just ['events']
		sweep_data = vc_data['event_sweep_table']	# we won't be able to load the whole sweeps table, but events gives 'start_sweep', 'end_sweep' indices.
		max_index = len(events)-1					# where the full sequence is sweeps[start_sweep:end_sweep], aka, end_sweep = {actual end sweep} + 1
		#
		# make some mpp queues:
		#jobs     = []
		#pipes    = []
		#cff_outs = []
		for ev_count, event in enumerate(events[1:]):
			# this row indexing should be cleaned up to allow for a recarray, or similar, here. so we should sweep over events, not events[1:]
			# get_event_blocks(sim_file=allcal_full_mks, event_number=None, block_table_name='block_info_table')
			#
			# "events" is a list type, but each row is some sort of numpy recarray (returns numpy.void), but acts like rec_array
			event_id=event['event_number']
			event_year = event['event_year']
			event_area = float(event['event_area'])
			event_magnitude = event['event_magnitude']
			#
			# yoder:
			start_sweep = vq_events['start_sweep_rec'][event_id]
			end_sweep   = vq_events['end_sweep_rec'][event_id]
			#
			blocks = sweep_data[start_sweep:end_sweep][()]		# ... and i think this should be the fastest way to do this... and new vq will probably alwo require
																# a query to the fault model data.
			#
			if verbose and ev_count%250==0:
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
			# note: this will be a change in the CFF calculation, namely that we're going to weight each CFF contribution by block area and then average
			# by dividing by total area. the previous calc_CFFs_from_blocks() averages over the count of blocks or sweeps... which is probably meaningless,
			# but we never really did anything with it either, so it's probably ok... but also, all the pre-calculated sets have bogus CFF measures (at the event level).
			#
			# ... but otherwise, this method of calculation is simple, stable, and fast.
			cff_initial = sum(blocks['area']*(blocks['normal_init']-blocks['shear_init']*blocks['mu']))/event_area
			cff_final   = sum(blocks['area']*(blocks['normal_final']-blocks['shear_final']*blocks['mu']))/event_area
			#
			#
			#names=['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'event_area']
			# ... and this is WAY faster...
			CFF += [[event_id, event_year, event_magnitude, cff_initial, cff_final, event_area]]
			#
			
		#
	# convert CFF to a structured array (basically, we're after named columns; in this case, defining the data types is
	# easy as well).
	# this is what we'll want to do, but for some reason it doesn't work worth a damn. just send back the list for now...
	#CFF = numpy.array(CFF, dtype=[('event_number', 'int32'), ('event_year', 'int32'), ('event_magnitude', 'float64'), ('CFF_init', 'float64'), ('CFF_final', 'float64')])
	#
	if verbose: print "finished CFF(t) for section: %d (%s: %f)" % (section_id, time.ctime(), time.time()-time_start)
	#
	# convert to structured array with named cols (note this syntax, because this is not as easy as it should be):
	CFF = numpy.core.records.fromarrays(zip(*CFF), names=['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'event_area'], formats=[type(x).__name__ for x in CFF[0]])
	#
	return CFF
#
def get_stress_on_section(sim_file=allcal_full_mks, section_id=None, n_cpus=None, fignum=0):
	# ... and "time_series" is implied.
	# get time series of stress on a block. plot 2 axes on one plot: normal_stress, shear_stress.
	# note: the CFF, what we'll ultimately want, is:
	# CFF_j = shear_stress_j(t) - mu_j*normal_stress_j(t)
	#
	events = get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=section_id, n_cpus=n_cpus, recarray=True)
	#col_dict = {}
	#map(col_dict.__setitem__, events[0], range(len(events[0])))
	col_dict = {col:j for j,col in enumerate(events.dtype.names)}
	#
	# each time step has an initial + final stress value.
	ts_single = []
	mags = []
	ts = []
	shear_stress  = []
	normal_stress = []
	ave_slip = []
	for rw in events[1:]:
		#ts += [rw[col_dict['event_year']], rw[col_dict['event_year']]]
		#ts_single += [rw[col_dict['event_year']]]
		#mags += [rw[col_dict['event_magnitude']]]
		#shear_stress += [rw[col_dict['event_shear_init']], rw[col_dict['event_shear_final']]]
		#normal_stress += [rw[col_dict['event_normal_init']], rw[col_dict['event_normal_final']]]
		ts += [rw['event_year'], rw['event_year']]
		ts_single += [rw['event_year']]
		mags += [rw['event_magnitude']]
		shear_stress += [rw['event_shear_init'], rw['event_shear_final']]
		normal_stress += [rw['event_normal_init'], rw['event_normal_final']]
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
	upper_shear_init = get_peaks(events[1:], col='event_shear_init', peak_type='upper')
	upper_shear_final = get_peaks(events[1:], col='event_shear_final', peak_type='upper')
	#print "lens: ", len(upper_shear_init), len(upper_shear_final)
	#print upper_shear_init[-5:]
	#print upper_shear_final[-5:]
	lower_shear_init = get_peaks(events[1:], col='event_shear_init', peak_type='lower')
	lower_shear_final = get_peaks(events[1:], col='event_shear_final', peak_type='lower')
	#
	upper_normal_init = get_peaks(events[1:], col='event_normal_init', peak_type='upper')
	upper_normal_final = get_peaks(events[1:], col='event_normal_final', peak_type='upper')
	lower_normal_init = get_peaks(events[1:], col='event_normal_init', peak_type='lower')
	lower_normal_final = get_peaks(events[1:], col='event_normal_final', peak_type='lower')
	#
	ax = myfig.axes[0]
	ax.plot(ts, shear_stress, '.-', label='shear stress')
	for tbl, col,lbl in [(upper_shear_init, 'event_shear_init', 'init_upper'), (upper_shear_final, 'event_shear_final', 'final_upper'), (lower_shear_init, 'event_shear_init', 'init_lower'), (lower_shear_final, 'event_shear_final', 'final_lower')]:
		#ax.plot(map(operator.itemgetter(col_dict['event_year'], upper_shear_init)), map(operator.itemgetter(col_dict['event_shear_init'], upper_shear_init)),
		#print tbl[0]
		#
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
#
# end CFF calculators and helpers...
#
##############################################
# technical bits and helper functions:
#
def xy_to_lat_lon(x, y, sim_file=allcal_full_mks, lat0=None, lon0=None, chi=111.1, return_format='dict'):
	# for now, a simple x,y conversion. we should probably use at least a spherical distance formula.
	# chi: degress to km conversion.
	#
	# we can build in conditionals later, but for now, assume we need y, chi, lon0, x, lat0
	#
	if lon0==None or lat0==None:
		# get these from sim file:
		if sim_file==None:
			print "nothing to work with..."
			return None
		#
		with h5py.File(sim_file) as vc_data:
			lat0 = vc_data['base_lat_lon'][0]
			lon0 = vc_data['base_lat_lon'][1]
	#
	deg2rad = 2.0*math.pi/360.
	#
	lat = lat0 + (y/1000.)/chi
	lon = lon0 + (x/1000.)/(math.cos(deg2rad*lat)*chi)		# x,y standard in m?
	#
	if return_format=='tuple':
		return lat, lon
	elif return_format=='list':
		return [lat, lon]
	elif return_format=='lonlat_list':
		return [lon, lat]
	else:
		return {'lat':lat, 'lon':lon}
#
def lat_lon_to_xy(lat, lon, sim_file=allcal_full_mks, lat0=None, lon0=None, chi=111.1, return_format='dict'):
	#
	# these should probably also include a z component?
	#
	if return_format == None: return_format = 'tuple'
	#
	if lon0==None or lat0==None:
		# get these from sim file:
		if sim_file==None:
			print "nothing to work with..."
			return None
		#
		with h5py.File(sim_file) as vc_data:
			lat0 = vc_data['base_lat_lon'][0]
			lon0 = vc_data['base_lat_lon'][1]
		#
	#
	deg2rad = 2*math.pi/360.
	#
	y = (lat-lat0)*chi * 1000.
	x = (lon-lon0)*math.cos(lat*deg2rad)*chi * 1000.		# xy, standard is in meters?
	#
	if return_format=='tuple':
		return x,y
	elif return_format=='list':
		return [x,y]
	else:
		return {'x':x, 'y':y}
#
def in_rec_array(rec_array, col_name, in_list, pipe=None):
	# wrapper for MPP finding sub-sets in list. note this requires more specialized preparation than a simple :in_row_list()
	# type function (to which we can map), but by containing the list comprehension, this should be pretty fast.
	# is this the same as find_in()? not quite. find_in() is designed to use a numerical index and use pipe() to
	# communicate back to a Process()
	#
	if pipe==None:
		return [rw for rw in rec_array if rw[col_name] in in_list]
	else:
		# you can use this with a Process() and communicate with Pipe()s, but it's pretty slow (same as
		# list comprehension in-line). somehow, porting this to a Pool() and using apply_async() runs 100x faster...
		# even with just the single cpu. in fact, the 1-cpu implementation appears to run faster than n_cpus.
		#
		pipe.send([rw for rw in rec_array if rw[col_name] in in_list])
#
def in_rec_test(section_ids=None, sim_file=allcal_full_mks, n_cpus=None):
	# results summary:
	# 1) filter(), with a lambda: function, is a little bit slower than the list comprehension
	# 2) Process() and Pipe() is pretty slow -- about the same as the in-line list comprehsnsion (if i recall, it's
	#    the Pipe() part that is really slow. might be able to improve this by using direct reference to a structured array.
	# 3) the list-comprehension approach is pretty slow as well (all of these 9-10 seconds for EMC faults.
	# 4) using Pool() and apply_async() is almost 10x faster (about 1.1 seconds).
	if n_cpus==None: n_cpus = mpp.cpu_count()
	#
	section_ids = (section_ids or emc_section_filter['filter'])
	print "section_ids: ", section_ids
	#
	with h5py.File(sim_file) as vc_data:
		if True:
			print "len_0: %d" % len(vc_data['block_info_table'])
			t0=time.time()
			print t0
			block_info = numpy.array([rw for rw in vc_data['block_info_table'] if rw['section_id'] in section_ids], dtype=vc_data['block_info_table'].dtype)
			print "bi_len: %d" % len(block_info)
			#
			print time.time(), " :: ", time.time()-t0
			t0=time.time()
			#
			#pipe_in, pipe_out = mpp.Pipe()
			#proc = mpp.Process(target=in_rec_array, kwargs={'rec_array':vc_data['block_info_table'], 'col_name':'section_id', 'in_list':section_ids, 'pipe':pipe_out})
			#proc.start()
			#proc.join()
			#
			#block_info = numpy.array(pipe_in.recv(), dtype=vc_data['block_info_table'].dtype)
			block_info = numpy.array(filter(lambda rw: rw['section_id'] in section_ids, vc_data['block_info_table']), dtype=vc_data['block_info_table'].dtype) 
			print "bi_len: %d" % len(block_info)
			#
			print time.time(), " :: ", time.time()-t0
		if True:
			t0=time.time()
			print t0
			pool = mpp.Pool(n_cpus)
			results = []
			tbl = vc_data['block_info_table']
			N_len = len(tbl)
			dN = int(N_len/n_cpus)
			#for i in xrange(n_cpus):
			#
			for i in xrange(n_cpus):
				#
				N_term = (1+i)*dN
				if i==N_len-1: N_term = N_len
				results+=[pool.apply_async(in_rec_array, args=(), kwds={'rec_array':tbl[(i)*dN:N_term], 'col_name':'section_id', 'in_list':section_ids})]
				#
			pool.close()
			pool.join()
			#
			if len(results)>1:
				block_info2 = numpy.array(reduce(numpy.append, [x.get() for x in results]), dtype=tbl.dtype)
			else:
				# reduce() will throw an error if you give it only one value.
				block_info2 = results[0].get()
			#
			print "bi_len: %d" % len(block_info)
			print time.time(), " :: ", time.time()-t0
	return block_info, block_info2
#
def make_structured_arrays(file_profile = 'data/VC_CFF_timeseries_section_*.npy'):
	# wrapper to convert a bunch of normal arrays or maybe lists to numpy structured arrays (numpy.recarray).
	# (aka, this is a one time, or at least special occasion, script. with the exception that it might be instructive
	# for developing new array handling scripts, it can probably be chucked).
	G=glob.glob(file_profile)
	#
	for g in G:
		print "fixing file: %s" % g
		try:
			z=make_structured_array_from_file(fname_in=g, fname_out=g)
			print "successful..."
		except:
			print "failed to 'fix' file %s. it might have already been converted..." % g
#
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
	#col_names=[]
	#for key in cols_dict.keys():
	#	col_names += [[key, cols_dict[key][1]]]
	#
	col_names = [[key,cols_dict[key][1]] for key in cols_dict.iterkeys()]
	col_names.sort(key=lambda x:x[1])
	#
	return [x[0] for x in col_names]

def _worker_in_table(src_table, test_table, test_col, target_Q):
	#
	# note we must be passed a 1-D vector test_table. 
	target_Q.put([x for x in src_table if x[test_col] in test_table])
	#return [x for x in src_table if src_table[test_col] in test_table]
#
#
def f_inv_weibull(P=None, t0=None, tau=None, beta=None):
	# inverse (solving for t or delta_t):
	#print "prams: ", P,t0,tau,beta
	P,t0,tau,beta = [float(x) for x in [P,t0,tau,beta]]
	#
	return tau*((t0/tau)**beta - math.log(1.-P))**(1./beta)
#
def f_weibull(x=None, chi=1.0, beta=1.0, x0=None):
	'''
	# weibull distribution (for fitting).
	# if different parameter ordering is necessary, as per specifics of fitting routine or something like that,
	# use an anonymous "lambda" function. aka, the "func" pram: fitter(f_wieibul, prams) -->
	# fitter(lambda x, chi, beta: f_weibull(x, chi, beta, my_x0_value), prams ...
	# so, "lambda" takes the parameters x, chi, beta and passes them to f_weibull along with the externally defined my_x0_value.
	# now, the fitting function looks like ([x], (variable prams) ), as expected by curve_fit().
	'''
	if x0==None: x0=0.0		# ... but we give it None so the fitting algorithm does not try to fit...
	#
	return 1.0 - numpy.exp( ((x0/chi)**beta) - ((x/chi)**beta))
#
def mcfitter(func=f_weibull, X=None, Y=None, prams_dict=None, nits=10**6, n_cpus=None):
	# quick mc fitter MC fitter for tough functions.
	# func: function fitting to
	# X,Y: X and Y data
	# prams dict is like {'pram_name':[min, max], ...}
	#
	if n_cpus==None:
		try:
			n_cpus = mpp.cpu_count()
		except:
			n_cpus = 1
	#
	if n_cpus>1:
		# do MPP. parse into SPP jobs. call this function (quasi-recursively) with n_cpus = 1.
		print "passing..."
		#
		print "doing mpp MC fit..."
		if n_cpus==None: n_cpus = mpp.cpu_count()
		#
		results=[]
		pool = mpp.Pool(n_cpus)
		for i in xrange(n_cpus):
			results+=[pool.apply_async(mcfitter, kwds={'func':f_weibull, 'X':X, 'Y':Y, 'prams_dict':prams_dict, 'nits':nits/n_cpus, 'n_cpus':1})]
		pool.close()
		#pool.join()
		#
		#Zmin, Z=[], []
		Zmin, Z = results[0].get()
		 
		for j, result in enumerate(results):
			if j==0: continue
			#
			a,b = result.get()
			Zmin = numpy.append(Zmin, a)
			Z    = numpy.append(Z, b)
			#
		Zmin.sort(order='chi_sqr')
		#
		return Zmin, Z
	#
	else:
		#
		# do SPP
		pram_handler = {}
		for key in prams_dict.keys():
			pram_handler[key] = {'min':prams_dict[key][0], 'delta':prams_dict[key][1]-prams_dict[key][0]}
			#
			if pram_handler[key]['delta']==0: 
				# choosing random numbers is expensive. if delta is zero, spoof random.Random() with "getter()" class
				# that quickly return 0.
				pram_handler[key]['rand']=getter(rand_val = 0.0)
			else:
				pram_handler[key]['rand']=random.Random()
		#
		this_prams = []
		min_pramses = []	# list of min. parameters
		pramses = []
		N_x = len(X)
		ndof = float(len(X) - len(prams_dict))
		#
		for i in xrange(nits):
			this_prams_dict = {key:pram_handler[key]['min'] + pram_handler[key]['delta']*pram_handler[key]['rand'].random() for key in prams_dict.keys()}
			#print this_prams_dict
			chi_sqr = numpy.sum([(func(X[j], **this_prams_dict)-Y[j])**2. for j in xrange(N_x)])/ndof
			#
			#print chi_sqr, i
			pramses += [[i] + [this_prams_dict[key] for key in this_prams_dict.keys()] + [chi_sqr]]
			if i==0 or chi_sqr < min_pramses[-1][-1]:
				min_pramses += [pramses[-1]]
			#
		# numpy.core.records.fromarrays(zip(*best_fit_array), names = ['section_id', 'tau', 'beta', 'sigma_tau', 'sigma_beta', 'mean_chi_sqr'], formats = [type(x).__name__ for x in best_fit_array[0]])
		#
		min_pramses = numpy.core.records.fromarrays(zip(*min_pramses), names = ['index'] + [key for key in prams_dict.keys()] + ['chi_sqr'], formats = [type(x).__name__ for x in min_pramses[0]])
		pramses = numpy.core.records.fromarrays(zip(*pramses), names = ['index'] + [key for key in prams_dict.keys()] + ['chi_sqr'], formats = [type(x).__name__ for x in pramses[0]])
		
		#
		return min_pramses, pramses	
#		
	
	
			

