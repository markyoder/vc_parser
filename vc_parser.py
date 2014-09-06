import matplotlib.pyplot as plt
plt.ion()
#
import math
import h5py
import numpy
import scipy
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
		event_blocks = [x for x in block_info if x['section_id'] in section_ids]
	#
	return [block_cols] + event_blocks
	
#def get_stress_time_series(sim_file=allcal_full_mks, block_id):
#	pass

def get_event_time_series_on_section(sim_file=allcal_full_mks, section_id=None):
	with h5py.File(sim_file,'r') as vc_data:
		view_ebs = vc_data['events_by_section']
		event_ids_tbl = view_ebs['section_%d' % section_id]
		#
		event_ids_ary = numpy.array(numpy.zeros(event_ids_tbl.len()), dtype='int64')
		#event_ids_ary.shape(numpy.zeros(event_ids_tbl.len(), len(event_ids_tbl[0])))
		event_ids_tbl.read_direct(event_ids_ary)
		#
		# now, fetch event data:
		# this is a choice, computationally. we can do it in a list comprehension [], or we can index the source list (event_id values
		# are sorted, so once we find one event in the master list, we can stop looking for it in the event_ids_ary list.
		# list comprehension is simpler...
		section_events = [x for x in vc_data['event_table'] if x['event_number'] in event_ids_ary]
		#
	#
	#return event_ids_ary
	return section_events

def get_event_time_series_on_section_mp(sim_file=allcal_full_mks, section_id=None, n_procs=2):
	if n_procs==None: n_procs = mpp.cpu_count
	#
	with h5py.File(sim_file,'r') as vc_data:
		# for now, single-process this prelimnary step (no point i think in mpp)...
		view_ebs = vc_data['events_by_section']
		event_ids_tbl = view_ebs['section_%d' % section_id]
		#
		event_ids_ary = numpy.array(numpy.zeros(event_ids_tbl.len()), dtype='int64')
		#event_ids_ary.shape(numpy.zeros(event_ids_tbl.len(), len(event_ids_tbl[0])))
		event_ids_tbl.read_direct(event_ids_ary)
		#
		# now, fetch event data:
		# this is a choice, computationally. we can do it in a list comprehension [], or we can index the source list (event_id values
		# are sorted, so once we find one event in the master list, we can stop looking for it in the event_ids_ary list.
		# list comprehension is simpler...
		src_table = vc_data['event_table']	# table from which we draw data
		test_table = event_ids_ary			# table we test against -- are dat in this list?
		#
		# now, we'll need to chop this up for multi-processing... and figure out which type of mpp we want to use...
		chunk_size = len(src_table)/n_procs
		#results_queue = multiprocessing.Queue()
		workers = mpp.Queue()
		section_events = mpp.Queue()
		worker_sections = []
		worker_pool = mpp.Pool(processes=n_procs)
		#
		wrker_sections=mpp.Queue()
		for i in xrange(n_procs):
			i0=i*chunk_size
			i1 = (i+1)*chunk_size
			if i==(n_procs-1): i1 = len(src_table)
			print "process_%d: (%d, %d)" % (i, i0, i1)
			#
			wrker_sections.put(src_table[i0:i1])
			
			#section_events = [x for x in vc_data['event_table'] if vc_data['event_table']['event_number'] in event_ids_ary]
			#workers.put(mpp.Process(target=_worker_in_table, args=(src_table[i0:i1], test_table, 'event_number', section_events)))
			#workers[-1].start()
			#print "worker %d started" % i
		resultses = worker_pool.map(_worker_in_table, [wrker_sections, test_table, 'event_number', section_events])
		#for i in xrange(len(workers)):
		#	workers[i].join()
			
		#
	#
	#return event_ids_ary
	return section_events
'''
class _table_sorter(mpp.Process):
	def __init__(self,src_table, test_table, test_col_src, target_tbl):
		self.src_table = src_table
		self.test_table = test_table
		self 
'''
def _worker_in_table(src_table, test_table, test_col, target_Q):
	#
	# note we must be passed a 1-D vector test_table. 
	target_Q.put([x for x in src_table if x[test_col] in test_table])
	#return [x for x in src_table if src_table[test_col] in test_table]
		
	
	
			

