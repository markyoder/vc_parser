'''
# some super simple scripts to get started reading VQ output data (using HDF5). 
# some of these scripts may be taylored towards an older data format, but the basic
# frameworks should apply.
#
# of course, also see the PyVQ tools available from GitHub. newer data sets will also be made
# available
'''

import h5py

def function_1(sim_file):
	# vq output are either CSV or HDF5. we'll focus on HDF5 since it's faster and better in every way.
	# the python library of interest is h5py. access a vq data file like this:
	#
	# contemporary "file open" way:
	# vc_data = h5py.File(sim_file)
	# # do stuff...
	# # more stuff...
	# vc_data.close()
	#
	# the problem with this approach is that if you forget to .close() the file or you have an error in your
	# script before .close() fires, the file is left open, locked, inaccessible, and makes a big mess (and 
	# unlike working with regular text files, you'll really experience this problem). instead, use "with":
	#
	with h5py.File(sim_file) as vc_data:
		# get the section_id values from block info:
		section_ids = set(vc_data['block_info_table']['section_id'])
	#
	print "section_ids in data set: %s" % ', '.join([str(s_id) for s_id in section_ids])
		
