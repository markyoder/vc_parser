import vc_parser
import vc_geodetic

import pylab as plt
import datetime as dtm
import numpy
import random


def vc_markov_lambda(sections=vc_parser.emc_sections, n_bins=250, start_year=10000., end_year=None, catalog=None):
	'''
	#
	# look for this makov-lambda distribution in vc time-series. basically, we're looking for asymmetry in the 
	# distribution (histogram) of interval-differences (delta_t_2 - delta_t_1) in time series.
	#
	# ... so far, this is inconclusive...
	#
	'''
	#
	# permit the option to provide a catalog (speed things up a bit...)
	if catalog==None:
		catalog = vc_parser.combine_section_CFFs(sections=sections, start_year=start_year, end_year=end_year)
	print "catalog(s) combined... ", len(catalog)
	#
	delta_ts = [[x,x-catalog['event_year'][i]] for i,x in enumerate(catalog['event_year'][1:])]
	#
	delta_delta_ts = [[x[0], x[1] - delta_ts[i][1]] for i,x in enumerate(delta_ts[1:])]
	#
	plt.figure(0)
	plt.clf()
	plt.plot([x[0] for x in delta_ts], [x[1] for x in delta_ts], '.-')
	#
	plt.figure(1)
	plt.clf()
	#clrs = ['r' if x[1]<0 else 'b' for x in delta_delta_ts]
	clrs = [0 if x[1]<0 else 1 for x in delta_delta_ts]
	#plt.plot([x[0] for x in delta_delta_ts], [x[1] for x in delta_delta_ts], '-', color='b')
	plt.plot([x[0] for x in delta_delta_ts if x[1]>0], [x[1] for x in delta_delta_ts if x[1]>0], '.', color='b')
	plt.plot([x[0] for x in delta_delta_ts if x[1]<=0], [x[1] for x in delta_delta_ts if x[1]<=0], '.', color='r')
	#
	plt.figure(2)
	plt.clf()
	plt.hist([x[1] for x in delta_delta_ts], bins=n_bins, log=True)
	
