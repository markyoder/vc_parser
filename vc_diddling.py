import vc_parser
import vc_geodetic

import h5py

import pylab as plt
import datetime as dtm
import numpy
import random
import math

import glob


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


def fix_section_arys():
	# quickly, fix section arays where event_number is assigned as a float, not an integer.
	G=glob.glob('data/VC_CFF_timeseries_section_*.npy')
	aryses=[]
	for g in G:
		X=numpy.load(g)
		#return X
		
		print "type: ", type(X['event_number'][0]), X['event_number'][0], X.dtype
		#
		new_event_num = [int(x) for x in X['event_number']]
		my_formats = ['int'] + [type(x).__name__ for j,x in enumerate(X[0]) if j>0]
		#print "formats: ", my_formats
		#print "names  : ", X.dtype.names
		cols = zip(*X)
		new_X = numpy.core.records.fromarrays([new_event_num]+cols[1:], names=X.dtype.names, formats = my_formats)
		#
		new_X.dump(g)
		#
		aryses+=[new_X]
	
	return aryses
	

def m8_probs(vc_data_file=vc_parser.default_sim_file, t0=158., m0=7., m1=8.0, mt=7.6, delta_t=30., beta=1.1):
	# get events:
	with h5py.File(vc_data_file) as vc_data:
		events = vc_data['event_table'][()]
	#
	mags = [[rw['event_year'], rw['event_magnitude']] for rw in events if rw['event_magnitude']>=m0]
	intervals = [[rw[0], rw[0]-mags[j][0]] for j,rw in enumerate(mags[1:])]
	#
	mean_int = numpy.mean(zip(*intervals)[1])	# sort of a stupid way to calculate this, but we'll get away with it because it's pretty small.
	#
	n_factor_1 = 10.**(m1-m0)
	n_factor_2 = 10.**(1.5*(m1-mt) + 1.0*(mt-m0))
	#
	mean_int_1 = mean_int*n_factor_1
	mean_int_2 = mean_int*n_factor_2
	#
	print "n_factors: ", n_factor_1, n_factor_2
	print "mean_int: ", mean_int, mean_int_1, mean_int_2
	#
	# now, assume Poisson statistics and tau = mean_int. we can use a weibull dist., but since beta_allcal~1, it won't make a big difference.
	#
	P1 = 1.0 - math.exp(-delta_t/(mean_int*n_factor_1))
	P2 = 1.0 - math.exp(-delta_t/(mean_int*n_factor_2))
	#
	# Baysean probabilities:
	Pb1 = 1.0 - math.exp((t0/mean_int_1)**beta - ((delta_t+t0)/(mean_int_1))**beta)
	Pb2 = 1.0 - math.exp((t0/mean_int_2)**beta - ((delta_t+t0)/(mean_int_2))**beta)
	#
	print "poisson probs: \n", P1, P2
	print "\nBaysean probabilities for t0=%f:\n %f, %f" % (t0, Pb1, Pb2)
	
	
	
	
	
	
	
	
	
