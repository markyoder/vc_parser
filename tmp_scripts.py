# some temporary maintenance scripts....
import numpy
import glob
import os
import pylab as plt

def fix_roc_files():
	# outcome: only one dupe; they look the same.
	# in some cases, roc optimizations were re-run with n=1001 instead of n=1000. some of these might be to check the original run (was it corrupted?).
	#
	# first, let's look at all these cases. then, if the runs are ok, let's just combine the two files (so these lucky few will actually
	# be 2001 runs, not 1000.
	#
	# first, get a list of all the 1001 files (note, we'll want to change both their "all_prams" and best-prams versions 
	#(though we could probably just get rid of the best prams since they're just the best entry of the full-prams). also, 
	# rename the figure if it exists.
	fignum=0
	#
	fls = glob.glob('dumps/gji_roc_lt_500/roc_sec_lt_*_1001_allprams.npy')
	#fls = glob.glob('dumps//gji_roc_gt_500/roc_sec_*_nits_1001_allprams.npy')
	#return fls
	
	for fl in fls:
		fl_2 = fl.replace('_1001_', '_1000_')
		print "fls: ", fl, fl_2
		if not os.path.isfile(fl_2):
			# we have just the one file...
			#os.system('mv %s %s' % (fl, fl_2))
			print "fls*: ", fl, fl_2
			continue
		#
		print "look at this one..."
		# otherwise, we need to compare them:
		f=plt.figure(fignum, figsize=(4,8))
		plt.clf()
		plt.title(fl)
		ax1 = f.add_axes([.1,.1, .4, .8])
		ax2 = f.add_axes([.55, .1, .4, .8])
		#
		datas_1 = numpy.load(fl)
		datas_2 = numpy.load(fl_2)
		#
		F,H = zip(*[[rw['F'], rw['H']] for rw in datas_1])
		ax1.plot(F,H, '.')
		ax1.plot(range(2), range(2), 'r-')
		#
		F,H = zip(*[[rw['F'], rw['H']] for rw in datas_2])
		ax2.plot(F,H, '.')
		ax2.plot(range(2), range(2), 'r-')
		#
		fignum+=1
		
	
	return fls
	
def rename_1001_files():
	#fl_path = 'dumps/gji_roc_lt_500'
	fls = glob.glob('dumps/gji_roc_lt_500/*_1001*')
	for fl in fls:
		alt_fl = fl.replace('_1001.', '_1000.')
		if not os.path.isfile(alt_fl):
			os.system('mv %s %s' % (fl, alt_fl))
		
	
