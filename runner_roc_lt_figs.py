# roc runner script. for now, simple and stupid...

# to run headless, use this:
import matplotlib
matplotlib.use('Agg')
import pylab as plt

import multiprocessing as mpp
import vc_parser

import vc_paper_emc_figs as vpp


def f1():
	a=vpp.create_ROC_figs_LT_data(section_ids = vpp.vc_parser.emc_sections, nits=2500, num_roc_points=100)

def f2():
	a=vpp.create_ROC_figs_data(section_ids = vpp.vc_parser.emc_sections, nits=2500, num_roc_points=100)

if __name__ == '__main__':
	#print "executing..."
	#
	jobs = []
	
	p1 = mpp.Process(target=f1)
	jobs.append(p1)
	p1.start()
	#
	p2 = mpp.Process(target=f2)
	jobs.append(p2)
	p2.start()
	#
	print "jobs started."
	p1.join()
	p2.join()
	#
	
