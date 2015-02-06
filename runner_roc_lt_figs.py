# roc runner script. for now, simple and stupid...

# to run headless, use this:
import matplotlib
matplotlib.use('Agg')
import pylab as plt

import multiprocessing as mpp
import vc_parser

import vc_paper_emc_figs as vpp


def f1():
	#a=vpp.create_ROC_figs_LT_data(section_ids = vpp.vc_parser.emc_sections, nits=500, num_roc_points=100, output_dir='dumps/gji_roc_lt_500')
	a=vpp.create_ROC_figs_LT_data(section_ids = vpp.vc_parser.emc_sections[0:20], nits=1000, num_roc_points=100, output_dir='dumps/gji_roc_lt_500')

def f1a():
	#a=vpp.create_ROC_figs_LT_data(section_ids = vpp.vc_parser.emc_sections, nits=500, num_roc_points=100, output_dir='dumps/gji_roc_lt_500')
	a=vpp.create_ROC_figs_LT_data(section_ids = vpp.vc_parser.emc_sections[20:], nits=1000, num_roc_points=100, output_dir='dumps/gji_roc_lt_500')


def f2():
	a=vpp.create_ROC_figs_data(section_ids = vpp.vc_parser.emc_sections[0:20], nits=1000, num_roc_points=100, output_dir='dumps/gji_roc_gt_500')

def f2a():
	a=vpp.create_ROC_figs_data(section_ids = vpp.vc_parser.emc_sections[20:], nits=1000, num_roc_points=100, output_dir='dumps/gji_roc_gt_500')


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
	p3 = mpp.Process(target=f1a)
	jobs.append(p3)
	p3.start()
	#
	p4 = mpp.Process(target=f2a)
	jobs.append(p4)
	p4.start()
	#
	print "jobs started."
	p1.join()
	p2.join()
	p3.join()
	p4.join()
	#
	print "jobs finished."
	#
	
