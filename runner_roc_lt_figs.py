# roc runner script. for now, simple and stupid...

# to run headless, use this:
import matplotlib
matplotlib.use('Agg')
import pylab as plt

import vc_paper_emc_figs as vpp

a=vpp.create_ROC_figs_LT_data(section_ids = vpp.vc_parser.emc_sections, nits=2500, num_roc_points=100)

