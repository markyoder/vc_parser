'''''''''''''''''''''''''''''''''''''''''''''
# vc_geodidic.py
# code property of University of California Davis, Dept. of Physics
# Rundle research group
# additional affiliations:
#  NASA Jet Propulsion Laboratory
#  NASA QuakeSim and E-DECIDER missions
#
# author: Mark R. Yoder, PhD
# email:  mark.yoder@gmail.com
#         mryoder@ucdavis.edu
#
# summary:
# produce various forms of (simulated) geodetic data, figures, etc. from the Virtual California (VC)
# earthquake simulator (or more generally, the Virtual Quake (VQ) earthquake simulator framework). The 
# basic strategy will be to infer surface displacement, analogous to GPS station displacements, from
# VC event-slip data.
#
'''''''''''''''''''''''''''''''''''''''''''''

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
import numpy
import scipy
import scipy.optimize as spo
import operator
import glob
import random
#
#import cStringIO
import sys
import json
import cPickle
import time
import os
import datetime as dtm
import pytz
#
#import imp
#import inspect		# use this module to determine a function's parameters and similar tasks.
import multiprocessing as mpp
#
#import ANSStools
#import BASScast



