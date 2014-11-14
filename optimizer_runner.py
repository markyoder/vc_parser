import vc_parser
import numpy

#datas = vc_parser.optimize_metric_1(b_min=-.15, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01, nits=100000)

datas = def optimize_metric_faultwise(b_min=-.15, b_max=.1, d_b=.01, nyquist_min=.1, nyquist_max=.8, d_nyquist=.01, nits=100000, dump_file='dumps/optimize_faultwise_trig_105')

# [b, nyq, mean, stdev]

#datas = numpy.core.records.fromarrays(zip(*datas), names=['b', 'nyq', 'mean', 'stdev'], formats = [type(x).__name__ for x in datas[0]])

datas.dump('optimization_test.npy')
