import vc_parser
import numpy

datas = vc_parser.optimize_metric_1(b_min=-.15, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01, nits=3000)

# [b, nyq, mean, stdev]

datas = numpy.core.records.fromarrays(zip(*datas), names=['b', 'nyq', 'mean', 'stdev'], formats = [type(x).__name__ for x in datas[0]])

datas.dump('optimization_test.npy')
