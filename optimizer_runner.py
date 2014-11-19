import vc_parser
import numpy
import sys

#datas = vc_parser.optimize_metric_1(b_min=-.15, b_max=.1, d_b=.01, nyquist_min=.2, nyquist_max=.8, d_nyquist=.01, nits=100000)

def main(b_min=-.15, b_max=.1, d_b=.01, nyquist_min=.1, nyquist_max=.8, d_nyquist=.01, nits=100000, dump_file='dumps/optimize_faultwise_trig_105'):
	#datas = vc_parser.optimize_metric_faultwise(b_min=-.15, b_max=.1, d_b=.01, nyquist_min=.1, nyquist_max=.8, d_nyquist=.01, nits=100000, dump_file='dumps/optimize_faultwise_trig_105')
	datas = vc_parser.optimize_metric_faultwise(b_min=b_min, b_max=b_max, d_b=d_b, nyquist_min=nyquist_min, nyquist_max=nyquist_max, d_nyquist=d_nyquist, nits=nits, dump_file=dump_file)

# [b, nyq, mean, stdev]

#datas = numpy.core.records.fromarrays(zip(*datas), names=['b', 'nyq', 'mean', 'stdev'], formats = [type(x).__name__ for x in datas[0]])

#datas.dump('optimization_test.npy')
def prnt(*args):
	for arg in args: sys.stdout.write(str(arg))

if __name__ == '__main__':
	# assume standard arg, kwarg format (known args in proper order, then kwargs as desired). parse kwargs based on "=".
	#
	expected_types = [float, float, float, float, float,float, int, str]
	expected_types_names = ['b_min', 'b_max', 'd_b', 'nyquist_min', 'nyquist_max', 'd_nyquist', 'nits', 'dump_file']
	expected_types_dict = {expected_types_names[i]:expected_types[i] for i,x in enumerate(expected_types)}
	#
	#arguments = [arg for arg in sys.argv]
	#
	these_args = []
	these_kwargs = {}
	status_ok=True
	print "\n\n"
	#
	for i,arg in enumerate(sys.argv[1:]):
		# note: first arguement is the module name.
		print "raw arg: %s" % arg
		if '=' not in arg:
			# it's a kwarg.
			if len(these_kwargs)>0:
				print "unnamed arg cannot follow kwarg. doh!"
				#return None
				status_ok=False
				break
			#
			these_args += [expected_types[i](arg)]
			#
		if '=' in arg:
			#key, val = [str.strip(x) for x in arg.split('=')]
			these_kwargs.update({key:expected_types_dict[key](val) for key, val  in [[str.strip(x) for x in arg.split('=')]]})
		#
	#
	if status_ok:
		print "args: ", [prnt(arg , " :: " , type(arg).__name__, "\n ") for arg in these_args]
		print "kwargs: ", [prnt(key, " :: ", val, " :/: ", type(val).__name__, "\n ") for key,val in these_kwargs.iteritems()]
		print "***"
		print these_args
		print "**"
		print these_kwargs
		print "/n*******/n"
		rval = main(*these_args, **these_kwargs)
	else:
		print 'exited with error...'
