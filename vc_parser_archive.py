


'''
# archived bits from vc_parser and related modules. i just can't get used to deleting stuff and counting on GitHub to git it back. 
# archive here first... will these functions work if we always run from vc_parser and use import *??
'''

def forecast_metric_1(ary_in=None, m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False, fnum=0, set_name=None, f_gt_lt=operator.gt, section_id=16, over_ride=False):
	'''
	# depricated:
	# this can be tossed in exchange for the more generalized format:
	# x = psa_forecast_1(CFF)		# (and other forecast metrics can be substituted of course)
	# y = evaluate_alert_segments(x,CFF)
	# so, until we get rid of this not the over_ride=False default parameter. if True, then run the old code. otherwise, wrap the new bits.
	# 
	#
	# forecast based on seismic acceleration. specifically, for a set of data (a fault-section time-series),
	# find the slopes of the inter-occurrence intervals. slope<0 implies acceleration, implies hazard.
	# the evaluation metric will be basically a ROC type:
	# 1) total time "at alert" serves as "false alarm"
	# 2) fraction of "predicted" events. this will vary. for starters:
	#     - slope at time t_i of event or time t_i-1 of the previous event is below b_0 (typically, b_0=0.0).
	#     - the "_i-1" condition is fair because 1) we don't know the slope is changin until the even happens
	#       and 2) in a real scenario, we'd have aftershocks, which woluld turn the post-seismic slope positive.
	#
	# parameters:
	# m0: recurrence magnidue
	# b_0: critical slope; issue alert for f_gt_lt(b,b_0)==True, aka, b>b_0 (see below)
	# nyquist_factor: fraction of inter-event sequence length to sample
	# do_plot: do plots ?
	# fnum: plot on figure(fnum)...
	# set_name: a name to give the set if ary_in is an array, not a string, or if we want to specify it.
	# (we should use the plotting routine from vc_parser.plot_CFF_ary(), namely the top subfig. )
	# f_gt_lt: function for greater/lesser evaluation... use operator.{gt, ge, lt, le}, or any two parameter input.
	'''
	#
	#
	if (ary_in==None or len(ary_in)==0) and section_id!=None:
		ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		# otherwise, assume we've been passed a proper CFF object:
		CFF = ary_in
		# give the set a name  so we don't return the whole input data object...
		if set_name==None: set_name='unnamed_CFF'
	#
	if set_name==None:
		# still none? must be a string type ary_in...
		set_name=ary_in
	#
	########################
	# depricating....
	if not over_ride:
		# ary_in=None, m0=7.0, b_0 = 0.0, nyquist_factor=.5, do_plot=False, fnum=0, set_name=None, f_gt_lt=operator.lt, section_id=16, detail=False
		x=psa_forecast(ary_in=ary_in, m0=m0, b_0=b_0, nyquist_factor=nyquist_factor, do_plot=do_plot, fnum=fnum, set_name=set_name, f_gt_lt=f_gt_lt, section_id=section_id, detail=False)	# though we might want detail=True...
		# evaluate_alert_segments(alert_segments=None, CFF=None, section_id=None, m0=7.0, do_plot=True, fnum=0)
		return evaluate_alert_segments(alert_segments=x,CFF=CFF, section_id=section_id, do_plot=do_plot, fnum=fnum)

	#
	recurrence_data = mean_recurrence(ary_in=CFF, m0=m0)
	nyquist_len = max(int(nyquist_factor*recurrence_data['mean_dN_fault']), 2)
	nyquist_time = nyquist_factor*recurrence_data['mean_dT']
	#
	trend_data = get_trend_analysis(ary_in=CFF, nyquist_len = nyquist_len, nyquist_time=nyquist_time)
	#trend_data_dict = {trend_data['event_year']:x for x in trend_data}
	CFF_dict = {x['event_year']:x for x in CFF}
	print "trend lengths: ", len(trend_data), len(CFF), nyquist_len
	max_n = len(trend_data)
	#
	# trend_data includes columns: ['event_number', 'event_year', 'lin_fit_a', 'lin_fit_b', 'rb_ratio', 'interval_rate_ratios']
	#
	# first, just get the total time under alert:
	alert_time = 0.0
	alert_segments = [[]]		# a collection of lists...
	running_b_sequence = []
	#
	for i, rw in enumerate(trend_data):
		# when we find b<b_0, we issue an alert until the next earthquake -- unless this one was 'big',
		# in this case m>=7.
		#
		if i>=(max_n-1): break
		#
		# note the comment evidence of variations on this metric, primarily involving some sort of mean-slope averaging.
		# a more exhaustive, automated MC approach, analysis is necessary to be certain, but preliminary analysis suggests that
		# we don't gain anything from the averaging... in fact, we tend to loose measurable, at least on fault 16, where most
		# of the prelim examination was done.
		
		#mean_b = numpy.mean(trend_data['lin_fit_b'][max(0, i-nyquist_len) : i+1])
		#
		this_b = rw['lin_fit_b']
		#
		#this_b = mean_b
		#running_b_sequence = (running_b_sequence + [this_b])[max(0, len(running_b_sequence)-nyquist_len):]
		#mean_b = numpy.mean(running_b_sequence)
		#this_b = mean_b
		#
		#if this_b >= b_0:
		if not f_gt_lt(this_b, b_0):		# abstracting the gt/lt option...
			# null case
			# if we've been collecting "alert" events, stop. if not, just troll along...
			if len(alert_segments[-1])>0:
				#
				alert_segments+= [[]]
		
		else:
			# f_gt_lt(this_b, b_0) is True
			# accelerating case (if f_gt_lt = operator.lt):
			#this_mag = CFF[-max_n:][i+1]['event_magnitude']
			this_mag = CFF_dict[rw['event_year']]['event_magnitude']
			#
			if len(alert_segments[-1])==0:
				alert_segments[-1]+=[[trend_data[i]['event_year'], trend_data[i]['lin_fit_b']]]
			#
			if this_mag<m0:
				# add the next event as the alert (aka, issue an alert until we have new data).
				#print "len alert_seg: %d" % len(alert_segments[-1])
				#
				# generalize language a bit:
				alert_segments[-1]+=[[trend_data[i+1]['event_year'], this_b]]
				#pass
				#
			if this_mag>=m0:
				# this is "the" earthquake. add this entry (it's probably already there) from the previous entry.
				#
				alert_segments[-1]+=[[trend_data[i]['event_year'], this_b]]
			#
		#
	#
	while len(alert_segments)>0 and len(alert_segments[-1])==0: alert_segments.pop()
	#
	# now, calculate total alert time:
	total_alert_time = 0.0
	total_total_time = CFF[-1]['event_year'] - CFF[0]['event_year']
	#
	for alert_segment in alert_segments:
		total_alert_time += (alert_segment[-1][0] - alert_segment[0][0])
	#
	# and prediction success:
	n_predicted = 0
	n_missed = 0
	n_total = 0
	#
	# was an alert issued at the time of an m>m0 event?
	# we should clean this up, but for now, make and use dictionaries...
	#alert_dict = {x[0]:x[1] for x in 
	#
	# this approach is problematic and i think might be producing some artifact effects.
	#alert_dict = {}
	#for seg in alert_segments:
	#	#
	#	for rw in seg:
	#		alert_dict[rw[0]] = rw[1]
	#
	j_alert_start=0
	for i, rw in enumerate(CFF):
		if rw['event_magnitude']<m0: continue
		n_total+=1
		#
		# this is a really dumb way to do this, so fix it:
		#
		# otherwise, its a big one. we need two conditions:
		# 1) the alert was active during the event.
		# 2) the alert had already been issued.
		CFF_index = rw['event_number']
		prev_year = CFF[i-1]['event_year']
		#if alert_dict.has_key(rw['event_year']):
		#if alert_dict.has_key(rw['event_year']) and alert_dict.has_key(prev_year):
		#
		## note: this is keying off a float index, so it has potential to be unreliable.
		#if alert_dict.has_key(rw['event_year']):						
		#	n_predicted += 1
		#else:
		#	n_missed += 1
		#
		try:
			while len(alert_segments)>j_alert_start and rw['event_year']>alert_segments[j_alert_start][-1][0] and j_alert_start<(len(alert_segments)-1): j_alert_start+=1
		except:
			print "alert indexing error. ", len(alert_segments)
		#
		for alert in alert_segments[j_alert_start:]:
			if rw['event_year'] > alert[0][0] and rw['event_year'] <= alert[-1][0]:
				n_predicted +=1
				#n_missed -=1
				break

	n_missed = n_total-n_predicted
	#
	if do_plot:
		# diagnostic plots of forecast metric:
		plt.figure(fnum)
		plt.clf()
		plt.title('VC PSA Hazard Metric')
		#
		ax_mags = plt.gca()
		ax_metric = ax_mags.twinx()		
		#
		ax_ints = ax_metric.twinx()
		ax_ints.set_yscale('log')
		ax_metric.set_yscale('linear')
		#
		min_mag = min(CFF['event_magnitude'])
		max_mag = max(CFF['event_magnitude'])
		intervals = [x-CFF['event_year'][i] for i,x in enumerate(CFF['event_year'][1:])]
		#
		#metric_pad_factor = min_mag
		#metric_pad_factor = 0.
		#min_metric = 0.
		min_metric = alert_segments[0][0][1]
		max_metric = alert_segments[0][0][1]
		#
		# do a quick spin to get min/max values and other useful stats:
		for segment in alert_segments:
			X,Y = zip(*segment)
			min_metric = min(min(Y), min_metric)
			max_matric = max(max(Y), max_metric)
		#
		for segment in alert_segments:
			X,Y = zip(*segment)
			min_metric = min(min(Y), min_metric)
			#
			# ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b']  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='m', zorder=1, alpha=.5)
			#ax_trend2.fill_between([x['event_year'] for x in trend_data], [1.  for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0.0 for x in trend_data], color='m', zorder=1, alpha=.25)
			#
			# show the metric value:
			#ln_metric = ax_metric.fill_between(X,[y for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [y<0. for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			ln_metric = ax_metric.fill_between(X,[y for y in Y],y2=[0.0 for y in Y], color='m', alpha=.3, where = [f_gt_lt(y,0.) for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			ax_metric.plot([x for i,x in enumerate(X) if Y[i]<0], [0.0 for y in Y if y<0], 'm-')
			ax_metric.plot([x for i,x in enumerate(X) if Y[i]<0], [y for y in Y if y<0], 'm-')
			#
			# and just an "alert!" box:
			#ln_metric = ax_metric.fill_between(X,[min_metric for y in Y],y2=[0.0 for y in Y], color='m', alpha=.15, where = [y<0. for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			ln_metric = ax_metric.fill_between(X,[min_metric for y in Y],y2=[0.0 for y in Y], color='m', alpha=.15, where = [f_gt_lt(y,0.) for y in Y], zorder=7, label='Hazard metric: $\\eta (b)$' )
			#
			#ax_mags.fill_between(X, [min_mag for x in X], [m0 for x in X], zorder=5, alpha=.2, color='c')
		#
		ln_ints = ax_ints.plot(CFF['event_year'][1:], intervals, 'b.-', lw=2, ms=7, label='Intervals $\Delta t = t_i - t_{i-1}$')
		ax_metric.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		ln_mags = ax_mags.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='g', alpha=.7, lw=1.5, label='magnitudes')
		X_big_mags, Y_big_mags = zip(*[[x['event_year'], x['event_magnitude']] for x in CFF if x['event_magnitude']>m0])
		ax_mags.vlines(X_big_mags, [min_mag for x in Y_big_mags], Y_big_mags, color='r', alpha=.9, lw=2.75)
		#
		# cosmetics:
		ax_mags.set_xlabel('Event Year $t$')
		ax_ints.set_ylabel('Inter-event interval $\Delta t$')
		ax_metric.yaxis.set_ticks([])
		ax_mags.set_ylabel('Earthquake magnitude $m$')
		#
		ax_metric.set_ylim(ymin=1.15*min_metric, ymax = -min_metric/10.)
	#
	print "preliminary report:"
	#H = float(n_predicted)/(float(n_predicted)+n_missed)
	H=float(n_predicted)/float(n_total)
	F = total_alert_time/total_total_time
	#print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, total_alert_time/total_total_time)
	#print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, float(n_predicted)/(float(n_predicted)+n_missed))
	#print "score: H-F: %f " % (float(n_predicted)/(float(n_predicted)+n_missed) - total_alert_time/total_total_time)
	print "alert time: %f / %f :: %f " % (total_alert_time, total_total_time, F)
	print "n_predicted: %d, n_missed: %d (%f )" % (n_predicted, n_missed, H)
	print "score: H-F, H/F: %f / %f " % (H-F, H/F)

	
	#
	#return alert_segments
	# , 'alert_segments':alert_segments
	#return {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'alert_segments':alert_segments, 'ary_in_name':ary_in, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
	return {'total_alert_time': total_alert_time, 'total_time':total_total_time, 'n_predicted':n_predicted, 'n_missed':n_missed, 'alert_segments':alert_segments, 'ary_in_name':set_name, 'b':b_0, 'm0':m0, 'nyquist_factor':nyquist_factor}
#
def plot_fc_metric_1(file_profile = 'data/VC_CFF_timeseries_section_*.npy', m0=7.0, b_0=0.0, nyquist_factor=.5, do_spp=False, do_plot=False, do_clf=True, n_cpus=None):
	'''
	# depricated? this script seems to be showing its age as well...
	#
	# simple ROC diagram. see the optimized one...
	# also note that the parameter list for forecast_metric_1() has changed, so generally speaking this 
	# function needs some maintenance or to be retired.
	#
	# scatter plot of hit_ratio vs alert_time_ratio for as many data as we throw at it.
	# note, this uses a single value(s) for (b_0, nyquist_factor). see optimized versions as well.
	#
	# ... and this works ok, but it is a bit limited in exactly what sets we can optimise. see simple_mpp_optimizer()
	# and simple_metric_optimzer() above.
	'''		
	# 
	G=glob.glob(file_profile)
	#
	X, Y = [], []
	#
	if do_spp:
		# SPP version (if we want to force this for some reason):
		for g in G:
			fc_data = forecast_metric_1(ary_in=g, m0=m0, b_0 = b_0, nyquist_factor=nyquist_factor)
			#
			Y+=[float(fc_data['n_predicted'])/(fc_data['n_predicted'] + float(fc_data['n_missed']))]
			X+=[fc_data['total_alert_time']/fc_data['total_time']]
	#
	else:
		# use MPP:
		if n_cpus==None: n_cpus = mpp.cpu_count()
		pool = mpp.Pool(n_cpus)
		if (m0, b_0, nyquist_factor)==forecast_metric_1.__defaults__[1:4]:
			print "defaults. use map_async()"
			# i'm guessing this is faster...
			result_set = pool.map_async(forecast_metric_1, G)
			pool.close()
			pool.join()
			#
			#print result_set
			resultses = result_set.get()	# will be a list of dictionary objects
			
		else:
			# add/"apply" each file to the pool.
			print "not default. use apply_async()"
			result_set_list = [pool.apply_async(forecast_metric_1, args = (g, m0, b_0, nyquist_factor)) for g in G]
			pool.close()
			pool.join()
			#
			#print "these prams: ", b_0, nyquist_factor
			#return result_set_list
			resultses = [x.get() for x in result_set_list]	# each entry in result_set_list is a dict; now we have a list of dicts.
			
		#
		# still MPP...
		#
		XY = [[x['total_alert_time']/x['total_time'] ,float(x['n_predicted'])/(x['n_predicted'] + float(x['n_missed']))] for x in resultses]
		X,Y = zip(*XY)
	# end MPP. now, either by SPP or MPP, we have X,Y
	#
	#
	#return fc_datas
	#
	if do_plot:
		#
		plt.figure(0)
		if do_clf: plt.clf()
		plt.plot(X, Y, '.')
		plt.plot([0., 1.], [0., 1.], '-')
		plt.ylabel('n_predicted/n_total')
		plt.xlabel('percent alert time')
	#
	#if resultses[0].has_key('alert_segments'): [x.pop('alert_segments') for x in resultses]
	resultses[0] = {key:val for key,val in resultses[0].iteritems() if key!='alert_segments'}	# a little bit faster and more elegant...
	return resultses

#
def plot_CFF_ary(ary_in=None, fnum=0, nyquist_factor=.5, gt_lt_eval=operator.lt, section_id=16):
	'''
	# (depricate? this bit can probably be removed. it was originally intended to investigate the usefulness of plotting CFF
	# time series... which appear to be mostly useless. subsequent figure-scripts plotting  interval \Delta t and d(\Delta t)/dt
	# apper more useful, and the alert system appears to be a little bit more complicated with nyquist_factor and b_0, so this script might 
	# be more clutter than useful at this stage.
	#
	# this script is for some of the earlier CFF numpy array types. newer arrays may require different scripting. (but it seems to work ok).
	# basics: plots the triple time-series figure with [intervals, magnitudes, alert], [CFF, Delta CFF (?)], and [mags,slopes] ].
	# this was a development figure where be basically decided not to use CFF. and the alert figure got moved to ???.
	# for older data sets (most likely):
	#	# 2 cols: event_number, CFF_initial
	#   # 3 cols: event_number, event_year, CFF_initial
	#   # 4 cols: event_number, event_year, CFF_initial, CFF_final
	#  later versions will include column names.
	#
	# ary_in: can be a numpy.recarray or a string/filename like: 'data/VC_CFF_timeseries_section_125.npy'
	#   numpy array would be like numpy.load(ary_in)
	#
	# gt_lt_eval: how do we evaluate greter/less than? submit gt/lt functions: operator.{lt, le, gt, ge} or
	# strings: {'le', 'leq', 'ge', 'geq'}, and we trap a few others as well.
	#
	# (this is primarily a development script and may need to be depricated).
	#
	'''
	#
	# are we going to be a greater to less than camp?
	#gl_lt_eval = operator.lt		# "less than", "less that or equal", "greater than", "greater than or equal"
	#gl_lt_eval = operator.le
	#gl_lt_eval = operator.gt
	#gt_lt_eval = operator.ge		# greater than or equal to; ge(4,5)=False, ge(4,4)=True, ge(4,3)=True
	# and allow for string inputs:
	if isinstance(gt_lt_eval, str):
		gt_lt_eval = gt_lt_eval.lower()
		if gt_lt_eval in ('lt', 'less', 'lessthan', 'less_than'): gl_lt_eval = operator.lt
		if gt_lt_eval in ('lte', 'leq', 'lessthanorequal', 'less_than_or_equal'): gl_lt_eval = operator.le
		if gt_lt_eval in ('ge', 'greater', 'greaterthan', 'greater_than'): gl_lt_eval = operator.gt
		if gt_lt_eval in ('gte', 'geq', 'greaterthanorequal', 'greater_than_or_equal'): gl_lt_eval = operator.ge
	#
	if section_id!=None and (ary_in==None or len(ary_in)==0):
		# guess the ary_in file name from section_id (and later on, trap for multiple section_id values):
		ary_in = 'data/VC_CFF_timeseries_section_%d.npy' % section_id
		
	if isinstance(ary_in, str):
		CFF = numpy.load(ary_in)
	else:
		CFF = ary_in
	#
	recurrence_data = mean_recurrence(ary_in=CFF, m0=7.0)
	nyquist_len = int(nyquist_factor*recurrence_data['mean_dN_fault'])
	nyquist_time = nyquist_factor*recurrence_data['mean_dT']
	#
	trend_data = get_trend_analysis(ary_in=CFF, nyquist_len = nyquist_len, nyquist_time=nyquist_time)
	#
	# what kind of array did we get?
	# this is not very efficient, in that we rewrite the whole enchilada for the new types, but i don't
	# expect that we'll be returning to these unformatted (unstructured) array types.
	# so, focus on newer structured arrays. leave the old thing in out of principle.
	#	
	if isinstance(CFF, numpy.recarray)==True:
		# it's a structured array.
		#cols = map(operator.itemgetter(0), CFF.dtype.descr)
		col = CFF.dtype.names
		# cols should be like: ['event_number', 'event_year', 'event_magnitude', 'cff_initial', 'cff_final', 'event_area']
		#
		f=plt.figure(fnum)
		f.clf()
		#
		# create two axes:
		# magnitudes plot
		# bottom:
		ax_mag = f.add_axes([.1, .05, .85, .25])
		#
		# CFF plot. (center):
		ax_CFF = f.add_axes([.1, .35, .85, .25], sharex=ax_mag)
		ax_dCFF = ax_CFF.twinx()	# over-plot stress (CFF) drop...
		#
		# top (intervals, etc.)
		ax_ints  = f.add_axes([.1, .65, .85, .25], sharex=ax_mag)
		ax_mag2  = ax_ints.twinx()
		ax_trend2 = ax_ints.twinx()
		#
		ax_trend = ax_mag.twinx()
		#
		X_init = CFF['event_year']
		X_finals = [x+.01 for x in X_init]	# this to show the stress drops after a mainshock (X_final is the time of the data point after the event).
		#
		Y0 = -1.*CFF['cff_initial']
		Y_final = -1.*(CFF['cff_final'])
		
		X = list(X_init) + list(X_finals)
		X.sort()
		#
		intervals = X_init[1:] - X_init[:-1]
		
		big_mags = zip(*[[x['event_year'], x['event_magnitude']] for x in CFF if x['event_magnitude']>7.0])
		#
		Y = []
		CFF_drops = [(x['cff_initial'] - x['cff_final'])**2. for x in CFF]
		for i, y in enumerate(Y0):
			Y += [Y0[i]]
			Y += [Y_final[i]]
		#
		# use "peak" values to cut through some noise.
		peaks = get_peaks(zip(*[X,Y]), col=1, peak_type='upper')
		X_peaks, Y_peaks = zip(*peaks)
		#
		# CFF Plot:
		#ax = plt.gca()
		ax_CFF.set_xscale('linear')
		ax_CFF.set_yscale('log')
		ax_CFF.set_ylabel('CFF')
		ax_dCFF.set_xscale('linear')
		ax_dCFF.set_yscale('log')

		# first, raw CFF (initial):
		ax_CFF.plot(X, Y, '.-', color='b', alpha=.2, zorder=4)
		ax_CFF.fill_between(X, Y, y2=min(Y), color='b', alpha=.2, zorder=4)
		ax_CFF.plot(X_peaks, Y_peaks, '-', zorder=5)
		ax_dCFF.plot(CFF['event_year'], CFF_drops, 'g.-', zorder=7, alpha=.9)
		#
		# Magnitude plot (by itself):
		ax_mag.set_xscale('linear')
		ax_mag.set_yscale('linear')
		ax_mag.set_ylabel('event magnitude $m$')
		ax_mag.set_xlabel('event year $t$')
		min_mag = min(CFF['event_magnitude']) - .5
		ax_mag.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='b', alpha=.9)
		#
		ax_trend.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], 'r-', zorder=5, alpha=.8)
		#ax_trend.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], y2=[0.0 for x in trend_data], where=[x['lin_fit_b']<0. for x in trend_data], color='r', zorder=2, alpha=.8)
		ax_trend.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], y2=[0.0 for x in trend_data], where=[gt_lt_eval(x['lin_fit_b'],0.) for x in trend_data], color='r', zorder=2, alpha=.8)
		ax_trend.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		ax_trend.set_ylabel('(log) interval slope $b$')
		#
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] for x in trend_data], 'r-', zorder=5, alpha=.8)
		#
		ax_trend2.fill_between([x['event_year'] for x in trend_data], [x['lin_fit_b']  for x in trend_data], y2=[0.0 for x in trend_data], where=[gt_lt_eval(x['lin_fit_b'],0.) for x in trend_data], color='m', zorder=1, alpha=.5)
		ax_trend2.fill_between([x['event_year'] for x in trend_data], [1.  for x in trend_data], y2=[0.0 for x in trend_data], where=[gt_lt_eval(x['lin_fit_b'],0.) for x in trend_data], color='m', zorder=1, alpha=.25)
		#		
		ax_trend2.plot([trend_data['event_year'][0], trend_data['event_year'][-1]], [0., 0.], 'k--')
		#
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['rb_ratio'] for x in trend_data], 'c--')
		#ax_trend2.plot([x['event_year'] for x in trend_data], [x['lin_fit_b'] + x['rb_ratio'] for x in trend_data], 'c-')
		
		#
		# Interval Plots:
		ax_ints.set_xscale('linear')
		ax_ints.set_yscale('log')
		ax_ints.set_ylabel('intervals $\\Delta t$')
		ax_ints.plot(X_init[1:], intervals, '.-', alpha=.9)
		ax_mag2.vlines(CFF['event_year'], [min_mag for x in CFF['event_magnitude']], CFF['event_magnitude'], color='g', alpha=.9, lw=2)
		ax_mag2.vlines(big_mags[0], [min_mag for x in big_mags[1]], big_mags[1], color='r', lw=2.5, alpha=.9)
		ax_mag2.set_ylabel('magnitude $m$')
	#	
	if isinstance(CFF, numpy.recarray)==False:
		# a regular, old-style, numpy.ndarray -- aka, no columns. guess the column structure from what we know...
		#
		# assume either [event_id/num, year, CFF] or [event_id/num, CFF]
		if len(CFF[0])==2:
			y_col=1
		if len(CFF[0])>=3:
			x_col=1
			y_col=2
		if len(CFF[0])>=5:
			x_col=1
			y_col=3
		#CFF_peaks = get_peaks(data_in=CFF, col=y_col, peak_type='lower')
		#
		zCFF = zip(*CFF)
		X=zCFF[x_col]
		Y=zCFF[y_col]
		Y=[-1*y for y in Y]
		CFF_peaks = get_peaks(data_in = zip(*[X,Y]), col=1, peak_type='upper')
		#
		zCFF_peaks = zip(*CFF_peaks)
		X_peaks = zCFF_peaks[0]
		Y_peaks = zCFF_peaks[1]
		#
		'''
		if len(CFF[0])==2:
			X=zCFF[0]
			Y=zCFF[1]
		if len(CFF[0])>=3:
			X=zCFF[1]
			Y=zCFF[2]
		'''
		#
		plt.figure(fnum)
		plt.clf()
		#
		ax = plt.gca()
		ax.set_xscale('linear')
		ax.set_yscale('log')
		plt.fill_between(X, Y, y2=min(Y), color='b', alpha=.2, zorder=4)
		plt.plot(X_peaks, Y_peaks, '.-', zorder=5)				
	#
	return CFF
	
