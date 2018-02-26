# Given a private linear estimator "estimator" with n samples in the
# prevalence form "prevalence" with maximum frequenccy maxfreq returns
# the estimated number of unseen species in n*t samples.
# Based off non-private estimator code from Orlitsky, Suresh, and Wu.

import numpy as np
from collections import Counter

def linear_estimator_Laplace(prevalence,estimator,n,t,maxfreq,r,eps):
	output = 0
	for i in range(0,maxfreq):
	    output += prevalence[i+1]*estimator[i]
	temp = 2.0 * (1+np.exp(r*(t-1)))
	a = np.random.laplace(0,temp/eps,1)
	return output + np.asscalar(a) # float(min(max(output,0),n*t))
