# Given a linear estimator "estimator" with n samples in the
# prevalence form "prevalence" with maximum frequenccy maxfreq returns
# the estimated number of unseen species in n*t samples.

import numpy as np
from collections import Counter

def linear_estimator(prevalence,estimator,n,t,maxfreq):
    output = 0
    for i in range(0,maxfreq):
        output += prevalence[i+1]*estimator[i]
    # print(output)
    return output # float(min(max(output,0),n*t))
