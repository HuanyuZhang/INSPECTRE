# The proposed estimator
# t = the horizon
# n = number of samples
# r = the parameter for the estimator
# maxfreq = maximum frequency of the underlying sample
# attentuation = 1 always
# Code from Orlitsky, Suresh, and Wu.

import numpy as np
import prob_coefficients as coeff
import math as mt

def proposed_estimator(t,n,r,maxfreq,attentuation):
    if (t > 0 and n > 0 ):
        poissontail = coeff.poissontail
        log = mt.log
        exp = mt.exp
        vectail = coeff.poissontailvecscipy(r,maxfreq+5)
        estimator = [0] * (maxfreq+1)
        logt = log(t)
        # attentuation = attentuation_poisson(r,t)
        for i in range(1,maxfreq+1):
            if vectail[i] > 0:
                estimator[i-1] = ((-1)**(i+1))*exp((i*logt) + log(vectail[i]))*attentuation
    else:
        estimator = [0] * (maxfreq+1)
    return estimator
