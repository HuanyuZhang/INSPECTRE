# Implements FGT estimator

import numpy as np
#from mpmath import *


def fgt_estimator(t,n):
    estimator = [t]
    for i in range(0,n-1):
        estimator.append(estimator[-1]*(-t))
    return estimator

