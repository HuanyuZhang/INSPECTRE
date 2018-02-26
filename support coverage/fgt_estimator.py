# Implements FGT estimator
# Based off code of Orlitsky, Suresh, and Wu.

import numpy as np


def fgt_estimator(t,n):
    estimator = [t]
    for i in range(0,n-1):
        estimator.append(estimator[-1]*(-t))
    return estimator

