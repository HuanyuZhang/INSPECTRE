# Returns the expected number of unseen symbols 
# for a distribution distribution
# with n samples and horizon t
# freq is the input frequency
# distribution is the underlying distribution.

import numpy as np
import math as mt
import generatesamples

def unseen_expected(n,t,freq,distribution):
    chance = [1 - (1-x)**(n*t) for x in distribution]
    return sum([a*(max(1-b,0)) for a,b in zip(chance,freq)])
    
