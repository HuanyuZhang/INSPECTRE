# Generates n samples from the "distribution"

import numpy as np


def generatesamples(distribution, n):
    values = np.random.rand(n)
    values.sort()
    cumulative = np.cumsum(distribution)
    z = 0
    k = len(distribution)
    freq = [0] * k
    for x in values:
        while x > cumulative[z]:
            z+=1
        freq[z]+=1         
    return freq

if __name__ == '__main__':
    generatesamples([0.2, 0.2, 0.2, 0.2, 0.2],3)
