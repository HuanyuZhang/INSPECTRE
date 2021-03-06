#!/usr/bin/env python3
"""
Main library for entropy estimation
Based largely off the implementation of Wu and Yang: https://github.com/Albuso0/entropy/blob/master/entropy.py
"""

import linecache
from math import log, floor
import numpy as np



class Entropy():
    """
    Entropy estimator
    """
    def __init__(self, k, L=None, M=None, N=None):
        """
        Args:
        k: int, required
        alphabet size

        L: int
        Polynoimal degree. Default c0*log(k)

        M: float
        M/n is the right-end of approximation interval. Default c1*log(k)

        N: int
        Threshold to apply polynomial estimator. Default c2*log(k)
        """
        self.k = k
        self.degree = L if L != None else floor(1.6*log(k)) #c0
        self.ratio = M if M != None else 3.5*log(k) #c1
        self.n_threhold = N if N != None else floor(1.6*log(k)) #c2


    def estimate(self, fin):
        """
        Polynomial estimator from a given fingerprint

        Args:
        fin: list of tuples (frequency, count)
        fingerprint of samples, fin[i] is the number of symbols that appeared exactly i times

        Return:
        the estimated entropy (bits)
        the sensitivity of the estimator (bits)
        """
        # get total sample size
        num = get_sample_size(fin)
        
        # get linear estimator coefficients
        a_coeffs = read_coeffs(self.degree)
        g_coeffs = np.empty(self.n_threhold+1)
        for j in range(self.n_threhold+1):
            start = self.degree if j > self.degree else j
            g_coeffs[j] = a_coeffs[start]
            for i in range(start, 0, -1):
                g_coeffs[j] = a_coeffs[i-1] + g_coeffs[j] * (j-i+1) / self.ratio
            g_coeffs[j] = (g_coeffs[j] * self.ratio + log(num/self.ratio) * j)/num

        # get estimate
        h_estimate = 0
        sym_num = 0

        sensitivity=[]
        for freq, cnt in fin:
            sym_num += cnt
            if freq > self.n_threhold: # plug-in
                p_hat = 1.*freq/num
                h_estimate += (-p_hat*log(p_hat)+0.5/num)*cnt
            else: # polynomial estimator
                h_estimate += g_coeffs[freq]*cnt
                # when we change one sample from this species, the entropy changes at most
                sensitivity.append(2*np.abs(g_coeffs[freq]-g_coeffs[min(freq+1,self.n_threhold)]))


        h_estimate += 1.0*g_coeffs[0] * (self.k-sym_num)
        h_estimate = h_estimate if h_estimate > 0 else 0
        # the sensitivity for the estimator is the max entropy change when we change only one sample
        return (h_estimate/log(2),max(sensitivity)/log(2))

    def estimate_Miller_Madow(self, fin):
        """
        Miller Madow estimator from a given fingerprint

        Args:
        fin: list of tuples (frequency, count)
        fingerprint of samples, fin[i] is the number of symbols that appeared exactly i times

        Return:
        the estimated entropy (bits)
        """
        num = get_sample_size(fin)

        h_estimate = 0
        sym_num = 0
        for freq, cnt in fin:
            if freq > 0:
                sym_num += cnt
                p_hat = 1.*freq/num
                h_estimate += (-p_hat*log(p_hat))*cnt
        return (h_estimate+1.0*(sym_num-1)/num/2)/log(2)

    def estimate_plug(self, fin):
        """
        Miller Madow estimator from a given fingerprint

        Args:
        fin: list of tuples (frequency, count)
        fingerprint of samples, fin[i] is the number of symbols that appeared exactly i times

        Return:
        the estimated entropy (bits)
        """
        num = get_sample_size(fin)

        h_estimate = 0
        sym_num = 0
        for freq, cnt in fin:
            sym_num += cnt
            if freq > 0:
                p_hat = 1.*freq/num
                h_estimate += (-p_hat*log(p_hat))*cnt

        return h_estimate/log(2)

def get_sample_size(fin):
    """
    get total sample size from a given fingerprint
    """
    num = 0
    for freq, cnt in fin:
        num += freq * cnt
    return num

def read_coeffs(degree):
    """
    read the coefficients from a given degree
    """
    line = linecache.getline("coeffs.txt", degree)
    return np.asarray(line.split()[1:]).astype(np.float)


def sample_to_fin(samples):
    """
    Return a fingerprint from samples
    """
    return hist_to_fin(sample_to_hist(samples))


def hist_to_fin(hist):
    """
    Return a fingerprint from histogram
    """
    fin = Counter()
    for freq in hist:
        fin[freq] += 1
    return fin.items()

def sample_to_hist(samples):
    """
    Return a histogram of samples
    """
    freq = Counter()
    for symbol in samples:
        freq[symbol] += 1
    return np.asarray(list(freq.values()))

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


class Counter(dict):
    """
    Class for counting items
    """
    def __missing__(self, key):
        return 0
