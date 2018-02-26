# iter = number of iterations
# maxindex = maximum value of t
# minindex = minimum value of t
# random seed = 222

# Reads a file app_c_total.csv" file, which is the histogram of the data and generates a plot on how well
# it estimates and sampling is done without replacement.

import numpy as np
import matplotlib.pyplot as pt
import math as mt
import csv

from generatesamples import generatesamples
from linear_estimator import linear_estimator
from proposed_estimator import proposed_estimator as prop
from fgt_estimator import fgt_estimator as fgt
from collections import Counter
from linear_estimator_Laplace import linear_estimator_Laplace

file_name = 'lastnames_total.csv'
eps_index = [2,1,0.5]

np.random.seed(222)
log = mt.log
maxindex = 10
minindex = 1
iter = 100

pt.figure()
figname = 'lastnames_total_eps.pdf'

for eps in eps_index:
    f = open(file_name, 'rt')
    xlab = 'Fraction of seen names'
    ylab = 'RMSE'

    ftemp = csv.reader(f)
    totalfreq = [int(y[0]) for y in ftemp]

    k = len(totalfreq) # number of species
    ntotal = sum(totalfreq) # number of samples

    index = np.arange(minindex,maxindex,1)
    x = [float(y)/maxindex for y in index]
    z = [ntotal*(1-y) for y in x] 

    Laplacemse = [0 for y in index]
    SGTmse = [0 for y in index]

    expectedvalue = [0 for y in index]

    for it in range(0,iter):
        for j in range(0,maxindex-minindex):
            n = int(x[j]*ntotal)
            t = float(ntotal-n)/n
            freq = [0] * k
            seensofar = 0
            freqsofar = 0

            #Sample without replacement
            for y in range(0,k):
                if(n > seensofar):
                    freq[y] += np.random.hypergeometric(totalfreq[y], ntotal - freqsofar-totalfreq[y], n - seensofar)
                seensofar += freq[y]
                freqsofar += totalfreq[y]
            maxfreq = max(freq)
            totalseen = sum(min(1,y) for y in freq)
            prevalence = Counter(freq)
            expectedvalue[j] += k 
                    
            # Parameters for SGT
            sc = 1.0
            tnew = max(t,1.001)
            propr = float(sc*log(n*(tnew+1)*(tnew+1)/(tnew-1)))/(float(2*tnew))

            # Estimator: when t>1 use SGT or use Good Turing estimator
            if(t > 1):
                temp = totalseen + linear_estimator(prevalence, prop(t,n,propr,maxfreq,1),n,t,maxfreq)
            else:
                temp = totalseen + linear_estimator(prevalence, fgt(t,n),n,t,maxfreq)
            SGTmse[j]+= (temp - k)**2
            
            # Estimator: our private estimator
            if(t > 1):
                temp =  linear_estimator_Laplace(prevalence, prop(t,n,propr,maxfreq,1), n, t,maxfreq,propr,eps)
            else:
                add_noise = np.random.laplace(0,max(t,0.01)/eps,1)               
                temp =  linear_estimator(prevalence, fgt(t,n),n,t,maxfreq)+np.asscalar(add_noise)
            temp = max(min(temp,n*t),0)  
            Laplacemse[j]+=((temp + totalseen- k)**2)


    biasSGTmse = np.array( [np.sqrt(float(y/iter)) for (yt,y) in zip(z,SGTmse)])
    biasLaplacemse = np.array( [np.sqrt(float(y/iter)) for (yt,y) in zip(z,Laplacemse)])

    temp_label = 'Private '+'eps='+str(eps)
    if(eps==eps_index[0]):           
        pt.plot(x,biasSGTmse, label = 'Non-private',linewidth = 2)
    pt.plot(x,biasLaplacemse, label = temp_label,linewidth = 2)

    pt.xlim([float(1)/maxindex,1])
    pt.xlabel(xlab,fontsize = 10)
    pt.ylabel(ylab, fontsize = 10)
    pt.gcf().subplots_adjust(left=0.17)
    pt.legend(bbox_to_anchor=(1.0,1.0))

    pt.savefig('%s' %(figname), format='pdf')
