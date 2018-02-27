# Generates synthetic data and creates plots
# from 6 distributions.
# n = number of seen samples
# k = support sie
# iter = number of iterations
# maxindex = max value for t
import numpy as np
import matplotlib.pyplot as pt
import math as mt

from generatesamples import generatesamples
from linear_estimator import linear_estimator
from proposed_estimator import proposed_estimator as prop
from unseen_expected import unseen_expected as expected
from collections import Counter
from time import time

log = mt.log
k = 100
n = k//2
maxindex = 10
xdimsubplots = 2
ydimsubplots = 3
iter = 100

fig, ax = pt.subplots(xdimsubplots, ydimsubplots)
figname = 'synthetic_plots2_total_eps.pdf'

eps_index = [10,2,1]
eps_length = len(eps_index)

index = np.arange(1,maxindex+1,1)
x = [float(y) for y in index]
z = [y*n for y in x]

# generate distributions
distarray = []
subtitle = []
raw_distribution = [1] * k
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Uniform') 
  
raw_distribution = [1] * (k//2) + [3] * (k//2)
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Two steps')


raw_distribution = [1/(float(i)**(0.5)) for i in range(1,k)]
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Zipf 1/2')

raw_distribution = [1/(float(i)**(1)) for i in range(1,k)]
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Zipf 1')
 
raw_distribution = [0] * k
for i in range(0,k):
    raw_distribution[i] = np.random.gamma(1,1)
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Dirichlet-1 prior')

raw_distribution = [0] * k
for i in range(0,k):
    raw_distribution[i] = np.random.gamma(0.5,1)
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Dirichlet-1/2 prior')

for firstindex in range(0,xdimsubplots):
    for secondindex in range(0,ydimsubplots):

        Laplacemse = [ [0 for y in index] for q_add in eps_index]
        SGTmse = [0 for y in index]

        distribution = distarray[firstindex*ydimsubplots + secondindex]
        for it in range(0,iter):
            # generate and count samples
            freq = generatesamples(distribution,n)
            maxfreq = max(freq)
            prevalence = Counter(freq)
            for j in range(0,maxindex):
                t = x[j]

                #expected number of new species in next n*t examples
                currexpected= expected(n,t,freq, distribution)
                
                # Parameters for estimators
                sc = 1.0
                tnew = max(t,1.001) # make sure t>1
                propr = float(sc*log(n*(tnew+1)*(tnew+1)/(tnew-1)))/(float(2*tnew))
                tnew = t
                atten = 1
                # SGT estimator
                temp = linear_estimator(prevalence, prop(t,n,propr,maxfreq,atten),n,t,maxfreq)                
                SGTmse[j]+= (temp - currexpected)**2

                for i in range(eps_length): #loop for different eps
                    # private estimator the sensitivity is bounded as 2.0 * (1+np.exp(r*(t-1)))
                    temp2 = temp+np.asscalar(np.random.laplace(0,(2.0 * (1+np.exp(propr*(t-1))))/eps_index[i],1))
                    # The estimator must be between 0 and n*t
                    temp2 = max(min(temp2,n*x[j]),0) 
                    Laplacemse[i][j]+= (temp2 - currexpected)**2    

        biasSGTmse = np.array( [np.sqrt(float(y/iter))/(n*yt) for (yt,y) in zip(x,SGTmse)])
        
        #plot rmse for non-private estimators
        ax[firstindex,secondindex].plot(x,biasSGTmse, label = 'Non-private',linewidth = 2)
        #plot rmse for private estimators
        for i in range(eps_length):   
            biasLaplacemse = np.array( [np.sqrt(float(y/iter))/(n*yt) for (yt,y) in zip(x,Laplacemse[i])])
            temp_label = 'Private '+'eps='+str(eps_index[i])
            ax[firstindex,secondindex].plot(x,biasLaplacemse, label = temp_label, linewidth = 2)

        # set parameter for plot
        ax[firstindex,secondindex].set_xlim([1,maxindex])
        ax[firstindex,secondindex].set_ylim(ymin=0)
        ax[firstindex,secondindex].set_title('%s' %(subtitle[firstindex*ydimsubplots+secondindex]), fontsize = 15)
        if (secondindex == 0):
            ax[firstindex,secondindex].set_ylabel('$Normalized-MSE$', fontsize = 17)
        if (firstindex == xdimsubplots-1):
            ax[firstindex,secondindex].set_xlabel('$t$', fontsize =17)

ax[xdimsubplots-1,ydimsubplots-1].legend(loc=4)
pt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
pt.setp([a.get_yticklabels() for a in ax[:, ydimsubplots-1]], visible=False)
pt.setp([a.get_yticklabels() for a in ax[:, ydimsubplots-2]], visible=False)
fig.set_size_inches(18.5, 10.5, forward=True)
fig.savefig('%s' %(figname), format='pdf')