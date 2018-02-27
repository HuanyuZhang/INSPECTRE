#!/usr/bin/env python3
"""
Main file for running entropy estimators
"""

from entropy import *
import math
import matplotlib.pyplot as pt

k=1000 # alphabet size
eps=0.5 # privacy parameter

iter = 100 # iteration time
maxindex = 10

xdimsubplots = 2
ydimsubplots = 3

# parameter for Laplace polynomial estimator
l_degree = math.ceil(1.2*log(k))
M_degree = (2.0)*log(k)
N_degree = math.floor(1.6*log(k))


## generate distributions
distarray = []
subtitle = []
entropy_value = []

raw_distribution = [1] * k
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Uniform') 
temp = [-y * log(y)/log(2.0) for y in distribution]
entropy_value.append(sum(temp))
  
raw_distribution = [1] * int(k/2) + [3] * int(k/2)
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Two steps')
temp = [-y * log(y)/log(2.0) for y in distribution]
entropy_value.append(sum(temp))


raw_distribution = [1/(float(i)**(0.5)) for i in range(1,k)]
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Zipf 1/2')
temp = [-y * log(y)/log(2.0) for y in distribution]
entropy_value.append(sum(temp))

raw_distribution = [1/(float(i)**(1)) for i in range(1,k)]
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Zipf 1')
temp = [-y * log(y)/log(2.0) for y in distribution]
entropy_value.append(sum(temp))

raw_distribution = [0] * k
for i in range(0,k):
    raw_distribution[i] = np.random.gamma(1,1)
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Dirichlet-1 prior')
temp = [-y * log(y)/log(2.0) for y in distribution]
entropy_value.append(sum(temp))

raw_distribution = [0] * k
for i in range(0,k):
    raw_distribution[i] = np.random.gamma(0.5,1)
sum_raw = sum(raw_distribution)
distribution = [float(y)/float(sum_raw) for y in raw_distribution]
distarray.append(distribution)
subtitle.append('Dirichlet-1/2 prior')
temp = [-y * log(y)/log(2.0) for y in distribution]
entropy_value.append(sum(temp))

entropy = Entropy(k,L=l_degree,M=M_degree,N=N_degree)
entropy_laplace= Entropy(k)

fig, ax = pt.subplots(xdimsubplots, ydimsubplots)
figname = 'entropy_'+'eps'+str(eps)+'k'+str(k)+'.pdf'

for firstindex in range(0,xdimsubplots):
    for secondindex in range(0,ydimsubplots):

        index = np.arange(1,maxindex+1,1)
        z = [y*(int)(k/10) for y in index]

        plug_mse = [0 for y in index]
        Miller_Madow_mse = [0 for y in index]
        poly_mse = [0 for y in index]
        plug_Laplace_mse = [0 for y in index]
        poly_Laplace_mse= [0 for y in index]
        real_value = [entropy_value[firstindex*ydimsubplots+secondindex] for y in index]

        for it in range(0,iter):
            for i in range(0,maxindex):
                n = z[i]
                freq = generatesamples(distarray[firstindex*ydimsubplots+secondindex], n)
                fin = hist_to_fin(freq)

                temp1 = entropy.estimate_plug(fin)
                plug_mse[i]+= (temp1 - real_value[i])**2 

                noise = np.random.laplace(0,2.0/log(2)*log(n)/(eps*n),1)
                plug_Laplace_mse [i] += (temp1 - real_value[i]+np.asscalar(noise))**2 

                Miller_Madow_mse[i] +=  (entropy.estimate_Miller_Madow(fin)-real_value[i])**2

                (temp2,sensitivity) = entropy_laplace.estimate(fin)
                poly_mse [i] += (temp2 - real_value[i])**2 

                (temp2,sensitivity) = entropy.estimate(fin)
                noise = np.random.laplace(0,sensitivity/(eps) ,1)

                poly_Laplace_mse [i] += (temp2 - real_value[i]+np.asscalar(noise))**2 

        plug_mse = [math.sqrt (float(y)/iter) for y in plug_mse]
        Miller_Madow_mse = [math.sqrt (float(y)/iter )for y in Miller_Madow_mse]
        poly_mse = [math.sqrt (float(y)/iter)for y in poly_mse]
        poly_Laplace_mse = [math.sqrt (float(y)/iter) for y in poly_Laplace_mse]
        plug_Laplace_mse = [math.sqrt (float(y)/iter) for y in plug_Laplace_mse]

        ax[firstindex,secondindex].plot(z,plug_mse, color = 'red', label = 'Plug-in', linewidth = 2)
        ax[firstindex,secondindex].plot(z,Miller_Madow_mse, color = 'blue', label = 'Miller',linewidth = 2)
        ax[firstindex,secondindex].plot(z,poly_mse, color = 'green', label = 'Poly',linewidth = 2)
        ax[firstindex,secondindex].plot(z,poly_Laplace_mse, color = 'black', label = 'Poly-Laplace',linewidth = 2)
        ax[firstindex,secondindex].plot(z,plug_Laplace_mse, color = 'orange', label = 'Plug-in-Laplace',linewidth = 2)



        ax[firstindex,secondindex].set_title('%s' %(subtitle[firstindex*ydimsubplots+secondindex]), fontsize = 15)
        if (secondindex == 0):
            ax[firstindex,secondindex].set_ylabel('RMSE', fontsize = 17)
        if (firstindex == xdimsubplots-1):
            ax[firstindex,secondindex].set_xlabel('Number of samples', fontsize =17)

ax[xdimsubplots-1,ydimsubplots-1].legend(bbox_to_anchor=(0.75,1.00))

pt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
pt.setp([a.get_yticklabels() for a in ax[:, ydimsubplots-1]], visible=False)
pt.setp([a.get_yticklabels() for a in ax[:, ydimsubplots-2]], visible=False)
fig.set_size_inches(18.5, 12, forward=True)
fig.savefig('%s' %(figname), format='pdf')

