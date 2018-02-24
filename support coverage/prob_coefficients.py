import numpy as np
import math as mt
import scipy.stats

def poissoncoefficient(lamba, n):
    log = mt.log
    exp = mt.exp
    output = sum((log(lamba) - log(i+1)) for i in range(0,n))
    return exp(output -lamba)

def binomialcoefficient(n,p,i):
    log = mt.log
    exp = mt.exp
    output = n*log(1-p)
    for i in range(0,i):
        output += log(n-i) + log(p/(1-p)) - log(i+1)
    return exp(output)

def poissontail(lamba,n):
    return sum(poissoncoefficient(lamba,i) for i in range(n,10+n+int(2*lamba)))

def poissontailvec(lamba,n):
    log = mt.log
    exp = mt.exp
    output = [0] * (n+2)
    running = exp(-lamba)
    for i in range(0,n):
        output[i] = running
        running  = running*(float(lamba)/float(i+1))
    output[n+1] = poissoncoefficient(lamba,n+1)
    for i in range(1,n+2):
        output[n+1-i] += output[n+2-i]
    return output

def poissontailvecscipy(lamba,n):
    output = [0] * (n+2)
    for i in range(0,n):
        output[i] = 1 - scipy.stats.poisson.cdf(i-1,lamba)
    return output

def binomialtail(n,p,i):
    output = 0
    for j in range(i,n+1):
        output += binomialcoefficient(n,p,j)
    return output

def logfactorial(n):
    log = mt.log
    if n==0:
        return 0
    else:
        return sum(log(i+1) for i in range(0,n))
