# INSPECTRE

We provide a Python implementation of our differentially private estimators of entropy and support coverage. This is a tutorial for anyone who wants to use or build upon these estimators.

Table of contents
=================
* [Table of contents](#table-of-contents)
* [Prerequisites](#prerequisites)
* [Entropy estimators](#entropy-estimators)
* [Support coverage estimator](#support-coverage-estimator)
* [References](#reference)

Prerequisites
=====
This project is implemented in Python3, using [Numpy](http://www.numpy.org) and [matplotlib](https://matplotlib.org/index.html). Before running the code, make sure Python3, Numpy and Matplotlib are installed. Note that our code for entropy estimation is incompatible with Python2.

* [Instruction for installing Python3](https://docs.python.org/3/using/index.html)
* [Instruction for installing Numpy](https://www.scipy.org/install.html)
* [Instruction for installing Matplotlib](https://matplotlib.org/users/installing.html) 


Entropy estimators
================
In this project, we implement our estimators for private estimation of entropy, including a private version of the plug-in estimator, and a privatized estimator based on the method of best-polynomial approximation.
We compare performance between these estimators and a number of non-private estimators, including the plug-in estimator, the Miller-Madow Estimator, and the sample optimal polynomial approximation estimator from the paper [Minimax Rates of Entropy Estimation on Large Alphabets via Best Polynomial Approximation](http://ieeexplore.ieee.org/abstract/document/7444171/). 

Some of our code (including most of ```entropy.py``` and ```coeffs.txt```) is based off of the non-private estimators from the project [entropy](https://github.com/Albuso0/entropy) by [Yihong Wu](http://www.stat.yale.edu/~yw562/) and [Pengkun Yang](https://sites.google.com/site/pyangece/). 

Comprehensive script
---------
We provide ```main_entropy.py``` as an example script for our private estimator. In this script, we compare performance for these entropy estimators on different distributions including uniform, a distribution with two steps, Zipf(1/2), a distribution with Dirichlet-1 prior, and a distribution with Dirichlet-1/2 prior. We use RMSE (root-mean-square error) to indicate the performance of the estimator.

### Program arguments

* ```k int```: Set alphabet size. 
* ```eps float```: Set privacy parameter eps.
* ```l_degree int```: Set polynomial degree for private poly. Default *L=1.2 log k*.
* ```M_degree float```: Set the right endpoint of approximation interval for private poly. Default *M=2.0 log k*.
* ```N_degree int```: Set the threshold to apply polynomial estimator for private poly. Default *M=1.6 log k*.

For the parameter of poly estimator, we just use the default values in their code, which are *L=1.6 log k, M=3.5 log k, N=1.6 log k*. Please see [entropy](https://github.com/Albuso0/entropy) for more information.


Support coverage estimator
================
In this project, we implement our estimator for private estimation of support coverage.
This is a privatized version of the Smoothed Good-Toulmin (SGT) estimator of [Alon Orlitsky](http://alon.ucsd.edu/), [Ananda Theertha Suresh](http://theertha.info/), and [Yihong Wu](http://www.stat.yale.edu/~yw562/), from their paper [Optimal prediction of the number of unseen species](http://www.pnas.org/content/113/47/13283?sid=c704d36c-5237-4425-84e4-498dcd5151b1).
We compare the performance of the private and non-private statistics on both synthetic data and real-world data, including US Census name data and a text corpus from Shakespeare's Hamlet.

Some of our code is based off the SGT implementation of Orlitsky, Suresh, and Wu, graciously provided to us by Ananda Theertha Suresh. Specific files used are indicated in comments, and ```hamlet_total.csv``` is a reformatted version of a provided file. ```lastnames_total.csv``` is a subsampling of a file of [Frequently Occurring Surnames from the Census 2000](https://www.census.gov/topics/population/genealogy/data/2000_surnames.html).

Synthetic data 
---------
We provide ```main_synthetic.py``` as an example for our private estimator. In this script, we compare performance for the estimators on different distributions including uniform, a distribution with two steps, Zipf(1/2), a distribution with Dirichlet-1 prior, and a distribution with Dirichlet-1/2 prior. We use RMSE (root-mean-square error) to indicate the performance of the estimator.

### Program arguments
* ```k int```: Set alphabet size for the distribution.
* ```eps_index list```: Set privacy parameter for the private estimators.
* ```n int```: The number of the seen samples

Real data
---------
We provide ```main_real.py``` as an example for our private estimator. In this script, we compare performance for the estimators on real data. We use RMSE (root-mean-square error) to indicate the performance of the estimator.

### Program arguments
* ```file_name string```: The name of the histogram file. The histogram file must be in ```.csv```, which has only one column and each row is the number of samples for each species. We provide ```hamlet_total.csv``` and ```lastnames_total.csv``` as some examples.
* ```eps_index list```: Set privacy parameter for the private estimators.

Reference
================
For detailed explanations and analysis of our estimators, please refer to our paper [INSPECTRE: Privately Estimating the Unseen]().
