# INSPECTRE

We provide implementation of our private entropy estimators and private support coverage estimator in Python. This is a tutorial for the people who want to use or develop upon these estimators.

Table of contents
=================
* [Table of contents](#table-of-contents)
* [Prerequisites](#prerequisites)
* [Entropy estimators](#entropy-estimators)
* [Support coverage estimator](#support-coverage-estimator)
* [References](#reference)

Prerequisites
=====
This project works on Python3, and we use the package of [Numpy](http://www.numpy.org), the fundamental package for scientific computing with Python. We also use [matplotlib](https://matplotlib.org/index.html) for the plot. Before running the code, make sure Python3, Numpy and Matplotlib are properly installed. 

* [Instruction for installing Python3](https://docs.python.org/3/using/index.html)
* [Instruction for installing Numpy](https://www.scipy.org/install.html)
* [Instruction for installing Matplotlib](https://matplotlib.org/users/installing.html) 


Entropy estimators
================
In this project, we compare the performance of our private entropy estimator with a number of non-private estimators, which include the plug-in estimator, the Miller-Madow Estimator, and the sample optimal polynomial approximation estimator from the paper [Minimax Rates of Entropy Estimation on Large Alphabets via Best Polynomial Approximation](http://ieeexplore.ieee.org/abstract/document/7444171/). 

We provide implementation of our private entropy estimators, including the private plug-in estimator and the private poly estimator. For the other non-private estimators, we simply borrow the code from the project [entropy](https://github.com/Albuso0/entropy) from github, Please see their project website for more information.

Comprehensive script
---------
We provide ```main_entropy.py``` as an example script for our private estimator. In this script, we compare performance for these entropy estimators on different distributions including uniform, a distribution with two steps, Zipf(1/2), a distribution with Dirichlet-1 prior, and a distribution with Dirichlet-1/2 prior. We use RMSE (root-mean-square error) to indicate the performance of the estimator.

### Program arguments

* ```k int```: Set alphabet size. 
* ```eps float```: Set privacy parameter eps.
* ```l_degree int```: Set polynomial degree for private poly. Default *L=1.2 log k*.
* ```M_degree float```: set the right endpoint of approximation interval for private poly. Default *M=2.0 log k*.


Support coverage estimator
================
In this project, we compare the performance of our private estimator with [Smoothed Good-Toulmin estimator (SGT)]() on both synthetic data and real data. We provide implementation of our support coverage estimator which is a privatized version of their statistic. For the SGT, we borrow the code from which [SGT]() used for the experiment part.

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
* ```file_name string```: The name of the histogram file. The histogram file must be in ```.crv```, which has only one column and each row is the number of samples for each species.
* ```eps_index list```: Set privacy parameter for the private estimators.

Reference
================
For detailed explanation of the estimator, please refer to our paper [INSPECTRE: Privately Estimating the Unseen]().
