# iSTAC
information-theoretic Spike-Triggered Average and Covariance (iSTAC)
estimator for neural receptive fields

**Description:** Estimates a set of linear filters that best capture a
 neuron's input-output properties, using an information-theoretic
 objective that optimally combines information from the
 spike-triggered average and spike-triggered covariance. The filters
 can be considered as the first stage in a linear-nonlinear-Poisson
 (LNP) model of the neuron's response. They are sorted by
 informativeness, providing an estimate of the mutual information
 gained by the inclusion of each filter.

**Relevant publication:**
Pillow & Simoncelli, *Journal of Vision*
2006. [[abs](http://pillowlab.princeton.edu/pubs/abs_iSTAC_JOV06.html),
      [pdf](http://pillowlab.princeton.edu/pubs/pillow_simoncelli_iSTAC_JOV06.pdf)]

Download
==========

* **Command line**: clone the repository from github (e.g., ```git
  clone git@github.com:pillowlab/iSTAC.git``` )
* **Browser**:  download zipped archive:  [iSTAC-master.zip](https://github.com/pillowlab/iSTAC/archive/master.zip)


Usage
=====

* Examine the script `test_iSTAC_script.m` for a line-by-line tutorial
on how to use the code contained in this package, which goes through
several simulated examples.

* The primary function used for estimating the filters is `compiSTAC.m`
