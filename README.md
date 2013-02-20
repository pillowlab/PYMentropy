PYM entropy estimator MATLAB reference implementation
=====================================================

version: $Id: README_PYM.md.txt 2719 2013-02-01 00:38:29Z memming $

This is a reference implementation of the entropy estimator based on Pitman-Yor mixture (PYM) prior. For the details of how we derive the estimator see the following papers:

- Evan Archer, Il Memming Park, Jonathan W. Pillow. Bayesian estimation of discrete entropy with mixtures of stick breaking priors. Neural Information Processing Systems [(NIPS) 2012](http://books.nips.cc/nips25.html)
- Evan Archer, Il Memming Park, Jonathan W. Pillow.  Bayesian Entropy Estimation for Countable Discrete Distributions. (submitted, soon to be in arXiv)

Quick example
=============
Let's estimate entropy from some sequence of natural numbers using the PYM estimator.

    >> [mm, icts] = multiplicitiesFromSamples([3 2 4 3 1 4 2 4 4]);
    >> [Hbls, Hvar] = computeH_PYM_v4(mm, icts)

    Hbls =
	2.1476

    Hvar =
	0.2715

where `Hbls` is the Bayes least squares estimate of entropy, and `Hvar` is the posterior variance of the estimate. The units of this toolbox is *nats* (natural logarithm); to convert to *bits*, divide the result by `log(2) = 0.6931...`.

Requirements and Installation
=============================
You must have Optimization toolbox (for `fmincon`).
To install, just add the package to your MATLAB path.
This package is developed under 7.13.0.564 (R2011b).

If using an older version of MATLAB, you may need [lightspeed](http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/) for fast digamma and polygamma implementations.

Run the **test script (TBD)** to check if your copy is working fine.

License
=======
This package is distributed under the BSD license. See LICENSE.txt for details.

Converting data
===============
There are three levels of representation of raw data we consider. All entropy estimators provided take a succint representation called *multiplicities*. The multiplicity representation consists of two vectors of same length: `mm` contains the number of symbols with the same number of occurrences, and `icts` contains the corresponding number of occurrences. In short, `mm(j)` is number of symbols with `icts(j)` occurrences. For example, for the symbol sequence `[a b b c c d d d d]`, the multiplicity representation is `mm = [1 2 1]` and `icts = [1 2 4]`, because there is only one symbol with a single occurrence (`a`), two symbols with two occurrences (`b` and `c`), and one symbol with 4 occurrences (`d`). The ordering in `mm` and `icts` is not meaningful. We provide the following utility functions to convert data to the multiplicities. However, for large data, it is recommended that you write your own code that can exploit your data structure to convert it to the multiplicity representation in a memory efficient manner.

Raw samples
-----------
Given a vector of symbols with unique numerical representation, use `multiplicitiesFromSamples`.

    >> [mm, icts] = multiplicitiesFromSamples([1 2 2 3.5 3.5 4 4 4 4])
    mm =
         1
         2
         1
    
    icts =
         1
         2
         4

Histogram representation
------------------------
Discrete entropy does not care about the symbol identity. Since we assume independent and identically distributed sample, the ordering is also irrelvant. Hence, a discrete histogram also contains all information about the data. To convert a histogram to multiplicities, use `multiplicitiesFromCounts`.

    >> [mm, icts] = multiplicitiesFromCounts([1 2 2 4])
    mm =
         1
         2
         1
    
    icts =
         1
         2
         4

Converting back to histogram
----------------------------
For some reason, if you need to see the histogram of your multiplicity representation, use `multiplicitiesToCounts`.

    >> hg = multiplicitiesToCounts(mm, icts)
    hg =
         4
         2
         2
         1

Spike train
-----------
For a specialized case when a cell array of spike timings from a simultaneously recorded neurons is given, the multiplicities can be extracted using `multisptimes2words.m`. [Nemenman](http://nsb-entropy.sourceforge.net/) also has a fast C++ implementation for extracting multiplicities from time series.

Entropy estimation
==================
Simply call:

    [Hbls, Hvar] = computeH_PYM_v4(mm, icts, prior)

where `mm` and `icts` are the multiplicity representation, and `prior` is an optional structure that rougly specifies the range of power-law tail behavior. If omitted, the default value is used.

The package also includes a few more entropy estimators. The functions in the form `computeH_*` all takes the multiplicity representation and returns an estimate of `H` and also a variance if it is supported.

Controlling the Prior
---------------------
Theoretically the PYM estimator has a degree freedom in specifying the prior in the gamma direction (details in the papers). Basically, any non-diverging, non-negative valued function `q(g = gamma)` defined on [0, Inf] can be used as a prior. The default is `q(g) = exp(-10/(1-g))` which gives more emphasis on light tailed or power-law tails with small-exponent distributions, and limit extremely heavy tailed power-law tails from the prior.

Use `prior = pymPriorFactory(priorName, param)` to generate prespecified alternative priors such as `q(g) = exp(-g)`. To use a custom prior, you need to specify the function handles for prior function `q(g)`, its first, and second derivatives in a structure (see pymPriorFactory.m for details).

**Enjoy!**
