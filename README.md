PYM entropy estimator MATLAB reference implementation
=====================================================

This code lives in github: https://github.com/pillowlab/PYMentropy

This is a reference implementation in MATLAB of the entropy estimator based on a Pitman-Yor mixture (PYM) prior. For the details behind how we derive the estimator, see the following papers:

- Evan Archer, Il Memming Park, Jonathan W. Pillow. Bayesian estimation of discrete entropy with mixtures of stick breaking priors. Neural Information Processing Systems [(NIPS) 2012](http://books.nips.cc/nips25.html)
- Evan Archer, Il Memming Park, Jonathan W. Pillow.  Bayesian Entropy Estimation for Countable Discrete Distributions. Journal of Machine Learning Research (JMLR), 15(81):2833−2868, 2014. http://arxiv.org/abs/1302.0328; http://jmlr.org/papers/v15/archer14a.html)

Quick example
=============
Let's estimate entropy from a sequence of natural numbers using the PYM estimator (in practice, the sequence would be samples drawn from some unknown distribution):

    >> [mm, icts] = multiplicitiesFromSamples([3 2 4 3 1 4 2 4 4]);
    >> [Hbls, Hvar] = computeH_PYM(mm, icts)

    Hbls =
	2.1476

    Hvar =
	0.2715

Here, `Hbls` is the Bayes least squares estimate of entropy, and `Hvar` is the posterior variance of the estimate. The units of this toolbox is *nats* (natural logarithm); to convert to *bits*, divide the result by `log(2) = 0.6931...`.

Requirements and Installation
=============================
To download go to github: https://github.com/pillowlab/PYMentropy/archive/master.zip

You must have Optimization toolbox (for `fmincon`).
To install, just add the package to your MATLAB path.
This package is developed under 7.13.0.564 (R2011b).

If using an older version of MATLAB, you may need [lightspeed](http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/) for fast digamma and polygamma implementations.

Run the unit test script `test_HPYM_randomized.m` to check if your copy is working fine.

License
=======
This package is distributed under the BSD license. See LICENSE.txt for details.

Converting data
===============
We represent raw data in three distinct ways, each at a different level of abstraction. The raw data itself are samples of discrete symbols from some unknown distribution. We do not operate directly upon the samples; rather, all entropy estimators provided take a succint representation called *multiplicities*. The multiplicity representation consists of two vectors of the same length: `mm` contains the number of symbols with the same number of occurrences, and `icts` contains the corresponding number of occurrences. In short, `mm(j)` is number of symbols with `icts(j)` occurrences. For example, the symbol sequence `[a b b c c d d d d]` has a multiplicity representation where `mm = [1 2 1]` and `icts = [1 2 4]`; this tells us that there is only one symbol with a single occurrence (`a`), two symbols with two occurrences (`b` and `c`), and one symbol with 4 occurrences (`d`). The ordering in `mm` and `icts` is not meaningful. We provide the following utility functions to convert data to the multiplicities. For large data, it is recommended that you write your own code that can exploit your data structure to convert it to the multiplicity representation in a memory efficient manner.

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
Discrete entropy does not care about the symbol identity. Since we assume independent and identically distributed samples, the ordering is also irrelvant. Hence, a discrete histogram also contains all information about the data. To convert a histogram to multiplicities, use `multiplicitiesFromCounts`.

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
If you need the histogram of your multiplicity representation, use `multiplicitiesToCounts`.

    >> hg = multiplicitiesToCounts(mm, icts)
    hg =
         4
         2
         2
         1

Spike train
-----------
For a specialized case when the data are contained in a cell array of spike timings from simultaneously recorded neurons, the multiplicities can be extracted using `multisptimes2words.m`. See also `fastWords2Counts`, `discreteTimeSeries2Words`, and `words2multiplicities.m`.
Alternatively, [Nemenman](http://nsb-entropy.sourceforge.net/) also has a fast C++ implementation for extracting multiplicities from time series.

Entropy estimation
==================
Simply call:

    [Hbls, Hvar] = computeH_PYM(mm, icts, prior)

where `mm` and `icts` form the multiplicity representation, and `prior` is an optional structure that rougly specifies the range of power-law tail behavior. If omitted, the default value is used.
Note that PYM estimator is not finite if there are less than 2 coincidences (1 coincidence = same sample observed twice).

The package also includes a few more entropy estimators. All functions of the form `computeH_*` take the multiplicity representation and return an estimate of `H` as well as a variance (if it is supported).

Controlling the Prior
---------------------
Theoretically, the PYM estimator has a degree freedom in specifying the prior in the gamma direction (details in the papers). Basically, any bounded, non-negative valued function `q(g = gamma)` defined on [0, Inf] can be used as a prior. The default prior is `q(g) = exp(-10/(1-g))`, which places more prior probability mass on distributions with light-tails, and limits the amount of prior mass placed on extremely heavy-tailed distributions.

Use `prior = pymPriorFactory(priorName, param)` to generate prespecified alternative priors, such as `q(g) = exp(-g)`. To use a custom prior, you must specify function handles for the prior function, `q(g)`, as well as its first and second derivatives. Once these function handles are placed in a structure (see pymPriorFactory.m for details), your custom prior may be used in a manner identical to that of the built-in priors.

**Enjoy!**
