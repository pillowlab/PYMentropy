function [Hbls, Hvar, SummaryStr] = computeH_PYM_v4(mm, icts, param, verbose)
% [Hbls, Hvar] = computeH_PYM_v4(mm, icts, param, verbose)
% Compute PYM estimate using either adaptive quadrature or an adaptive grid.
% The area of integration is chosen by finding the MAP (peak of the posterior)
% and expanding on the directions of the eigenvector of the estimated Hessian.
%
% INPUT:
%       mm - multiplicites
%     icts - counts for each multiplicty
%    param - (optional) parameter structure
%	  param.pymPrior: (optional) This is a structure obtained from pymPriorFactory.
%         param.returnEB: (optional) If true, return the empirical bayes
%                       estimate of entropy. That is, return the posterior
%                       mean of entropy given the parameters (alpha,d) at the MAP.
%	  param.verbose: (optional) if true, outputs debugging info
%	  param.nagrid: (optional) Integer > 0. Specify number of gridpoints
%	                to use in integral over alpha. 
%	  param.ndgrid: (optional) Integer > 0. Specify number of gridpoints
%	                to use in integral over d. 
%	  param.method: (optional) Specify integration method.
%                       = 1:  Use Chebyshev spectral integration
%                       = 2:  Use Gauss-Legendre quadrature
%  verbose -  (optional) if true, outputs debugging info
%
% OUTPUT:
%  Hbls      - posterior mean estimate of entropy
%	     (Inf or NaN is returned in the pathological case that the data
%	     have no coincidences, or when only a single symbol is observed.)
%  Hvar      - posterior variance
% SummaryStr - a structure containing variables relevant to computation
%
% Outline of algorithm:
%   - Find the MAP estimate of (d, alpha) given data and prior
%   - Estimate the Hessian of the posterior over (d, alpha)
%   - IF the hessian indicates the posterior is concentrated, 
%         * Estimate the region of integration to be a rectangular region around the
%	    MAP estimate and extends proportional to the Hessian.
%         * Use Gauss-Legendre quadrature to integrate.
%         * If param.nagrid and param.ndgrid are not specified, by default:
%                       nagrid = ndgrid = 25.
%   - ELSE the posterior is not well-concentrated
%         * Use Chebyshev spectral integration method (Boyd 1986) for
%           integration
%         * If param.nagrid and param.ndgrid are not specified, by default:
%                       nagrid = 250
%                       ndgrid = 25
%   - Use selected integration method to compute the following:
%	- Z: normalizer
%	- Hbls: expected entropy given (d, alpha)
%	- Hvar: variance of entropy given (d, alpha)
%
% See also: pymPriorFactory
% Requires: optimization toolbox
%
% Copyright 2012 Pillow lab. All rights reserved.

warning('PYM:deprecated', 'Deprecated. Use computeH_PYM.m instead');

[Hbls, Hvar, SummaryStr] = computeH_PYM(mm, icts, param, verbose);
