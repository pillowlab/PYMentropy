function [Hpy,Hvar] = computeHpyPrior(alphas,ds)
% Compute prior mean of P(H|alpha,d) under a Pitman-Yor process prior
% with parameter alpha and d.
%
% Accepts same size "alphas" and "ds" and returns posterior mean
% at each alpha and d with the same size.
%
% Ouptuts:
%   Hpy = mean entropy at each alpha
%   Hvar = variance of entropy at each alpha
% 
% Copyright 2012 Pillow Lab. All rights reserved.

Hpy = digamma(1+alphas) - digamma(1-ds);
Hvar = (alphas+ds) ./ ((alphas+1).^2 .* (1-ds)) + (1-ds)./(alphas+1) .* trigamma(2 - ds) - trigamma(2 + alphas);
