function [alpha,d] = converttoPYparams(eta,gam)
% [alpha,d] = converttoPYparams(eta,gam);
%
% Convert from entropy-separable to "natural" parameters for the Pitman
% process.
%
% INPUT:
%   eta (-  [0, \infty]  (controls expected entropy)
%   gam (-  [0, 1]       (controls shape of tails)
% 
% OUTPUT
%   alpha (-  [0, \infty]  (Dirichlet parameter; total prior measure).
%   d     (-  [0, 1]       (discount parameter; controls tails).
%
% Dependencies: digamma.m, invdigamma.m
%

alpha = invdigamma(eta.*(1-gam)+digamma(1))-1;
d = 1-invdigamma(digamma(1)-eta.*gam);

if any(isinf(alpha(:))) || any(isnan(alpha(:)))
    warning('alpha too large: returning inf / nan');
end
