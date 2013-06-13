function [eta,gam] = convertfromPYparams(alpha,d)
% [eta,gam] = convertfromPYparams(alpha,d)
%
% Convert natural Pitman-Yor process parameters to entropy-separable
% parameters
%
% % INPUT:
%   alpha (-  [0, \infty]  (Dirichlet parameter; total prior measure).
%   d     (-  [0, 1]       (discount parameter; controls tails).
% 
% OUTPUT
%   eta (-  [0, \infty]  (controls expected entropy)
%   gam (-  [0, 1]       (controls shape of tails)
%
% See also: pymPriorFactory, computeH_PYM_quad
%
%
% Dependencies: digamma.m

assert(all(alpha(:) > 0), 'Concentration param should be non-negative');
assert(all(d(:) >= 0 & d(:) < 1), 'Discount must be in [0,1)');

eta = digamma(alpha+1)-digamma(1-d);
gam = (digamma(1)-digamma(1-d))./eta;
