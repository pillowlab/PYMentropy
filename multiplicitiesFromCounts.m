function [mm,icts] = multiplicitiesFromCounts(nn)
% [mm,icts] = multiplicitiesFromCounts(nn);
%
% Compute the "multiplicities" from a histogram of counts nn.
%
% Creates a histogram of the unique (integer count) values in
% nn, and return just the non-zero entries.
%
% INPUT: 
%    nn - vector of counts; nn(j) is the number of samples in the j'th bin
%
% OUTPUT:
%    mm    - multiplicities (mm(j) is number of bins with icts(j) samples)
%    icts  - unique sample counts
%
% $Id: multiplicitiesFromCounts.m 1201 2012-04-17 07:54:53Z evan $

nn(nn==0) = [];

icts = unique(nn);
if(length(icts) == 1)
    mm = length(nn);
else
    mm = hist(nn,icts)';
end

if size(icts,2)>1
    icts = icts';
end
