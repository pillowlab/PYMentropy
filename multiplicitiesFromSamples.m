function [mm,icts] = multiplicitiesFromSamples(x)
% [mm,icts] = multiplicitiesFromSamples(x);
%
% Compute the "multiplicities" in a sample x.
%
% Creates a histogram of the unique (integer count) values in
% the histogram of x, and returns just the non-zero entries.
%
% INPUT: 
%    x - samples from some distribution over the integers
%
% OUTPUT:
%    mm    - multiplicities (mm(j) is number of bins with icts(j) samples)
%    icts  - unique sample counts
%
% $Id: multiplicitiesFromSamples.m 2864 2013-02-19 22:39:16Z memming $ 

ux = unique(x);

if length(ux) == 1
    icts = length(x);
    mm = 1;
else
    if isinteger(nn)
    else
	nn = hist(x, ux);
	[mm, icts] = histc(nn, unique(nn));
    end
end

if size(icts,2) > 1
    icts = icts';
end

if size(mm,2) > 1
    mm = mm';
end
