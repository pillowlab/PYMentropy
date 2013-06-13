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

ux = unique(x);

if length(ux) == 1
    icts = length(x);
    mm = 1;
else
    nn = histc(x, ux);
    [mm, icts] = hist(nn, unique(nn));
end

if size(icts,2) > 1
    icts = icts';
end

if size(mm,2) > 1
    mm = mm';
end
