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

nn(nn==0) = [];

icts = unique(nn);
if(length(icts) == 1)
    mm = double(length(nn));
else
    if isinteger(nn)
	% MATLAB hist cannot convert int8, uint32, etc
	mm = zeros(length(icts), 1);
	for k = 1:length(icts)
	    mm(k) = sum(nn == icts(k));
	end
	% convert to double because gammaln/etc cannot take integers
	icts = double(icts);
    else
	mm = hist(nn,icts)';
    end
end

if size(icts,2)>1
    icts = icts';
end

mm = double(mm);
icts = double(icts);
