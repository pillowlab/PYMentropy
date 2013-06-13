function nn = multiplicitiesToCounts(mm,icts)
% nn = multiplicitiesToCounts(mm,icts)
%
% Compute the counts (i.e., histogram bins) from the "multiplicities"
% representation of a dataset
%
% INPUT: 
%    mm    - multiplicities (mm(j) is number of bins with icts(j) samples)
%    icts  - unique sample counts
%
% OUTPUT:
%    nn - nn(j) is the number of samples in the j'th (sorted) bin
%

nbins = sum(mm);
nn = zeros(nbins,1);

ibin = 0; 
for j = 1:length(icts)
    nn(ibin+1:ibin+mm(j))=icts(j);
    ibin = ibin+mm(j);
end

% flip so that largest counts come first
nn = flipud(nn);
