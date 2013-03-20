function [mm icts nn bin] = words2multiplicities(words,Nbase)
% [mm icts nn wordvec] = words2multiplicities(words,Nbase)
%
% Compute multiplicities (# of bins with each unique number of samples)
% representation from spike words
% 
% INPUT: 
%   words - [Nwords x Ncells] matrix (each row is a word)
%   Nbase - [default: 2] Spike counts are considered digits in integer
%   system of base Nbase, with a binary system assumed by default. Any
%   spike count > (Nbase-1) is set to (Nbase - 1);
%
% OUTPUT: 
%       mm - multiplicities (mm(j) is number of bins with icts(j) samples)
%            note: mm is sorted so that mm(j) >= mm(j+1)
%     icts - unique sample counts
%       nn - sorted count of word frequencies
%      bin - [Ncellsx1] vector mapping word to its index in nn
%
% CAUTION: double can take up to 51 bits, but will break down after that!

if nargin < 2; Nbase = 2; end

wordlen = size(words, 2);

if wordlen * log2(Nbase) > 53
    error('Word too long. Try using fastWords2Counts.');
end

% If we have more spikes in a bin that can be represented by our chosen
% base, set to the largest digit representable by our base
words( words > Nbase-1 ) = Nbase - 1;

% Create vector describing our base.
basevec = Nbase.^(0:wordlen-1)';

% Make a vector of words
wordvec = words*basevec;

% histogram of data
[nn bin] = histc(wordvec, unique(wordvec));

% Compute multiplicities from histogram
[mm icts] = multiplicitiesFromCounts(nn);

% Sort the resulting multiplicities
[~, vv] = sort(mm, 'descend');
mm = mm(vv); 
icts = icts(vv);
   
obin = bin;
if nargout > 2
    [nn I] = sort(nn,'descend');
    for idx = 1:length(I)
        bin(obin == I(idx)) = idx;
    end
end
