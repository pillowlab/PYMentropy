function [mm icts wordvec] = multisptimes2words(sptimes, dt, T, Nbase)
% out = multisptimes2words(sptimes, dt, Nbase)
%
% Return unique integer symbols representing spike-time occurence 
% sequences in spike times from multiple cells.
% 
% INPUT: 
%   sptimes - [Ncells x1] cell array of spike times 
%        dt - [1x1] length of bins in times ( #bins will be T/dt)
%         T - [1x1] length of spike time recording 
%               * When T is empty (T = []), T is the maximum time
%                 occurring in sptimes. 
%     Nbase - [default: 2] Spike counts are considered
%             digits in integer system of 
%             base Nbase, with a binary system assumed by
%             default. Any spike count > (Nbase-1) is set to 
%             (Nbase - 1);
% OUTPUT: 
%       mm - multiplicities (mm(j) is number of bins with icts(j) samples)
%            * mm is sorted so that mm(j) >= mm(j+1)
%     icts - unique sample counts
%  wordvec - [T/dt x 1] vector of integer words (ie, the word wordvec(i) 
%            occurrs at time i).
%
%

    if(nargin < 4)
        Nbase = 2;
    end

    Ncell = length(sptimes);
    assert(T>dt, 'Bin size dt must be less than total time T.')
    % Generate appropriately-sized bins
    xx = 0:dt:T;

    Spvec = nan(length(xx),Ncell);

    % now bin all cells
    for cdx = 1:Ncell
        thecell = sptimes{cdx}; % this is a vector of spike times, in ms
        Spvec(:,cdx) = histc(thecell, xx);
    end

    % If we have more spikes in a bin that can be represented by our chosen
    % base, set to the largest digit representable by our base
    Spvec( Spvec > Nbase-1 ) = Nbase - 1;

    % Create vector describing our base.
    basevec = Nbase.^(0:Ncell-1)';

    % Finally, make a vector of words
    wordvec = Spvec*basevec;
    
    nn = hist(wordvec, unique(wordvec));
    [mm icts] = multiplicitesFromCounts(nn);
    
    [~, vv] = sort(mm, 'descend');
    mm = mm(vv); icts = icts(vv);
    
end