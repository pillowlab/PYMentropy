function [v moment2 m] = varianceEntropy(alpha, n)
% Analytical computation of variance of entropy under Dirichlet prior
% Input
%   alpha: (1xM) parameter of Dirichlet distribution
%   n: (Kx1) number of observation per bin
%
% Output
%   v: (1xM) variance of entropy
%   moment2: (1xM) second moment of entropy
%   m: (1xM) mean of entropy
%
% Note that this function doesn't call digamma/trigamma of lightspeed but uses psi of MATLAB
% For computing part of PY posterior, alpha is -d.
%
% See Also: reduced_varianceEntropy
%
% Copyright Pillow Lab 2010. All rights reserved.

n = n(:);
alpha = alpha(:)';
K = length(n);
N = sum(n);
aK = alpha * K;

% Mean entropy squared
nAlpha = n * ones(1, length(alpha)) + ones(K, 1) * alpha;
m = psi(0, N + aK + 1) - sum((nAlpha) .* psi(nAlpha + 1)) ./ (N + aK);
m2 = m.^2;

dgNK = psi(0, N + aK + 2);
tgNK = psi(1, N + aK + 2);
dgAlpha = psi(0, nAlpha + 1);
tgAlpha = psi(1, nAlpha + 2);

c = (N + aK + 1) .* (N + aK);
%crossTerms = sum((nAlpha).^2) .* ...
    %sum((dgAlpha - dgNK).^2) - tgNK;
for ka = 1:length(alpha)
    % CAUTION! below statement produces outer product of a potentially high
    % dimensional vector n.
    tempCrossTerms = nAlpha(:, ka) * nAlpha(:, ka)' .* ...
	((dgAlpha(:, ka) - dgNK(:, ka)) * (dgAlpha(:, ka) - dgNK(:, ka))' - tgNK(:, ka));
    tempCrossTerms = tempCrossTerms - diag(diag(tempCrossTerms));
    crossTerms(ka) = sum(tempCrossTerms(:));
end

diagTerms = (nAlpha + 1) .* (nAlpha) .* ...
    ((dgAlpha + 1./(nAlpha + 1) - ones(K,1) * dgNK).^2 ...
    + tgAlpha - ones(K,1) * tgNK);
diagTerms = sum(diagTerms);

moment2 = (crossTerms + diagTerms) ./ c;
v = moment2 - m2;

if( v < -2*eps )
    warning('\nvarianceEntropy.m: Produced negative variance values.\n')
end

v = max(0,v); moment2 = max(moment2, 0);
    
%assert(all(v >= 0))

% n = 0 case
% nv0 = (K - 1) * alpha / (alpha * K + 1) * ((psi(0, alpha + 1) - psi(0, aK + 2))^2 - psi(1, aK + 2)) + (alpha + 1) / (aK + 1) * ((psi(0,alpha + 2) - psi(0, aK + 2))^2 + psi(1, alpha + 2) - psi(1, aK + 2)) - m2

% nv0true = (alpha + 1)/(aK + 1) * psi(1, alpha + 1) - psi(1, aK + 1)
