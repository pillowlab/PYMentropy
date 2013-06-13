function v = reduced_varianceEntropy(mm, icts, alphas, K, flag)
% v = reduced_varianceEntropy(mm, icts, alphas, K, flag)
% Analytical computation of variance of entropy under Dirichlet prior
% Input
%   mm  - 
%   icts  - 
%   alpha - parameter of Dirichlet distribution (scalar or vector)
%   K     - # bins in distribution
%   flag  - 
%           0: return variance (default)
%           1: return 2nd moment
% 
% Output
%   v     - variance (or second moment) of entropy
%
% TBD: Vectorize multiple alpha (vector alpha). Currenly uses for loop.
%
% Copyright Pillow Lab 2011. All rights reserved.

v = zeros(numel(alphas), 1);

if(~exist('flag','var'))
    flag = 0;
end

N = dot(mm,icts);

% modify mm and icts to account for #zeros (finite bin case)
nZ = K-sum(mm);
if( nZ > 0 && ~sum(icts==0))
   mm(length(mm) +1) = nZ;
   icts(length(icts)+1) = 0;
end

% make sure everyone's a column vector
mm = mm(:);
alphas = alphas(:);
icts = icts(:);

dgNKs = psi(0, N + alphas * K + 2);
tgNKs = psi(1, N + alphas * K + 2);
dgAlphas = psi(0, icts * ones(1, numel(alphas)) + ones(numel(icts), 1) * alphas' + 1);
tgAlphas = psi(1, icts * ones(1, numel(alphas)) + ones(numel(icts), 1) * alphas' + 2);

for idx = 1:numel(alphas)
    alpha = alphas(idx);
    aK = alpha * K;
    dgNK = dgNKs(idx); dgAlpha = dgAlphas(:, idx);
    tgNK = tgNKs(idx); tgAlpha = tgAlphas(:, idx);

    zcross = mm;
    z1 = mm == 1;
    
    zcross(~z1) = mm(~z1).*(mm(~z1)-1)/2;
    zcross(z1) = 0;

    c = (N + aK + 1) .* (N + aK);
    crossTerms = (icts + alpha) * (icts + alpha)' .* ...
	((dgAlpha - dgNK) * (dgAlpha - dgNK)' - tgNK);
    q = crossTerms;
    crossTerms = (crossTerms - diag(diag(crossTerms))).*(mm*mm');
    crossTerms = sum(crossTerms(:)) + sum(sum(q.*diag(zcross)))*2;

    diagTerms = (icts + alpha + 1) .* (icts + alpha) .* ...
	((dgAlpha + 1./(icts+alpha+1) - dgNK).^2 + tgAlpha - tgNK);
    diagTerms = sum(diagTerms.*mm);

    v(idx) = (crossTerms + diagTerms) / c;
end

if(~flag)
    Sdir2 = computeHdir(mm,icts,K,alphas).^2;
    v = v - Sdir2;
end
