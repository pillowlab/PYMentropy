function [Hdir,Hvar] = computeHdir(mm,icts,K,alphas)
% [Hdir,Hvar] = computeHdir(mm,icts,K,alphas);
%
% Compute posterior mean of P(H|n,alpha), the expected entropy under a
% fixed Dirichlet prior with Dirichlet parameter alpha
%
% Accepts a vector "alphas" and returns posterior mean at each alpha
%
% Inputs:
%   mm = multiplicities (mm(j) is # bins with icts(j) elements)
%   icts = vector of unique counts
%   K = # total bins in distribution
%   alphas = scalar (or vector) of Dirichlet parameters 
%
% Ouptuts:
%   Hdir = mean entropy at each alpha
%   Hvar = variance of entropy at each alpha
% 
% Copyright Pillow Lab 2011. All rights reserved.

if isempty(mm) || isempty(icts)
    icts = 0;
    mm = K;
end


% modify mm and icts to account for #zeros (finite bin case)
nZ = K-sum(mm);
if( nZ > 0 && ~sum(icts==0))
   mm(length(mm) +1) = nZ;
   icts(length(icts)+1) = 0;
end

% Make alphas,icts,mm column vectors
if size(alphas,1)==1
    alphas = alphas';  % column vec
end
if size(icts,1)==1;
    icts = icts';  
end
if size(mm,1)==1;
    mm = mm'; 
end

N = icts'*mm; % number of samples
A = N+K*alphas; % number of effective samples (vector)
aa = bsxfun(@plus,alphas,icts'); % posterior Dirichlet priors

% Compute posterior mean over entropy
Hdir = digamma(A+1) - (1./A).*((aa.*digamma(aa+1))*mm);

%  Compute variance
if nargout > 1
    % NOTE: this is only prior variance!  (Need correct formula given data)
    %Hvar = - trigamma(K*alphas+1) + (alphas+1)./(K*alphas+1).*trigamma(alphas+1); 
    Hvar = reduced_varianceEntropy(mm,icts,alphas,K);
    Hvar = max(Hvar,0);
end
