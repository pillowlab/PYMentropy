function [Hpy, Hvar, dbgstr] = computeHpy(mm,icts,alphas,ds)
% [Hpy,Hvar] = computeHpy(mm,icts,alphas,ds)
%
% Compute posterior mean of P(H|n,alpha,d) under a Pitman-Yor process prior
% with parameter alpha and d.
%
% Accepts same size "alphas" and "ds" and returns posterior mean
% at each alpha and d with the same size.
%
% Posterior process is Dir(n_1-d, ..., n_K-d, alpha+d*K)
% So the first size biased sample is distributed as a mixture of Beta's:
%
% sum_{k=1}^K (n_k-d)/(n+alpha) Beta(1+n_k-d, n+alpha-n_k+d*(K-1))
% + (alpha+d*K)/(n+alpha) Beta(1-d, alpha+n+d*K) % 2->1 +n
%
% Hence the entropy is,
% sum_{k=1}^K (n_k-d)/(n+alpha) 
% (digamma(1+n+alpha+d*(K-2)) - digamma(1+n_k-d))
% + (alpha+d*K)/(n+alpha) 
% (digamma(1+alpha+n+d*(K-1)) - digamma(1-d))
%
% Inputs:
%   mm = multiplicities (mm(j) is # bins with icts(j) elements)
%   icts = vector of unique counts
%   K = # total bins in distribution
%   alphas = scalar (or vector) of concentration parameter
%   ds = scalar (or vector) of discount parameter
%
% Ouptuts:
%   Hpy = mean entropy at each alpha
%   Hvar = variance of entropy at each alpha
% 
% Copyright Pillow Lab 2011-2012. All rights reserved.

if isempty(mm) || isempty(icts)
    icts = 0;
    mm = 0;
end

if nargin < 4
    % We haven't passed any d's - we want Dirichlet Process posterior
    ds = zeros(size(alphas));
end
    
% make alphas,icts,mm column vectors
originalSize = size(alphas);
alphas = alphas(:); ds = ds(:); icts = icts(:); mm = mm(:);

N = icts'*mm; % number of samples
K = sum(mm(icts>0)); % number of tables
% assert(K > 0);

Hp = zeros(size(ds)); 
% Hp := E[h(p)]  := A (in notes)
% Hpi := E[h(pi)] := B
% Hpstar := E[h(p_*)]
[Hpi,HvarPi] = computeHpyPrior(alphas+K*ds,ds);  

if N == 0
    % if we have no data, return prior mean & variance
    Hpy = Hpi;
    Hvar = HvarPi;
    return
end

% compute E[p_*] and E[(1-p_*)]
oneminuspstarmean= (N-K*ds)./(alphas+N);
pstarmean  = (alphas+K*ds)./(alphas+N);
% compute E[(p_*).^2] and E[(1-p_*)^2]
pstarsq         = (alphas+K*ds).*(alphas+K*ds+1)./( (alphas+N).*(alphas+N+1) );
oneminuspstarsq = (N-K*ds).*(N-K*ds+1)./( (alphas+N).*(alphas+N+1) );
% compute E[h(p_*)]
Hpstar = digamma(alphas+N+1) - ...
	(alphas+K*ds)./(alphas+N).*digamma(alphas+K*ds+1) - ...
	(N-K*ds)./(alphas+N).*digamma(N-K*ds+1);

for k = 1:numel(ds)    
    Hp(k) = digamma(N-K*ds(k)+1) - ...
	    mm'*((icts-ds(k)) ./ (N-K*ds(k)) .* digamma(1+icts-ds(k)));
end

Hpy = oneminuspstarmean .* Hp + pstarmean .* Hpi + Hpstar;
Hpy = reshape(Hpy, originalSize);

dbgstr.Hp = Hp;
dbgstr.Hpi = Hpi;
dbgstr.Hpstar = Hpstar;
dbgstr.Hpy = Hpy;
dbgstr.pstarmean = pstarmean;
dbgstr.oneminuspstarmean = oneminuspstarmean;

%  Compute variance
if nargout > 1
    %% First term term of variance decomposition
    % First, compute pstarHpstar
    % compute E[p_*h(p_*)]
    pstarHpstar = (alphas + K*ds) .* (alphas + K*ds + 1) ./ ...
	( (alphas+N) .* (alphas+N+1) ) .* (digamma(alphas+N+2) - ...
	digamma(alphas+K*ds+2)) + ...
        (alphas+K*ds).*(N-K*ds) ./ ( (alphas + N) .* (alphas+N+1) ) .* ...
	(digamma(alphas+N+2) - digamma(N-K*ds+1));

    % E[p_*^2]
    Hmom2 = zeros(size(ds));
    % When alpha is huge , the second moment of the entropy should just be
    % zero. This is producing an error as a result of numerical troubles.
    for k = 1:numel(ds)
	[~, Hmom2(k)] = varianceEntropy(0, [alphas(k)+K*ds(k); N-K*ds(k)]);
    end

    % E[K]^2 (NOTE: Bad notation! This is not #bins, this is K(p_*) used in the notes!)
    Kmeansq = (oneminuspstarmean.*Hp + pstarmean.*Hpi + Hpstar).^2;
    % E[K^2]
    Ksqr = 2*pstarHpstar .* (Hpi - Hp) + 2*Hp.*Hpstar + Hmom2 + ...
        pstarsq.*(Hpi.^2 - 2*Hpi.*Hp) + ...
        2*(alphas+K*ds)./(alphas+N).*Hp.*Hpi + oneminuspstarsq.*Hp.^2;
    
    Trm1 = Ksqr - Kmeansq;

    %% Second term of variance decomposition
    % HvarPi is computed above. 
    HvarP = reduced_varianceEntropy(mm, icts, -ds, K);

    Trm2 = (N-K*ds).*(N-K*ds+1)./((alphas+N).*(alphas+N+1)).*HvarP + ...
	    (alphas+K*ds).*(alphas+K*ds+1)./((alphas+N).*(alphas+N+1)).*HvarPi;

    %% Add two variance terms
    Hvar = Trm1 + Trm2;

    dbgstr.pstarHpstar = pstarHpstar;
    dbgstr.pstarsq = pstarsq;
    dbgstr.oneminuspstarsq = oneminuspstarsq;
    dbgstr.Hmom2 = Hmom2;
    dbgstr.HvarP = HvarP;
end
