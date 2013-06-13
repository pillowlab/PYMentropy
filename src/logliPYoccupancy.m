function [logp dlogp ddlogp] = logliPYoccupancy(alphas,ds,mm,icts)
% logp \propto logliPYoccupancy(alpha,d,mm,icts)
% OR
% logp \propto logliPYoccupancy(alpha,d,x)
%
% Compute log probability of multiplicities mm
% on tables with counts icts after after receiving dot(mm,icts)
% customers, under Pitman-Yor Process with parameters alpha, d.
%
% Given by Antoniak 1974 (proposition 3) for DP
%
% Inputs:
%  alpha = concentration parameter
%     d  = discount paramter
%    mm  = multiplicities
%  icts  = ranks
% OR
%  alpha = concentration parameter
%      d = discount paramter
%      x = table of counts
%
% Output:
%  pk = unnormalized conditional probabiliy P(alpha, d | data)
%

if(nargin == 3)
    [mm icts] = multiplicitesFromCounts(mm);
end
if( isvector(alphas) && isscalar(ds) )
    ds = ds*ones(size(alphas));
elseif( isvector(ds) && isscalar(alphas) )
    alphas = alphas*ones(size(ds));
else
    assert(length(alphas) == length(ds))
end

N = dot(mm,icts);

mm = mm(icts>0);
icts = icts(icts>0);

K = sum(mm);

logp = zeros(size(alphas));
for adx = 1:length(alphas)
    alpha = alphas(adx); d = ds(adx);
    logp(adx) = -gammalndiff(alpha+1, N-1) ... % gamma(1+alpha)/gamma(n+alpha)
        + sum(log(alpha+(1:K-1)*d )) + ...
        dot(mm, gammalndiff(1-d, icts-1));
end

if(nargout>1)
    KK=1:(K-1);
    % we should do this in terms of digammas, but oh well for now.
    dlogp = zeros(length(alphas),2);
    for adx = 1:length(alphas)
        alpha = alphas(adx); d = ds(adx);
        Z = 1./(alpha+(KK)*d);
        dlogp(adx,1) = sum(Z) - digamma(alpha+N) + digamma(alpha+1);
        dlogp(adx,2) = dot(KK,Z) - dot(mm, digamma(icts-d) - digamma(1-d));
    end
%     KK=1:(K-1);
%     % we should do this in terms of digammas, but oh well for now.
%     dlogp = zeros(length(alphas),2);
%     for adx = 1:length(alphas)
%         alpha = alphas(adx); d = ds(adx);
%         dlogp(adx,1) = sum(1./(alpha+(KK)*d)) - sum(1./(alpha+(1:(N-1))));
%         for idx = 1:length(mm)
%             dlogp(adx,2) = dlogp(adx,2) ...
% 			    - mm(idx)*sum(1./((1:icts(idx)-1)-d)); %d
%         end
%         dlogp(adx,2) = dlogp(adx,2) + sum(KK./(KK*d+alpha)); 
%     end
end

if(nargout>2)
   ddlogp = zeros(length(alphas),3);
   for adx = 1:length(alphas)
       alpha = alphas(adx); d = ds(adx);
       Z = 1./(alpha+KK*d);
       Z2 = Z.^2; 
       % aa
       ddlogp(adx,1) = trigamma(1+alpha) - trigamma(alpha+N) - sum(Z2);
       % ad
       ddlogp(adx,2) = -sum(KK.*Z2);
       % dd
       ddlogp(adx,3) = -sum((KK.^2).*Z2) + dot(mm, trigamma(icts-d) - trigamma(1-d));
       if(any(isinf(ddlogp)))
           keyboard
       end
   end
   if(length(alphas) == 1)
       q = ddlogp;
       ddlogp = zeros(2);
       ddlogp(1,1) = q(1);
       ddlogp(1,2) = q(2); ddlogp(2,1) = q(2);
       ddlogp(2,2) = q(3);
   end
end
