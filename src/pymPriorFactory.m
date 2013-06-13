function pymPrior = pymPriorFactory(priorName, param)
% pymPrior = pymPriorFactory(priorName, param)
%
% PYM has a free parameter which is a distribution over gamma.
% gamma controls the amount of probability mass on d vs a for the same
% mean expected entropy contour. This function is sometimes denoted as
% $q(\gamma)$ in our manuscripts.
%
% Input:
%   priorName: string indicating the prior
%	'triangle', 'exponential', 'evan', 'dp', 'identity'
%   param: additional parameters specific for each prior
%	'exponential' and 'evan' requires a time constant
% Output:
%   pymPrior.desc: @(pymPrior) human readable name of the prior
%	    .f: @(pymPrior,a,d) function handle to prior on gamma 
%		(and its derivative) transformed through the Jaccobian
%	    .param: records extra parameters used to create the function handle
%
% Copyright 2012 Pillow lab. All rights reserved.

pymPrior.param = [];
switch lower(priorName)
case 'identity'
    pymPrior.desc = @(pyp) 'identity';
    pymPrior.fgamma = @(pyp, gg) ones(size(gg));
    pymPrior.dfgamma = @(pyp, gg) zeros(size(gg));
    pymPrior.ddfgamma = @(pyp, gg) zeros(size(gg));
case 'triangle'
    pymPrior.desc = @(pyp) 'triangle';
    pymPrior.fgamma = @(pyp, gg) (1 - gg);
    pymPrior.dfgamma = @(pyp, gg) -ones(size(gg));
    pymPrior.ddfgamma = @(pyp, gg) zeros(size(gg));
case 'exponential'
    if nargin < 2; param = 1; end
    assert(param > 0);
    pymPrior.param = param;
    pymPrior.desc = @(pyp) sprintf('exponential [%.2f]', pyp.param);
    pymPrior.fgamma = @(pyp, gg) exp(-gg ./ pyp.param);
    pymPrior.dfgamma = @(pyp, gg) -exp(-gg ./ pyp.param) / pyp.param;
    pymPrior.ddfgamma = @(pyp, gg) exp(-gg ./ pyp.param) / pyp.param.^2;
case {'evan', 'default'}
    if nargin < 2; param = 0.1; end
    assert(param > 0);
    pymPrior.param = param;
    pymPrior.desc = @(pyp) sprintf('evan [%.2f]', pyp.param);
    pymPrior.fgamma = @(pyp, gg) exp(-1 ./(1-gg)/pyp.param);
    pymPrior.dfgamma = @(pyp, gg) -exp(-1 ./(1-gg)/pyp.param)./(1-gg).^2/pyp.param;
    pymPrior.ddfgamma = @(pyp, gg) exp(-1 ./(1-gg)/pyp.param)./(1-gg).^4/pyp.param^2.*(1+2*(gg-1)*pyp.param); 
    case 'pillow'
    pymPrior.desc = @(pyp) 'pillow';
    pymPrior.fgamma = @(pyp, gg) 1 - gg.^2;
    pymPrior.dfgamma = @(pyp, gg) -2 * gg;
    pymPrior.ddfgamma = @(pyp, gg) - 2*ones(size(gg));
case 'dp'
    pymPrior.desc = @(pyp) 'DP';
    pymPrior.f = @(pyp,a,d) a.* trigamma(a+1) .* (d == 0); % TODO
    error('Would you like to implement me? Please?');
otherwise
    error('Unknown prior specification');
end
pymPrior.f = @(pyp,a,d) computePrior(pymPrior, a, d);

end

function [p dp ddp] = computePrior(pymPrior, a, d)
    [~, gg_t, hh_t, iSz, dga, dgd, dha, dhd] = common_prior_Jacobian(a,d);
    p = reshape(pymPrior.fgamma(pymPrior, gg_t), iSz);
    if nargout > 1
	% for optimization purpose only
	dg = pymPrior.dfgamma(pymPrior, gg_t);
	dp(:,1) = dg .* dga;
	dp(:,2) = dg .* dgd;
    end
    if nargout > 2
	% for optimization purpose only
        ddp = zeros(2);
        ddg = pymPrior.ddfgamma(pymPrior, gg_t);
        ddp(1,1) = ddg.*dga.^2 + dg.*(2*gg_t./hh_t.^2.*dha.^2 - gg_t./hh_t.*psi(2,1+a));
        ddp(1,2) = ddg.*dga.*dgd - dg.*(dhd.*dha./hh_t.^2 - (2*gg_t .* dha .* dhd ./ hh_t.^2));
        ddp(2,1) = ddp(1,2);
        ddp(2,2) = ddg.*dgd.^2 + dg.*(psi(2,1-d).*(digamma(1) - digamma(a+1))./(hh_t.^2) + 2*trigamma(1-d).^2.*(digamma(1) - digamma(1+a))./hh_t.^3);
    end
end

function [J, gg_t, hh_t, iSz, dga, dgd, dha, dhd] = common_prior_Jacobian(a,d)
    iSz = size(a);
    a = a(:);  d = d(:);
    [hh_t, gg_t] = convertfromPYparams(a,d);
    if any(isnan(gg_t))
	warning('pymPrior:NaNGamma', 'NaN in gamma');
    end
    
    % compute jacobian
    dha = trigamma(a+1); dhd = trigamma(1-d);
    dga = -(gg_t./hh_t).*dha; dgd = dhd.*(digamma(a+1)-digamma(1))./(hh_t.^2);
    J = dha.*dgd - dhd.*dga; % compute 2x2 determinant
    assert(all(J>0)) % jacobian should always be invertible
end
