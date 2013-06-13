function [Hbls, Hvar, SummaryStr] = computeH_PYM(mm, icts, param, verbose)
% [Hbls, Hvar] = computeH_PYM(mm, icts, param, verbose)
% Compute PYM estimate using either adaptive quadrature or an adaptive grid.
% The area of integration is chosen by finding the MAP (peak of the posterior)
% and expanding on the directions of the eigenvector of the estimated Hessian.
%
% INPUT:
%       mm - multiplicites
%     icts - counts for each multiplicty
%    param - (optional) parameter structure
%	  param.pymPrior: (optional) This is a structure obtained from pymPriorFactory.
%         param.returnEB: (optional) If true, return the empirical bayes
%                       estimate of entropy. That is, return the posterior
%                       mean of entropy given the parameters (alpha,d) at the MAP.
%	  param.verbose: (optional) if true, outputs debugging info
%	  param.nagrid: (optional) Integer > 0. Specify number of gridpoints
%	                to use in integral over alpha. 
%	  param.ndgrid: (optional) Integer > 0. Specify number of gridpoints
%	                to use in integral over d. 
%	  param.method: (optional) Specify integration method.
%                       = 1:  Use Chebyshev spectral integration
%                       = 2:  Use Gauss-Legendre quadrature
%  verbose -  (optional) if true, outputs debugging info
%
% OUTPUT:
%  Hbls      - posterior mean estimate of entropy
%	     (Inf or NaN is returned in the pathological case that the data
%	     have no coincidences, or when only a single symbol is observed.)
%  Hvar      - posterior variance
% SummaryStr - a structure containing variables relevant to computation
%
% Outline of algorithm:
%   - Find the MAP estimate of (d, alpha) given data and prior
%   - Estimate the Hessian of the posterior over (d, alpha)
%   - IF the hessian indicates the posterior is concentrated, 
%         * Estimate the region of integration to be a rectangular region 
%           around the MAP estimate and extends proportional to the Hessian.
%         * Use Gauss-Legendre quadrature to integrate.
%         * If param.nagrid and param.ndgrid are not specified, by default:
%                       nagrid = ndgrid = 25.
%   - ELSE the posterior is not well-concentrated
%         * Use Chebyshev spectral integration method (Boyd 1986) for
%           integration
%         * If param.nagrid and param.ndgrid are not specified, by default:
%                       nagrid = 250
%                       ndgrid = 25
%   - Use selected integration method to compute the following:
%	- Z: normalizer
%	- Hbls: expected entropy given (d, alpha)
%	- Hvar: variance of entropy given (d, alpha)
%
% See also: pymPriorFactory
% Requires: optimization toolbox
%
% Copyright 2012-2013 Pillow lab. All rights reserved.

if nargin < 4; verbose = false; end
if nargin < 3; param = []; end
if ~isfield(param, 'pymPrior'); param.pymPrior = pymPriorFactory('default'); end
if isfield(param, 'verbose'); verbose = param.verbose; end
if isfield(param, 'nagrid'); Na = param.nagrid; else Na = []; end
if isfield(param, 'ndgrid'); Nd = param.ndgrid; else Nd = []; end
% NOTE: parameter 'method' is handled below, before the integral weight computation. 

if ~any(icts > 1) % no coincidence, PYM estimator is infinite
    if verbose; fprintf('PYM:No coincidence. Returning infinity\n'); end
    Hbls = Inf; Hvar = Inf;
    return;
end

minAlpha = 10 * eps; % alpha == 0 is very bad for gamma parameterization

% Compute a statistic that relates to the amount of data relative to # of bins
N = dot(mm,icts); K = sum(mm);
if verbose; fprintf('K/N = %g/%g = %g\n', K, N, K/N); end
if K == 1
    warning('PYM:OneSymbol', 'Pathological case of only one symbol');
    Hbls = NaN; Hvar = Inf;
    return
end

%% First, find the MAP estimate of alpha and d using fmincon
% Note that the posterior evidence p(alpha, d|data) concentrates as more
% data is given, and the precision of MAP needs to be increased in such cases

% cost function: negative log posterior over (alpha, d), given data and prior
nlpy = @(x) nlogPostPYoccupancy(x(1),x(2),mm,icts,param.pymPrior);
precisionMAP = 1e-8; % we will increase the precision if needed
mpt = [1; 0.01]; % alpha = 1, d = 0.01 as initial point
% linear constraints for optimization boundaries
if isfield(param, 'allowNegativeAlpha') && param.allowNegativeAlpha
    A = [-1 -1; 0 1; 0 -1]; b = [0;1;0]; % interval constraints (alpha > -d)
else
    A = [-1  0; 0 1; 0 -1]; b = [0;1;0]; % alpha > 0 version
end

options = optimset('Algorithm', 'interior-point', 'Display', 'off', ...
    'GradObj', 'on', 'Hessian', 'user-supplied', 'HessFcn', @hessianfcn, 'TolX', precisionMAP, 'TolFun', 1e-10);
[mpt, fval, ~, ~, ~, ~, hessian] = fmincon(@(x) nlpy(x), mpt(:), A, b,[],[],[],[],[],options);

if verbose
    fprintf('MAP (alpha = %g, d = %g) [fval = %g]\n', mpt(1), mpt(2), fval);
end

if isfield(param, 'returnEB') && param.returnEB
    if verbose; fprintf('Returning empirical bayes estimate.'); end
    [Hbls,Hvar] = computeHpy(mm,icts,mpt(1),mpt(2));
    return
end

fval = -fval;

%% Verify numerical Hessian
if verbose; for idx = 1:2; for jdx = idx:2
    fprintf('(%d,%d)', idx,jdx)
    HessCheck_Elts(@(x) nlpy([x(1) x(2)]), [idx jdx], mpt(:));
end; end; end

if(nargout > 2) % if we are going to pass back summary structure
    SummaryStr.fval = fval;
end

if any(isnan(hessian(:))) || any(isinf(hessian(:)))
    warning('PYM:naninfHessian', 'Hessian contains nan or inf values. Data may be very unlikely under the chosen prior.')
    p = true;
else
    [~, p] = chol(hessian); % use chol to see if Hessian is positive definite
end

if mpt(1) < minAlpha
    warning('PYM:smallAlphaMAP', 'MAP alpha is very small. Reparameterization may be numerically unstable');
end

if p || rank(hessian) < 2  % if we failed to get a good Hessian
    warning('PYM:hessian', 'Hessian not positive definite. Computing integral on full semi-infinite interval. This will be slower.');
    method = 1; % TODO: rather than indexing the integration methods by integers, add variables.
    dl = eps; du = 1-eps;
else % use the Hessian to make the region of integration
    method = 2; % Posterior appears to be adequately concentrated. Approximate integral on rectangular region.
    invHessian = inv(hessian); % CAUTION: small eigenvalues can explode!
    invHessian = 6*sqrt(invHessian); % 6 times "std"
    %% Using the invHessian, determine the region of integration
    al = max(mpt(1) - invHessian(1,1), eps);
    au = mpt(1) + invHessian(1,1);
    dl = max(mpt(2) - invHessian(2,2), eps);
    du = min(1-eps,max(mpt(2) + invHessian(2,2), eps));
end

if verbose; fprintf('[al,au,dl,du] = [%g,%g,%g,%g]\n',al,au,dl,du); end
if verbose; fprintf('== Hessian ==\n'); disp(hessian); end

if(nargout > 2) % if we are going to pass back summary structure
    SummaryStr.al = al; SummaryStr.au = au;
    SummaryStr.dl = dl; SummaryStr.du = du;
    SummaryStr.lpy = @(x) exp(-nlogPostPYoccupancy(x(1),x(2),mm,icts,param.pymPrior) - fval);
end

%% Compute weights for integral

% If method field is specified, use the specified method!
if isfield(param, 'method'); method = param.method; end;

switch method % do the semi-infinite integral
    case 1 %% TODO: Use this method when the hessian fails and we cannot choose a reasonable region over which to integrate
        if isempty(Na); Na = 250; end
        if isempty(Nd); Nd = 25; end
	% Chebyshev spectral integration method (Boyd 1986)
        aw = zeros(1,Na);
        %
        tt = pi*(1:Na)/(Na+1);
        jj = 1:Na;
        ax = cot(tt/2).^2; % the gridpoints at which we evaluate the integrand
        L = 1; % This parameter controls the amount of weight the numerical rule gives to the tail. It is chosen by heuristic.
        ax = L*ax;
        z0 = (1-cos(pi*jj))./jj;
        for idx = 1:Na
            aw(idx) =  sum(sin(idx*tt).*z0);
        end
        aw = L*( 2*sin(tt)./(1-cos(tt)).^2 ) .* aw .* ( 2/(Na+1) );
        aw = aw'; % make aw into a column vector
        %         [dx dw Nd]   = gq100(dl, du);
        [dx dw] = lgwt(Nd, dl, du);
    case 2 % do gauss quadrature on the finite interval
        [ax aw Na]  = gq100(al, au, Na); % If param.nagrid is not specified, input Na will be []. Default Na is specified in gq100.
        [dx dw Nd]  = gq100(dl, du, Nd);
    otherwise
        error('Unrecognized integration method: %d. Method parameter must be either 1 or 2; please see documentation.', method);
end

if verbose 
    fprintf('Integrating alpha with %d gridpoints.\n Integrating d with %d gridpoints.\n', Na,Nd)
end

% likelihood evaluated at the MAP
likMapVal = logliPYoccupancy(mpt(1), mpt(2), mm, icts);

%% Compute function at gridpoints
if Nd * Na < 1e4
    [aa dd] = meshgrid(ax,dx);
    loglik = logliPYoccupancy(aa(:),dd(:),mm,icts);
    lik = exp(loglik - max(loglik));
    prior = param.pymPrior.f(param.pymPrior,aa(:),dd(:));
    [mc vc] = condH(aa(:), dd(:), mm, icts);
    A = bsxfun(@times, lik.*prior, vec(dw * aw'));
    Z = sum(A);
    Hbls = sum(A.*mc);
    Hvar = (Z*sum(sum(A.*vc)) - Hbls.^2)/Z.^2;
    Hbls = Hbls/Z;
else
    % On a very large grid, it will be more efficient to compute the marginals
    % rather than store the entire grid at once.
    a_Z = zeros(1,Na);
    a_Hbls = zeros(1,Na);
    a_Hmom2 = zeros(1,Na);
    for adx = 1:Na
        aa = ax(adx)*ones(size(dx));
	if verbose; fprintf('%0.1f\r', adx/Na); end
        loglik = logliPYoccupancy(aa,dx,mm,icts);
        lik = exp(loglik - max(loglik));
        prior = param.pymPrior.f(param.pymPrior,aa,dx);
        [mc vc] = condH(aa, dx, mm, icts);
        %% Compute integral
        qq = dw.*lik.*prior;
        a_Z(adx) = sum(qq);
        a_Hbls(adx) = sum(qq.*mc);
        a_Hmom2(adx) = sum(qq.*vc);
    end
    Z = dot(a_Z,aw);
    Hbls = dot(a_Hbls, aw);
    Hvar = (Z*dot(a_Hmom2,aw) - Hbls^2)/Z^2;
    Hbls = Hbls/Z;
end

if sum(mm(icts>1)) < 2
    if verbose; fprintf('Less than 2 coincidences: infinit var [originally: %g]\n', Hvar); end
    Hvar = Inf;
elseif nargout > 1 && Hvar < 0
    warning('PYM:negativeVariance', 'Negative variance obtained due to suspected numerical errors');
    save([mfilename '_' datestr(now, 30) '_error']);
    if verbose; beep; keyboard; end
end

%%%%%%%%
%% Utility functions below
%%%%%%%%

%% Code for computing the conditional mean and variance
% given concentration and discount parameters
    function [m v] = condH(a,d,mm,icts)
        [m,v] = computeHpy(mm,icts,a(:),d(:));
        v = v + m.^2;
        v = reshape(v, size(a));
        m = reshape(m, size(a));
    end

    function hessian = hessianfcn(x, lambda)
        [~,~,hessian] = nlpy(x);
    end
end% end computeH_PYM_quad

%% code for computing likelihood and its derivatives.
function [nlogp ndlogp nddlogp] = nlogPostPYoccupancy(a,d,mm,icts,pymPrior)
    if a < 0 || d < 0 || d >= 1; nlogp = inf; ndlogp = -[a,d]; return; end
    if nargout > 2
	[logp dlogp, ddlogp] = logliPYoccupancy(a,d,mm,icts);
	[prior, dprior, ddprior] = pymPrior.f(pymPrior, a, d);
	ndlogp = -dlogp - dprior ./ prior;
	nddlogp = -ddlogp - (prior*ddprior - dprior'*dprior) ./prior^2;
    elseif nargout > 1
	[logp dlogp] = logliPYoccupancy(a,d,mm,icts);
	[prior, dprior] = pymPrior.f(pymPrior, a, d);
	ndlogp = -dlogp - dprior ./ prior;
    else
	logp = logliPYoccupancy(a,d,mm,icts);
	prior = pymPrior.f(pymPrior, a, d);
    end
    nlogp = -logp - log(prior);

    nlogp(a < 0) = inf;
    nlogp(d < 0 | d >= 1) = inf;
    %fprintf('(a,d) = [%g, %g], nlppyo = [%g]\n', a, d, nlogp);
end

function [x w N] = gq100(a,b, ngrid)
    %% This code adapted from the code 'lgwt' by Greg von Winckel.
    if nargin < 3 || isempty(ngrid); N = 25; else N = ngrid; end
    persistent N_persistent y_persistent Lp_persistent

    if isempty(N_persistent) || N_persistent ~= N
	[~,~,y,Lp] = lgwt(N,-1,1);
	% Compute the weights

	N_persistent = N;
	y_persistent = y;
	Lp_persistent = Lp;
    end
    
    N = N_persistent;
    y = y_persistent;
    Lp = Lp_persistent;

    N1 = N; N2 = N+1;
    x = (a*(1-y) + b*(1+y))/2;
    w = (b-a) ./ ((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
