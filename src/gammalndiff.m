function [f,checkf] = gammalndiff(x,dx)
% [f,checkf] = digammalndiff(x,dx)
%
% Compute the difference gammaln(x+dx)-gammaln(x), which is accurate when x
% is very large
%
% Three algorithms are possible:
% (1) for x<1e10, gammaln(x+dx)-gammln(x)
% (2) for x>1e10, digamma(x+dx/2)*dx
% (3) for checking accuracy: sum(log(x+(0:dx-1))) 
%     (passed back as optional second argument, to check accuracy).
%

TOL = 1e10; % use approximate formula for values above this

if length(x)==1  
    % scalar x
    if x<TOL
	f = gammaln(x+dx)-gammaln(x);
    else
	f = digamma(x+dx/2)*dx;
    end
else
    % vector x
    f = zeros(size(x));
    ii1 = (x<TOL);
    ii2 = (x>=TOL);
    if length(dx)==1  % scalar dx
	f(ii1) = gammaln(x(ii1)+dx)-gammaln(x(ii1));
	f(ii2) = digamma(x(ii2)+dx/2)*dx;
    else % vector dx
	f(ii1) = gammaln(x(ii1)+dx(ii1))-gammaln(x(ii1));
        if any(ii2)
            f(ii2) = digamma(x(ii2)+dx(ii2)/2).*dx(ii2);
        end
    end
end

% Code for checking accuracy (not really that much slower unless dx>1e3)
if nargout>1
    checkf = sum(log(x+(0:dx-1)));
end
