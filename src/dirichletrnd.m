function p = dirichletrnd(alphas,K,N)
% p = dirichletrnd(alphas,K,N)
% This function draws N K-vectors from a dirichlet distribution with
% parameters alphas. If alphas is not a vector, it is assumed that all K
% parameters equal the scalar alphas. 
% inputs:
%       alphas : dirichlet parameters; scalar or K-vector
%           K : number of bins
%           N : number of draws to return
% returns: 
%       p     : [K x N] matrix (random draws along columns)
% 
% Evan Archer
% Copyright Pillow Lab 2010-2013. All rights reserved.

    assert( min( alphas(:) ) >= 0, 'Alphas should be >= 0.');
    if any(alphas(:) == 0)
       	warning('PYM:dirchletrnd', 'Some alphas are zero.');
    end

    if(nargin < 3)
        N = 1;
    end
    
    if(nargin < 2 || isempty(K))
        K = length(alphas);
    end
    
    alphas = alphas(:); 
    
    I = ones(K,N);
    
    scalar_alphas = (length(unique(alphas)) == 1);
        
    alphas = bsxfun(@times, I, alphas);
    
    y = gamrnd(alphas, 1);
    
    sy = sum(y);
    
    % if we've sampled an empty vector, 
    % uniformly place a random one someplace in the vector
    nii = find(sy<eps);
    if ~isempty(nii)
        if ~scalar_alphas
           warning('PYM:dirchletrnd:toosmall', 'I drew distributions empty up to numerical precision. In this case, I place all mass in a single bin chosen uniformly at random.') 
        end
        ii0 = round((K-1)*rand(1,length(nii)))+1;
        for kdx = 1:length(nii)
            y(ii0(kdx), nii(kdx))  = 1;
        end
        sy = sum(y);
    end
    
    p = bsxfun(@times, y, 1./sy);
end
