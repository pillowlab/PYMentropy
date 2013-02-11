function p = dirichletrnd(alpha,K,N)
% p = dirichletrnd(alpha,K,N)
% This function draws N K-vectors from a dirichlet distribution with
% parameters alpha. If alpha is not a vector, it is assumed that all K
% parameters equal the scalar alpha. 
% inputs:
%       alpha : dirichlet parameters; scalar or K-vector
%           K : number of bins
%           N : number of draws to return
% returns: 
%       p     : [K x N] matrix (random draws along columns)
% 
% $Id: dirichletrnd.m 708 2011-11-09 20:47:44Z pillow $
% Evan Archer
% Copyright Pillow Lab 2010/2011. All rights reserved.


    assert( min( alpha(:) ) > 0, 'Alphas should be >0.')

    if(nargin < 3)
        N = 1;
    end
    
    if(nargin < 2 || isempty(K))
        K = length(alpha);
    end
    
    alpha = alpha(:); 
    
    I = ones(K,N);
    
    scalar_alpha = (length(unique(alpha)) == 1);
        
    alpha = bsxfun(@times, I, alpha);
    
    y = gamrnd(alpha, 1);
    
    sy = sum(y);
    
    % if we've sampled an empty vector, 
    % uniformly place a random one someplace in the vector
    nii=find(sy<eps);
    if ~isempty(nii)
        if(~scalar_alpha)
           warning('I drew distributions empty up to numerical precision. I am placing all mass in a single bin chosen uniformly at random; this may not be dirichlet in this case.') 
        end
        ii0= round((K-1)*rand(1,length(nii)))+1;
        for kdx = 1:length(nii)
            y(ii0(kdx), nii(kdx))  = 1;
        end
        sy = sum(y);
    end
    
    p = bsxfun(@times, y, 1./sy);
end
%     
%     
%     for idx = 1:N
%         sy = 0;
% %         while sy == 0
%             y = gamrnd(alpha, I);
%             sy = sum(y);
%             if (sy)
%                 p(:,idx) = y/sy;
%             else
%                 fprintf('*');
%                 p(:,round(K*rand)) = 1;
%             end
% %         end
%     end
    
    
    
    
% % Discard probabilities beneath numerical precision (to avoid problems in sums)     
%     p( p < 1e-15) = 0;
%     p = bsxfun(@times, p,1./sum(p));
%     
