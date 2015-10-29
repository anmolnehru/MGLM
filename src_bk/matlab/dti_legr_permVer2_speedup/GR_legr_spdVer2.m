function [p, V, E, Y_hat logY Yv] = GR_legr_spdVer2(X, Y, p, logY, Yv, varargin)
%Geodesic regression for 3d circle and 1d tangent space
%    X is input variables
%    Y is values on manifolds. Symmetric positive matrices in Y(:,:,1) ..
%    Y(:,:,N)
%    V is a set of tangent vectors. V(:,:,1) ... V(:,:,ndimX).
% This is for symmetric positive matrices.

    ndimX = size(X,1);
    ndimY = size(Y,1);
    ndata =  size(X,2);
    
    if ndata ~= size(Y,3)
        error('Different number of input variables and response variables')
    end
    
    if nargin >=6
        niter = varargin{1};
    else
        niter = 100;
    end
    if isempty(p);
        p = karcher_mean_spd(Y, [], niter);
    end
    if isempty(logY)
        logY = logmap_pt2array_spd(p, Y);
    end  
    % xx xy xz yy yz zz
    Xc = X - repmat(mean(X,2),1,ndata);
    %% Embedding in R6
    if isempty(Yv)
        Yv = embeddingR6_vecs(p,logY);
    end
    L = Yv/Xc;
    logYv_hat = L*Xc;
    
    %% Think about algorithm carefully.
    
    V_hat = invembeddingR6_vecs(p,logYv_hat);
    V = invembeddingR6_vecs(p,L);
    
    Y_hat = expmap_pt2array_spd(p,V_hat);

    %% Assume no invalid prediction.
%      for i = 1:ndata
%          if ~isspd(Y_hat(:,:,i))
%              Y_hat(:,:,i) = proj_M_spd(Y_hat(:,:,i));
%              disp('projection')
%          end
%      end


    E = 0 ;
    for i = 1:ndata
        E = E + dist_M_spd_ver2(Y_hat(:,:,i),Y(:,:,i))^2;
    end
end
