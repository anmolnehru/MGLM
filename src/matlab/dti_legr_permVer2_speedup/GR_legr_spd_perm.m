function ErrMx = GR_legr_spd_perm(X, Y,idx_dti, varargin)
%GR_LEGR_SPD_PERM

    ndata =  size(X,2);
    nperms = size(idx_dti,1);
    
    if nargin >=4
        niter = varargin{1};
    else
        niter = 100;
    end
    p = karcher_mean_spd(Y, [], niter);
    logY = logmap_pt2array_spd(p, Y);
    %% Embedding in R6
    Yv = embeddingR6_vecs(p,logY);

    Yhat_perm = zeros([size(Y) nperms]);
    
    %% For each permutation
    for iperm = 1:size(idx_dti,1)
        Xc = X(:,idx_dti(iperm,:));
        L = Yv/Xc; %% speed up ?
        logYv_hat = L*Xc;

        %% Think about algorithm carefully.
        V_hat = invembeddingR6_vecs(p,logYv_hat);
        Yhat_perm(:,:,:,iperm) = expmap_pt2array_spd(p,V_hat);
    end



%% Optimize this part
%% Error calculation
ErrMx = zeros(1,nperms);

for iperm =1:nperms
    
    E = 0 ;
    for i = 1:ndata
        E = E + dist_M_pt2array_spd(Y_hat(:,:,i),Y(:,:,i))^2;
    end
end
