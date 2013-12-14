function d = dist_M_spd(X,Y)
    V = logmap_spd(X,Y);
    d = sqrt(innerprod_TpM_spd(V,V,X));
    
    %% This version will be faster.
    invrtX = inv(sqrtm(X));
    tmp = logm(invrtX*Y*invrtX);
    d = sqrt(trace(tmp*tmp));
    
    
    %%
%     [U D ] = eig(X);
%     invg = diag(1./sqrt(diag(D)))*U';
%     
%     [U D ] = eig(invg*Y*invg);
%     d = sum(diag(log(D)).^2);
    
end