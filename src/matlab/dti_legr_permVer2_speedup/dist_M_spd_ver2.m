function d = dist_M_spd_ver2(X,Y)  
    [U D ] = eig(X);
    invg = U*diag(1./sqrt(diag(D)))*U';
    [U D ] = eig(invg*Y*invg);
    d = sqrt(sum(diag(log(D)).^2));
    
end