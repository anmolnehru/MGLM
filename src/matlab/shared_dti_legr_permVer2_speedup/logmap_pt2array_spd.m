function V = logmap_pt2array_spd(p,X)
%LOGMAP_PT2ARRAY_SPD maps X from GL(n)/O(n) to TpM at xi.



[U D] = eig(p);
g = U*sqrt(D);
invg = inv(g);


%% For each data
for i = 1:size(X,3)
    if norm(P-X) < 1e-18
        v = zeros(size(P));
        return 
    end
    y = invg*X*invg';
    [V S] = eig(y);
    H = g*V;
    v = H*diag(log(diag(S)))*H';
end