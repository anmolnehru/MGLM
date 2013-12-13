function W_new = paralleltranslate(p,V,W)
%PARALLELTRANSLATE sends W from TpM to Texp(p,V)M.
if size(V,2) < size(W,2)
    V = V*ones(1,size(W,2));
end
vnorm = sqrt(diag(V'*V));
scale_factor = diag(V'*W)./vnorm.^2;
scale_factor(isnan(scale_factor)) = 0;
W_orth = W - V*diag(scale_factor);
if size(p,2) == 1
    V_par = p*((-sin(vnorm)).*vnorm)' + V*diag(cos(vnorm));
else
    V_par = p*diag((-sin(vnorm)).*vnorm) + V*diag(cos(vnorm));
end
W_new = V_par*diag(scale_factor) + W_orth;
