function E = fevalGR_mv(p,V,X,P,expmap,logmap)
% Sum of squared geodesic error.
%
ndata = length(X);
P_hat = prediction(p,V,X,expmap);
E = 0 ;

for i = 1:ndata
    E = E + norm(logmap(P_hat(:,i),P(:,i)))^2;
end