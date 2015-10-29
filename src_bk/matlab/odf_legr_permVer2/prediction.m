function phat = prediction(p,V,X,expmap)
%prediction predicts values based on one or multiple tangent vectors
ndata = length(X);
phat = zeros(size(p,1),ndata);

for i = 1:ndata
    phat(:,i) = expmap(p,V*X(:,i));
end
