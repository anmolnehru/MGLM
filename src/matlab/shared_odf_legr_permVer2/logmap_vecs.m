function V = logmap_vecs(xi,X)
% LOGMAP maps X on d-sphere to a tangent space of manifold M at xi.
% X is column vectors

V = zeros(size(X));
for j = 1:size(X,2)
    xj = X(:,j);
    V(:,j) = logmap(xi,xj);
end

