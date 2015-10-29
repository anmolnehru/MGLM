function X = expmap_vecs(xi,V)
% EXPMAP_VECS exponential map V from tangent space of M at xi to X
% on d-sphere.
X = zeros(size(V));
for j = 1:size(V,2)
    X(:,j) = expmap(xi,V(:,j));
end