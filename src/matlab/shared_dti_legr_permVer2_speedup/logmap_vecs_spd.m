function V = logmap_vecs_spd(X,Y)
% LOGMAP maps X on GL(n)/O(n) to a tangent space of manifold M at xi.
% X is column vectors

V = zeros(size(Y));
if size(X,3) ==1 
    for i = 1:size(Y,3)
        yi = Y(:,:,i);
        V(:,:,i) = logmap_spd(X,yi);
    end
else
    for i = 1:size(X,3)
        xi = X(:,:,i);
        yi = Y(:,:,i);
        V(:,:,i) = logmap_spd(xi,yi);
    end
end

