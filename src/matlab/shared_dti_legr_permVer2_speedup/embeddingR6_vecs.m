function Vnew = embeddingR6_vecs(p,V)
%EMBEDDINGR6_VECS embeds V in TpM onto R6.
%    
%    p is a base point. V is tangent vectors.
%    V is 3-by-3-n. n is the number of matrices.

if length(size(V)) == 2
    Vnew = embeddingR6(p, V);
    return;
elseif length(size(V)) == 3
    nmx = size(V,3);
    Vnew = zeros(6,nmx);
    for i=1:nmx
        Vnew(:,i) = embeddingR6(p, V(:,:,i));
    end
else
    error('V is wrong input');
end
