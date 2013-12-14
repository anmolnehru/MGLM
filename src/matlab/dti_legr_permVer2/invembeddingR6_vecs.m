function Vnew = invembeddingR6_vecs(p, V)
%INVEMBEDDINGR6 converts v in R6 to vnew in TpM.
%    
%    p is a base point. 
%    vnew is a tangent vector.
%    v is a vector in R6.
%    step 1
%    v = [Sxx Sxy Sxz Syy Syz Szz]';    

if size(V,2) == 1
    Vnew = invembeddingR6(p, V);
    return;
elseif length(size(V)) == 2
    nmx = size(V, 2);
    Vnew = zeros(3,3, nmx);
    for i=1:nmx
        Vnew(:,:,i) = invembeddingR6(p, V(:,i));
    end
else
    error('V is wrong input');
end