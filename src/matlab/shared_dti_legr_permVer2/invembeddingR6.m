function vnew = invembeddingR6(p, v)
%INVEMBEDDINGR6 converts v in R6 to vnew in TpM.
%    
%    p is a base point. 
%    vnew is a tangent vector.
%    v is a vector in R6.
%    step 1
%    v = [Sxx Sxy Sxz Syy Syz Szz]';    
    w = [1  sqrt(2) sqrt(2)   1 sqrt(2) 1]';
    S = vec2symmx(v./w);
%    step 2    
    sqrtp= sqrtm(p);
    vnew = sqrtp*S*sqrtp;
end