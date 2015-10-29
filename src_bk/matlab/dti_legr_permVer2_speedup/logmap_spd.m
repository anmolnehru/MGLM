function v = logmap_spd(P,X)
%LOGMAP_SPD maps X on GL(n)/O(n) to a tangent space of manifold M at P.
%
%    The tangent space is a set symmetric matrices(n)

%% Common
if norm(P-X) < 1e-18
    v = zeros(size(P));
    return
end
[U D] = eig(P);
g = U*sqrt(D);

%invg = inv(g);
invg = diag(1./sqrt(diag(D)))*U';

%% For each data
y = invg*X*invg';
[V S] = eig(y);
H = g*V;
v = H*diag(log(diag(S)))*H';


%rtX = sqrtm(X);
%invrtX = inv(rtX);
%v = rtX*logm(invrtX*Y*invrtX)*rtX;
% v = (v+v')/2;