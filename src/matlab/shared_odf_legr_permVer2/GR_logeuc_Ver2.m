function [p, V, E, Y_hat, logY, Yu, U] = GR_logeuc_Ver2(X, Y, p, logY, Yu, U, varargin)
%GR_LOGEUC is logeuclidean geodesic regression. 
%    X input variable in Euclidean space.
%    Y output variable on Sphere
%    p is a base point for tangent space.
%    V is a set of tangent vectors. Column vectors
%    E is geodesic error.
%    Y_hat is prediction
%    U orthogonal bases of the tangent space

if nargin >=7
    niter = varargin{1};
else
    niter = 100;
end

[ ndim ndata] = size(X);

if isempty(p)
    p = karcher_mean(Y, ones(ndata,1)/ndata, niter);
end
if isempty(logY)
    logY = logmap_vecs(p, Y);
end

% Linear transform

%Xc = X - repmat(mean(X,2),1,ndata);

if isempty(Yu)
    % Get orthogonal bases
    U = null(ones(size(p,1),1)*p');
    Yu = U'*logY; %logY is represented by U
end
% Yu = L*X
L = Yu/X;
V = U*L;
logY_hat = V*X;
Y_hat = expmap_vecs(p,logY_hat);
E = geodesic_squred_error(Y, Y_hat, @logmap);