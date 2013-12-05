function [xbar xbar_history] = karcher_mean(X, W, niter)
% KARCHER_MEAN calculates weighted means on d-sphere.
% W is weight
% X is column vectors
% niter is max iteration
xbar = X(:,1);
if isempty(W)
   W = ones(size(X,2),1)/size(X,2);
else
   W = W/norm(W,1);
end

% r = 1; % learning rate 
xbar_history = xbar;
for iter = 1:niter
    phi = logmap_vecs(xbar,X)*W;
    xbar = expmap(xbar, phi);
    xbar_history = [xbar_history, xbar];
    if phi < 1e-10
        break
    end
end