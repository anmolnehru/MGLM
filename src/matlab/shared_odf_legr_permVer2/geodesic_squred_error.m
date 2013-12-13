function gsr = geodesic_squred_error(X, X_hat, logfunc)
[ndim, ndata] = size(X);
tanvec = zeros(ndim, ndata);
for idata = 1:ndata
    tanvec(:,idata) = logfunc(X(:,idata), X_hat(:,idata));
end
%tanvec
gsr = trace(tanvec'*tanvec);