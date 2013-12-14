function v = logmap(xi,xj)
% LOGMAP maps xj on d-sphere to a tangent space of manifold M at xi.
%w = xj-xi'*xj*xi;
%if norm(w) < 1e-10
%    v = zeros(size(xi));
%    return 
%end
%w = w/norm(w);
%v = w*acos(xi'*xj);

v = (xj-xi'*xj*xi)/zero2one(sqrt(1-(xi'*xj)^2))*zero2one(acos(xi'*xj));