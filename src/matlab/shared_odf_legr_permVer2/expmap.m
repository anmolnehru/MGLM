function exp_x_y = expmap(x,y)
% EXPMAP is for exponential map onto d-sphere.
% EXPMAP maps Y in TxM (tangent space of M at X) onto d-sphere.
%normy = norm(y);
%if normy < eps
%    exp_x_y = x;
%else
%    exp_x_y = x*cos(normy)+ y/norm(y)*sin(normy);
%end

exp_x_y = cos(norm(y))*x+sin(norm(y))*y/zero2one(norm(y));