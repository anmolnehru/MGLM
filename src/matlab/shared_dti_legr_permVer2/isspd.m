function T = isspd(mx)
% Check matrices are symmetric positive definite.
T = zeros(size(mx,3),1);
for i=1:size(mx,3)
   T(i) = (sum(eig(mx(:,:,i)) <= 0 ) ==0);
end