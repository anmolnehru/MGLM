function w_new = paralleltranslateAtoB(a, b, w, logmap)
%PARALLELTRANSLATEATOB translate vector w from TaM to TbM

if size(a,2) < size(b,2)
    a = a*ones(1,size(b,2));
elseif size(a,2) > size(b,2)
    b = b*ones(1,size(a,2));
end

v = zeros(size(a));
for i = 1:size(a,2)
    v(:,i) = logmap(a(:,i),b(:,i));
end

w_new = paralleltranslate(a, v, w);

