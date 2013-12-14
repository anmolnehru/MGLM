% Comparison
t1 =clock;
for i =1:10000
    [U D] = eig(p);
    g = U*sqrt(D);
    invg1 = inv(g);
end
t2 = clock;
elpasedtime1 = etime(t2,t1)


t1 =clock;
for i =1:10000
    [U D] = eig(p);
    g = U*sqrt(D);
    invg2 = diag(1./sqrt(diag(D)))*U';
end
t2 = clock;
elpasedtime2 = etime(t2,t1)