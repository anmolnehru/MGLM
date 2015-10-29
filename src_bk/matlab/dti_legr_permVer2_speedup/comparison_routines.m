%
t1 =clock;
for i =1:10000
    L1  = (Xc'\Yv')';
end
t2 = clock;
elpasedtime1 = etime(t2,t1)

t1 =clock;
for i =1:10000
    L2  = (Yv/Xc);
end
t2 = clock;
elpasedtime2 = etime(t2,t1)

t1 =clock;
for i =1:10000
    [Lx Ux ]  = lu(Xc');
    L3  = (Ux\(Lx\Yv'))';
end
t2 = clock;
elpasedtime3 = etime(t2,t1)

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