% Comparison
t1 =clock;
for i =1:1000
    d1 = dist_M_spd(p,p2);
end
t2 = clock;
elpasedtime1 = etime(t2,t1)


t1 =clock;
for i =1:1000
    d2 = dist_M_spd_ver2(p,p2);
end
t2 = clock;
elpasedtime2 = etime(t2,t1)