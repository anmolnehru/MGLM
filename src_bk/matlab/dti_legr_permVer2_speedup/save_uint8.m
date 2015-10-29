function save_uint8(fname,mx)
fid = fopen(fname,'w');

for irow = 1:size(mx,1)
    fprintf(fid, '%d ', mx(irow,:));
    fprintf(fid,'\n');
end

fclose(fid);