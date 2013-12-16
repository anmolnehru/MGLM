load data
save('X_arma.mat','Xs','-ascii');

ui_mask_job = uint8(mask_job);
save_uint8('mask_job_arma.mat', ui_mask_job)


%%
Dxx = 1;Dxy = 2;Dxz = 3;Dyy = 5;Dyz = 6;Dzz = 9;
ndim =6;
Dset = [Dxx, Dxy Dxz Dyy Dyz Dzz]';
nROIs = length(mxstack_job);
nrows = ndim*nROIs;
ncols = size(mxstack_job{1},3);
Ys = zeros(nrows,ncols);

for iROI=1:nROIs
    tmp = reshape(mxstack_job{iROI},9,[]);
    Ys(((iROI-1)*ndim+1):(iROI*ndim) ,:) = tmp(Dset,:);
end

save('Ys_arma.mat','Ys','-ascii');
load idx_dti_test

%%
fid = fopen('idx_dti_test_int_arma.mat','w');
for irow = 1:size(uidx_dti,1)
    fprintf(fid, '%d ', uidx_dti(irow,:));
    fprintf(fid,'\n');
end
fclose(fid);
%%
for irow = size(uidx_dti,1)
    for icol = size(uidx_dti,2)
        fprintf(fid,[repmat('%g ',1,size(uidx_dti,2)),'\n'], uidx_dti);
    end
end
fclose(fid);