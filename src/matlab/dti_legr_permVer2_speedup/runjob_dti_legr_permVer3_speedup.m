% runjob_dti_legr
% shared job code
%clear;
%load('data');
%load('idx_dti_test.mat');
% exp_name
% X
% mxstack
% mask

mxstack = mxstack_job;
disp('Data Loaded.')
t1 = clock;
nvoxels = size(mask_job,1);
ErrMx = zeros(nvoxels,size(idx_dti,1));

for imask = 1:nvoxels
    ErrMx(imask,:) = GR_legr_spd_perm(X, Y,idx_dti);
end

save(strcat([exp_name,'_result_legr_dti_Ver2']),'ErrMx','mask_job');
t2 = clock;
elapsedtime = etime(t2,t1)