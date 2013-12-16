% runjob_dti_legr
% shared job code
%clear;
load('data');
load('idx_dti.mat');


disp('Data Loaded.')
nvoxels = size(mask_job,1);
ErrMx = zeros(nvoxels,size(idx_dti,1));
t1 = clock;
for imask = 1:nvoxels
    ErrMx(imask,:) = GR_legr_spd_perm(Xs, mxstack_job{imask},idx_dti);
end
t2 = clock;
elapsedtime = etime(t2,t1)
save(strcat([exp_name,'_result_legr_dti_Ver2_5']),'ErrMx','mask_job','elapsedtime');

