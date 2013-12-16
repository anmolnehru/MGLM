% runjob_dti_legr
% shared job code
%clear;
load('data');
load('idx_dti_test_100.mat')
%load('idx_dti_test_100.mat');

% exp_name
% X
% mxstack
% mask

mxstack = mxstack_job;
disp('Data Loaded.')

nvoxels = size(mask_job,1);
ErrMx = zeros(nvoxels,size(idx_dti,1));
t1 = clock;
for imask = 1:nvoxels
%    tic
%    fprintf('%d/%d\n ',imask,nvoxels);
    Y = mxstack{imask};
    p =[];
    logY = [];
    Yv = [];
    for iperm = 1:size(idx_dti,1)
%        iperm
        Xs_perm = Xs(:,idx_dti(iperm,:));
        [p V E Y_hat logY Yv] = GR_legr_spdVer2(Xs_perm,Y,p,logY, Yv);
        ErrMx(imask,iperm) = E(end);
    end
%    toc;
end
t2 = clock;
save(strcat([exp_name,'_result_legr_dti_Ver2']),'ErrMx','mask_job');
elapsedtime = etime(t2,t1)