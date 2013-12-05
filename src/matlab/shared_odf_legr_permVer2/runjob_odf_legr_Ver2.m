% shared job code
load('data');
load('idx_odf.mat');
% exp_name
% X
% mxstack
% mask

mxstack = mxstack_tmp;
disp('Data Loaded.')

results = cell(size(roi_voxel_indices,1),1);
ErrMx = zeros(size(roi_voxel_indices,1),size(idx_odf,1));
tic
for imask = 1:size(roi_voxel_indices,1)
%   fprintf('%d/%d\n ',imask,size(roi_voxel_indices,1));
    Y = mxstack{imask};
    p = [];
    logY = [];
    Yu = [];
    U = [];
    
    for iperm = 1:size(idx_odf,1)
        Xs_perm = Xs(:,idx_odf(iperm,:));
         [p, V, E, Y_hat, logY, Yu, U] = GR_logeuc_Ver2(Xs_perm, Y, p, logY, Yu, U);
         ErrMx(imask,iperm) = E(end);
    end
    
end
toc;
save(strcat([exp_name,'_result_legr_odf']),'ErrMx','roi_voxel_indices');
