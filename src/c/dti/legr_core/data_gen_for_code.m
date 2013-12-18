%% Save files
fh = fopen('idx_dti_test_int_arma.mat','w');
for i = 1:size(idx_dti,1)
    for j = 1:size(idx_dti,2)
        fprintf(fh, '%d ',idx_dti(i,j));
    end
    fprintf(fh,'\n');
end
fclose(fh)