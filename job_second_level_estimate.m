function job = job_second_level_estimate(spmfile)
matlabbatch{1}.spm.stats.fmri_est.spmmat = spmfile;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);