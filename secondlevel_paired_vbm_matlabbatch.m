function matlabbatch = secondlevel_paired_vbm_matlabbatch(groups,outdirs,covars,par)
%%

  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_spec_secondlevel_RS_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;

par = complet_struct(par, defpar);

%%
ijob = 1;
for igroup = 1 : length(groups)
    iout = 1;
%icon
    jobs{ijob}.spm.stats.factorial_design.dir = outdirs{igroup}(iout);

    for ipatient = 1 : length(groups{igroup})
%             jobs{ijob}.spm.stats.factorial_design.des.pt.pair(ipatient).scans = {
%                                                                           '/home/anna.skrzatek/data/nifti_test/firstlevel_RS/PARKGAMEII_001_NB_a/Caudate_L_V1/con_0001.nii,1'
%                                                                           '/home/anna.skrzatek/data/nifti_test/firstlevel_RS/PARKGAMEII_001_NB_a/Caudate_L_V2/con_0001.nii,1'
%                                                                           };
        scans = [groups{igroup}{ipatient}(1),groups{igroup}{ipatient}(2)];
        jobs{ijob}.spm.stats.factorial_design.des.pt.pair(ipatient).scans = spm_select('expand',scans(:));
    end
        jobs{ijob}.spm.stats.factorial_design.des.pt.gmsca = 0;
        jobs{ijob}.spm.stats.factorial_design.des.pt.ancova = 0;
        jobs{ijob}.spm.stats.factorial_design.cov(1).c = covars{1,igroup};
        jobs{ijob}.spm.stats.factorial_design.cov(1).cname = 'Age';
        jobs{ijob}.spm.stats.factorial_design.cov(1).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(1).iCC = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(2).c = covars{2,igroup};
        jobs{ijob}.spm.stats.factorial_design.cov(2).cname = 'Gender';
        jobs{ijob}.spm.stats.factorial_design.cov(2).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(2).iCC = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(3).c = covars{3,igroup};
        jobs{ijob}.spm.stats.factorial_design.cov(3).cname = 'TIV';
        jobs{ijob}.spm.stats.factorial_design.cov(3).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(3).iCC = 1;

%         jobs{ijob}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        jobs{ijob}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        jobs{ijob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        jobs{ijob}.spm.stats.factorial_design.masking.im = 1;
        jobs{ijob}.spm.stats.factorial_design.masking.em = {''};
        jobs{ijob}.spm.stats.factorial_design.globalc.g_omit = 1;
        jobs{ijob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        jobs{ijob}.spm.stats.factorial_design.globalm.glonorm = 1;

        spm('defaults','FMRI')

        ijob = ijob +1;
        iout = iout +1;
    end
    

[ jobs ] = job_ending_rountines( jobs, skip, par );

end