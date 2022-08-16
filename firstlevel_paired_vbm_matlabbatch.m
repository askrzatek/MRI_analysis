function matlabbatch = firstlevel_paired_vbm_matlabbatch(subj,outdir,covars,par)
%%

    skip = [];
    if ~exist ('par','var')
       par = '';
    end

    %% defpar
    defpar.jobname                            = 'spm_spec_firstlevel_paired_VBM';
    defpar.walltime                           = '04:00:00';
    defpar.run                                = 0;

    par = complet_struct(par, defpar);

    %%
    ijob = 1;
    for iS = 1: length(subj)
        iout = 1;
        jobs{ijob}.spm.stats.factorial_design.dir = outdir(iS);
        scans = [subj{iS}(1),subj{iS}(2)];
        jobs{ijob}.spm.stats.factorial_design.des.t2.scans1 = spm_select('expand',scans(1));
        jobs{ijob}.spm.stats.factorial_design.des.t2.scans2 = spm_select('expand',scans(2));
        %     {
    %                                                                    '/home/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mri/smwp1rmv_PARKGAMEII_001_NB_18_07_2018_V1_S3_t1mpr_S256_0_8iso_p2_a.nii,1'
    %                                                                    '/home/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mri/smwp1rmv_PARKGAMEII_001_NB_29_08_2018_V2_S3_t1mpr_S256_0_8iso_p2_a.nii,1'
    %                                                                    };
        jobs{ijob}.spm.stats.factorial_design.des.t2.dept = 0;
        jobs{ijob}.spm.stats.factorial_design.des.t2.variance = 1;
        jobs{ijob}.spm.stats.factorial_design.des.t2.gmsca = 0;
        jobs{ijob}.spm.stats.factorial_design.des.t2.ancova = 0;
        jobs{ijob}.spm.stats.factorial_design.cov(1).c = covars{1,iS};
        jobs{ijob}.spm.stats.factorial_design.cov(1).cname = 'Age';
        jobs{ijob}.spm.stats.factorial_design.cov(1).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(1).iCC = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(2).c = covars{2,iS};
        jobs{ijob}.spm.stats.factorial_design.cov(2).cname = 'Gender';
        jobs{ijob}.spm.stats.factorial_design.cov(2).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(2).iCC = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(3).c = covars{3,iS};
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