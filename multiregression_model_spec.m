function jobs = multiregression_model_spec(outdir, cons_a, cons_c, covars, target_regressor, par)
%%

  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_spec_multiregression_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;
defpar.nb_cond                            = 1;
defpar.nb_cons                            = 1;

par = complet_struct(par, defpar);
ijob = 1;
for i = 1 : par.nb_cond % a loop for each target regressor
    for iout = 1 : par.nb_cons % a loop for each contrast // con
%         jobs{ijob}.spm.stats.factorial_design.dir = {'/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/test_multiple_regression/Axial_IL'};
%         jobs{ijob}.spm.stats.factorial_design.des.mreg.scans = {
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_001_NB_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_002_BM_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_007_SD_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_008_JR_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_025_CA_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_039_KM_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_043_PD_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_048_SB_a/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_003_SM_c/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_023_LJ_c/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_028_PC_c/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_033_DD/con_0017.nii,1'
%                                                                     '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/PARKGAMEII_044_CK_c/con_0017.nii,1'
%                                                                     };

        %%
        jobs{ijob}.spm.stats.factorial_design.dir = outdir{i}(iout);
        cons{iout} = [cons_a{1,iout};cons_c{1,iout}];
        %cons{2,iout} = cons_c{iout};
        
        jobs{ijob}.spm.stats.factorial_design.des.mreg.scans = spm_select('expand',cons{iout});
        %%
        %%
        jobs{ijob}.spm.stats.factorial_design.des.mreg.mcov.c = target_regressor.value{i};
        %%
        jobs{ijob}.spm.stats.factorial_design.des.mreg.mcov.cname = target_regressor.name{i};
        jobs{ijob}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
        jobs{ijob}.spm.stats.factorial_design.des.mreg.incint = 1;
        %%
        jobs{ijob}.spm.stats.factorial_design.cov(1).c = covars{1};
        %%
        jobs{ijob}.spm.stats.factorial_design.cov(1).cname = 'Age';
        jobs{ijob}.spm.stats.factorial_design.cov(1).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(1).iCC = 1;
        %%
        jobs{ijob}.spm.stats.factorial_design.cov(2).c = covars{2};
        %%
        jobs{ijob}.spm.stats.factorial_design.cov(2).cname = 'Gender';
        jobs{ijob}.spm.stats.factorial_design.cov(2).iCFI = 1;
        jobs{ijob}.spm.stats.factorial_design.cov(2).iCC = 1;
        jobs{ijob}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        jobs{ijob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        jobs{ijob}.spm.stats.factorial_design.masking.im = 1;
        jobs{ijob}.spm.stats.factorial_design.masking.em = {''};
        jobs{ijob}.spm.stats.factorial_design.globalc.g_omit = 1;
        jobs{ijob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        jobs{ijob}.spm.stats.factorial_design.globalm.glonorm = 1;
        spm('defaults','FMRI')
        ijob = ijob +1;
    %    job = do_cmd_sge(jobs{ijob},par);
    end
end

[ jobs ] = job_ending_rountines( jobs, skip, par );


end