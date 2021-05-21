function matlabbatch = secondlevel_RS_matlabbatch(groups,outdirs,par)
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
    for icon = 1 : length(groups{igroup})
        % jobs{ijob}.spm.stats.factorial_design.dir = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/secondlevel_RS_mar21/RS_PARKGAME_c/SMA_R_c/V2'};
        % jobs{ijob}.spm.stats.factorial_design.des.t1.scans = {
        %                                                           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_10_26_PARKGAMEII_003_SM_26_10_2018_V2_c/S05_RS/model/Modele_VOI__Supp_Motor_Area_R/con_0001.nii,1'
        %                                                           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2019_07_12_PARKGAMEII_023_LJ_12_07_2019_V2_c/S06_RS/model/Modele_VOI__Supp_Motor_Area_R/con_0001.nii,1'
        %                                                           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2019_07_31_PARKGAMEII_028_PC_31_07_2019_V2_c/S06_RS/model/Modele_VOI__Supp_Motor_Area_R/con_0001.nii,1'
        %                                                           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2019_11_06_PARKGAMEII_033_DD_06_11_2019_V2_c/S06_RS/model/Modele_VOI__Supp_Motor_Area_R/con_0001.nii,1'
        %                                                           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2020_07_27_PARKGAMEII_044_CK_27_07_2020_V2_c/S06_RS/model/Modele_VOI__Supp_Motor_Area_R/con_0001.nii,1'
        %                                                           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2020_11_18_PARKGAMEII_047_BF_18_11_2020_V2_c/S06_RS/model/Modele_VOI__Supp_Motor_Area_R/con_0001.nii,1'
        %                                                           };
        jobs{ijob}.spm.stats.factorial_design.dir = outdirs{igroup}(icon);
        jobs{ijob}.spm.stats.factorial_design.des.t1.scans = spm_select('expand',groups{igroup}{icon});
        jobs{ijob}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        jobs{ijob}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        jobs{ijob}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        jobs{ijob}.spm.stats.factorial_design.masking.im = 1;
        jobs{ijob}.spm.stats.factorial_design.masking.em = {''};
        jobs{ijob}.spm.stats.factorial_design.globalc.g_omit = 1;
        jobs{ijob}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        jobs{ijob}.spm.stats.factorial_design.globalm.glonorm = 1;

        spm('defaults','FMRI')

        ijob = ijob +1;
    end
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

end