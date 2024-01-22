function matlabbatch = ANOVA_RS_matlabbatch(A_groups,A_outdir,covars,par)

%%

  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_spec_secondlevel_RS_ANOVA_auto';
defpar.walltime                           = '04:00:00';
defpar.run                                = 0;
defpar.con_auto                           = 1;

par = complet_struct(par, defpar);

%%
ijob = 1;
for iroi = 1 : length(A_outdir)
    jobs{ijob}.spm.stats.factorial_design.dir = A_outdir(iroi); %{'/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/resliced_RS_ANOVA_V1V2/Cereb6_L_V2_V1'};
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(1).name = 'SESSION';
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(2).name = 'GROUP';
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
    jobs{ijob}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                        1];
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(1).scans =   spm_select('expand',A_groups{1,1}{iroi});
%     {
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_001_NB_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_002_BM_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_007_SD_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_008_JR_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_025_CA_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_039_KM_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_043_PD_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_048_SB_a/Cereb6_L_V1/con_0001.nii,1'
%                                                                        };
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(2).levels = [1
                                                                        2];
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(2).scans =   spm_select('expand',A_groups{2,1}{iroi});
%     {
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_003_SM_c/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_023_LJ_c/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_028_PC_c/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_033_DD/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_044_CK_c/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_047_BF_c/Cereb6_L_V1/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_052_HJ_c/Cereb6_L_V1/con_0001.nii,1'
%                                                                        };
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(3).levels = [2
                                                                        1];
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(3).scans =   spm_select('expand',A_groups{1,2}{iroi});
%     {
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_001_NB_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_002_BM_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_007_SD_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_008_JR_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_025_CA_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_039_KM_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_043_PD_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_048_SB_a/Cereb6_L_V2/con_0001.nii,1'
%                                                                        };
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                        2];
    jobs{ijob}.spm.stats.factorial_design.des.fd.icell(4).scans =   spm_select('expand',A_groups{2,2}{iroi});
%     {
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_003_SM_c/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_023_LJ_c/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_028_PC_c/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_033_DD/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_044_CK_c/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_047_BF_c/Cereb6_L_V2/con_0001.nii,1'
%                                                                        '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/firstlevel_RS/PARKGAMEII_052_HJ_c/Cereb6_L_V2/con_0001.nii,1'
%                                                                        };
    jobs{ijob}.spm.stats.factorial_design.des.fd.contrasts = par.con_auto;
    %%
    jobs{ijob}.spm.stats.factorial_design.cov(1).c = vertcat(covars{1,:},covars{1,:});
    %%
    jobs{ijob}.spm.stats.factorial_design.cov(1).cname = 'Age';
    jobs{ijob}.spm.stats.factorial_design.cov(1).iCFI = 1;
    jobs{ijob}.spm.stats.factorial_design.cov(1).iCC = 1;
    %%
    jobs{ijob}.spm.stats.factorial_design.cov(2).c = vertcat(covars{2,:},covars{2,:});
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
    
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

end