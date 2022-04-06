%%
function jobs = vbm_preproc_long (par)

%%

  skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.jobname                            = 'spm_vbm_long_preproc';
defpar.walltime                           = '04:00:00';
defpar.run                                = 1;

par = complet_struct(par, defpar);

%%%%%%%%%%%%%%%%%
jobs{1}.spm.tools.cat.cat_simple_long.datalong.subjects = {
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_001_NB_18_07_2018_V1_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_001_NB_29_08_2018_V2_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_002_BM_25_07_2018_V1_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_002_BM_05_09_2018_V2_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_003_SM_19_09_2018_V1_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_003_SM_26_10_2018_V2_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_007_SD_09_01_2019_V1_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_007_SD_15_02_2019_V2_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_008_JR_18_01_2019_V1_S4_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_008_JR_18_01_2019_V2_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_023_LJ_24_05_2019_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_023_LJ_12_07_2019_S02_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_025_CA_29_05_2019_V1_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_025_CA_10_07_2019_S02_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_028_PC_21_06_2019_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_028_PC_S02_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_033_DD_11092019_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_033_DD_S02_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_039_KM_S01_S4_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_039_KM_S02_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_043_PD_S01_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_043_PD_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_044_CK_S01_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_044_CK_27_07_2020_V2_S3_t1mpr_S256_0_8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_047_BF_S1_S3_t1mpr_S256_0.8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_047_BF_S4_t1mpr_S256_0.8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_048_SB_S3_t1mpr_S256_0.8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_048_SB_06_01_2021_V2_S3_t1mpr_S256_0.8iso_p2.nii,1'
                                                                  }
                                                                  {
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V1/mv_PARKGAMEII_052_HJ_S1_S3_t1mpr_S256_0.8iso_p2.nii,1'
                                                                  '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/VBM/firstlevel_vbm/V2/mv_PARKGAMEII_052_HJ_19_07_2021_V2_S3_t1mpr_S256_0.8iso_p2.nii,1'
                                                                  }
                                                                  }';
%%
jobs{1}.spm.tools.cat.cat_simple_long.mods = {};
jobs{1}.spm.tools.cat.cat_simple_long.tpm = 'adults';
jobs{1}.spm.tools.cat.cat_simple_long.admin.experimental = 0;
jobs{1}.spm.tools.cat.cat_simple_long.admin.new_release = 0;
jobs{1}.spm.tools.cat.cat_simple_long.admin.lazy = 0;
jobs{1}.spm.tools.cat.cat_simple_long.admin.ignoreErrors = 1;
jobs{1}.spm.tools.cat.cat_simple_long.admin.verb = 2;
jobs{1}.spm.tools.cat.cat_simple_long.admin.print = 2;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.neuromorphometrics = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.lpba40 = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.cobra = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.hammers = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.thalamus = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.ibsr = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.aal3 = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.mori = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.anatomy3 = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.julichbrain = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 1;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
jobs{1}.spm.tools.cat.cat_simple_long.ROImenu.atlases.ownatlas = {''};
jobs{1}.spm.tools.cat.cat_simple_long.fwhm_vol = 6;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.Desikan = 1;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.Destrieux = 1;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.HCP = 0;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.Schaefer2018_100P_17N = 0;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.Schaefer2018_200P_17N = 0;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.Schaefer2018_400P_17N = 0;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.Schaefer2018_600P_17N = 0;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.sROImenu.satlases.ownatlas = {''};
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.fwhm_surf1 = 12;
jobs{1}.spm.tools.cat.cat_simple_long.surface.yes.fwhm_surf2 = 20;

%jobs{1}.spm.tools.cat.cat_simple_long.registration = 0;

jobs{1}.spm.tools.cat.cat_simple_long.ignoreErrors = 1;
jobs{1}.spm.tools.cat.cat_simple_long.nproc = 1;
jobs{1}.spm.tools.cat.cat_simple_long.debug = 1;


spm('defaults','FMRI')
global cat; cat_defaults; cat.extopts.expertgui=0;

[ jobs ] = job_ending_rountines( jobs, skip, par );

end