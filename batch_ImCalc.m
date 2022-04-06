function matlabbatch = batch_ImCalc(imgs, exp, output, par)
%%

matlabbatch{1}.spm.util.imcalc.input = {
                                        '/home/anna.skrzatek/data/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2018_07_25_PARKGAMEII_002_BM_25_07_2018_V1_a/S14_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2018_08_29_PARKGAMEII_001_NB_29_08_2018_V2_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2018_09_05_PARKGAMEII_002_BM_05_09_2018_V2_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2018_09_19_PARKGAMEII_001_SM_19_09_2018_V1_c/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2018_10_26_PARKGAMEII_003_SM_26_10_2018_V2_c/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_01_09_PARKGAMEII_007_SD_09_01_2019_V1_a/S12_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_01_18_PARKGAMEII_008_JR_18_01_2019_V1_a/S15_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_02_15_PARKGAMEII_007_SD_15_02_2019_V2_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_03_06_PARKGAMEII_008_JR_06_03_2019_V2_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_05_24_PARKGAMEII_023_LJ_24_05_2019_V1_c/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_05_29_PARKGAMEII_025_CA_29_05_2019_V1_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_06_21_PARKGAMEII_028_PC_21_06_2019_V1_c/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_07_10_PARKGAMEII_025_CA_10_07_2019_V2_a/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_07_12_PARKGAMEII_023_LJ_12_07_2019_V2_c/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_07_31_PARKGAMEII_028_PC_31_07_2019_V2_c/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_09_11_PARKGAMEII_033_DD_11_09_2019_V1_c/S12_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2019_11_06_PARKGAMEII_033_DD_06_11_2019_V2_c/S12_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_01_08_PARKGAMEII_039_KM_08_01_2020_V1_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_02_21_PARKGAMEII_043_PD_21_02_2020_V1_a/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_02_26_PARKGAMEII_039_KM_26_02_2020_V2_a/S11_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/S12_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_05_25_PARKGAMEII_043_PD_25_05_2020_V2_a/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_07_27_PARKGAMEII_044_CK_27_07_2020_V2_c/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2020_11_04_PARKGAMEII_048_SB_04_11_2020_V1_a/S10_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        '/home/anna.skrzatek/data/nifti_test/2021_01_06_PARKGAMEII_048_SB_06_01_2021_V2_a/S12_ACTIVATION/wbet_Tmean_vtde1_mask.nii,1'
                                        };
%%
matlabbatch{1}.spm.util.imcalc.output = '/home/anna.skrzatek/data/nifti_test/wmean_mask.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26)/26';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;