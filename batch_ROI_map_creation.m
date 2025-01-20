function matlabbatch = batch_ROI_map_creation(imgs, exp, output, par)
%%
matlabbatch{1}.spm.util.imcalc.input = {
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg1.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg10.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg12.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg13.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg14.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg15.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg155.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg156.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg16.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg19.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg2.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg3.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg33.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg34.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg37.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg38.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg4.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg5.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg6.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg61.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg62.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg65.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg67.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg68.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg7.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg71.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg72.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg74.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg78.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg8.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg85.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg86.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg87.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg88.nii,1'
                                        '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01/CAREN_RSN01_AAL_Reg9.nii,1'
                                        };
%%
matlabbatch{1}.spm.util.imcalc.output = 'CAREN_Sallience_Network';
matlabbatch{1}.spm.util.imcalc.outdir = {'/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN'};
matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26+i27+i28+i29+i30+i31+i32+i33+i34+i35';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = -1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;