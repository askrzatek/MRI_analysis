[job] = doublerun_con_batch(fspm, par)

matlabbatch{1}.spm.stats.con.spmmat = {'/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/firstlevel_sts_tapas_doublerun_jan21/SPM.mat'};
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'F-all';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 0 0 0 0 0
                                                        0 1 0 0 0 0
                                                        0 0 1 0 0 0
                                                        0 0 0 1 0 0
                                                        0 0 0 0 1 0
                                                        0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Rest';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'REAL_Left';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 1 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'REAL_Right';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 1 0 0 0];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'IMAGINARY_Left';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 1 0 0];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'IMAGINARY_Right';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 1 0];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'Instruction';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.delete = 0;