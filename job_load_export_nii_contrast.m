%-----------------------------------------------------------------------
% Job saved on 22-Apr-2020 00:02:51 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
function jobs = test_batch()

par.subdir = 'ANOVA2x2_LxT';
fsub = {'/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti/secondlevel_ACTIVATION_PARK_S1/'};
fspm = char(get_subdir_regex(fsub(1), par.subdir));
SPM_search = char(addsuffixtofilenames(fspm, '/SPM.mat'));
if exist(SPM_search, 'file')
    SPM_found = get_subdir_regex_files (fspm, 'SPM.mat');
end

load (char(SPM_found));

spm('defaults', 'FMRI');

    jobs{1}.spm.stats.results.spmmat = {'/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti/secondlevel_ACTIVATION_PARK_S1/ANOVA2x2_LxT/SPM.mat'};
    jobs{1}.spm.stats.results.conspec.titlestr = '';
    jobs{1}.spm.stats.results.conspec.contrasts = [1 2];
    jobs{1}.spm.stats.results.conspec.threshdesc = 'none';
    jobs{1}.spm.stats.results.conspec.thresh = 0.01;
    jobs{1}.spm.stats.results.conspec.extent = 5;
    jobs{1}.spm.stats.results.conspec.conjunction = 1;
    jobs{1}.spm.stats.results.conspec.mask.none = 1;
    jobs{1}.spm.stats.results.units = 1;
    jobs{1}.spm.stats.results.export{1}.nary.basename = 'clusters_ROI1';

    spm_jobman('run', jobs);
end