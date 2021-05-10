function [jobs] = spm_job_voi(fspm, fmask, par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function used to create Volume Of Interest files necessary to later VOI analysis      %%
%% INPUTS :                                                                               %%
%%%%%%%%% fspm - list of SPM.mat model paths                                             %%
%%%%%%%%% fmask - list of model mask.image(1) paths : 1 per SPM model                                 %%
%%%%%%%%% par - roi_dir, jobname, walltime, run                                                       %%
%
%% INTER_VARS:
%%%%%%%%% roi_struct - structure containing lists of voi.name & mask.image               %unnecessary %%
%%%%%%%%%%%%%%%%%%%% roi_name - voi.name string : x per SPM model                                     %%
%%%%%%%%%%%%%%%%%%%% roi_froi - mask.image(2) path for each voi : x per SPM model        %unnecessary %%
%
%% OUTPUT :                                                                                           %%
%%%%%%%%% jobs - SPM.mat model files                                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                     %%
%% EXAMPLES : %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                                                               %
%fspm            = {'/home/anna.skrzatek/data/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/SPM.mat'}
%fspm            = addsuffixtofilenames(gpath(e(i)),'SPM.mat');                                            %
%
%fmask           = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/mask.nii,1'}      %
%fmask           = {sprintf('%s,1',string(addsuffixtofilenames(gpath(e(i)),'mask.nii')))};
%
%fmask           = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/wbet_Tmean_vtde1_mask.nii,1'}      %
%fmask           = addsuffixtofilenames(gpath(e.gser('run_RS'),'wbet_Tmean_vtde1_mask.nii')))};
%fmask           = e.gser('run_RS').gvol('wmask') .toJob ;
%
%roi_name = 'Putamen_Right'                                                                                                                                                      %
%roi_name = {'Putamen_R','Putamen_L'};                                                                                                                                           %
%
%% NOT NEEDED anymore for the roi_name and par.roi_dir are enough to find the roi files
%roi_froi = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_RestingState/Putamen_Right.nii,1'}
%roi_struct_path{iroi} = sprintf('%s/%s.nii,1',par.roi_dir,roi_name{iroi}) ;
%roi_froi = roi_struct_path{:}
%
%par.jobname     = 'spm_voi_ts_extract_loop_test'                                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~exist ('par','var')
        par = '';
    end
    
    %% defpar
    defpar.jobname                            = 'spm_voi_ts_extract';
    defpar.walltime                           = '04:00:00';
    defpar.run                                = 0;
    defpar.roi_dir                            = '/home/anna.skrzatek/data/nifti_test/firstlevel_sts_tapas_doublerun_jan21/wbet_mask/ROI_aal_pariet_mot_premot_cereb_BG';
    defpar.sge                                = 1;
    
    par = complet_struct(par, defpar);
    
    %% SPM : util.voi
    idx = 1;
    %% subject loop
    if iscell(fspm(1)) %% not sure what to do if iscell = 0
        nrSubject = length(fspm);
    end
    for i = 1 : nrSubject
        roi_froi = cellstr(char(get_subdir_regex_files(par.roi_dir,'.*'))).';
        for j = 1 : length(roi_froi)
            [~,roi_name] = fileparts(roi_froi{j}); 
        %     jobs{1}.spm.util.voi.spmmat                = {'/home/anna.skrzatek/data/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/SPM.mat'};
            jobs{idx}.spm.util.voi.spmmat                = fspm(i);
            jobs{idx}.spm.util.voi.adjust                = NaN;
            jobs{idx}.spm.util.voi.session               = 1;
        %     jobs{1}.spm.util.voi.name                  = 'Putamen_Right';
            jobs{idx}.spm.util.voi.name                  = roi_name;
        %    jobs{1}.spm.util.voi.roi{1}.mask.image     = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/mask.nii,1'};
            jobs{idx}.spm.util.voi.roi{1}.mask.image     = spm_select('expand',fmask(i));
            jobs{idx}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        %     jobs{1}.spm.util.voi.roi{2}.mask.image     = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_RestingState/Putamen_Right.nii,1'};
%             jobs{idx}.spm.util.voi.roi{2}.mask.image     = spm_select('expand',cellstr(sprintf('%s/%s.nii',par.roi_dir,roi_name{j})));
            jobs{idx}.spm.util.voi.roi{2}.mask.image     = spm_select('expand',cellstr(roi_froi{j}));
            
            jobs{idx}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            jobs{idx}.spm.util.voi.expression            = 'i1&i2';
            idx = idx +1;
        end
    end
    skip = [];

    spm('defaults','FMRI')
    [ jobs ] = job_ending_rountines( jobs, skip, par );
end