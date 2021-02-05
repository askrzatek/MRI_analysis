function [job] = spm_job_voi(fspm, roi_struct, fmask, par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function used to create Volume Of Interest files necessary to later VOI analysis      %%
%% INPUTS :                                                                               %%
%%%%%%%%% fspm - list of SPM.mat model paths                                             %%
%%%%%%%%% fmask - list of model mask.image(1) paths : 1 per SPM model                    %%
%%%%%%%%% roi_struct - structure containing lists of voi.name & mask.image               %%
%%%%%%%%%%%%%%%%%%%% roi_struct.name - voi.name string : x per SPM model                 %%
%%%%%%%%%%%%%%%%%%%% roi_struct.froi - mask.image(2) path for each voi : x per SPM model %%
%%%%%%%%% par - jobname, walltime, run                                                   %%
%                                                                                        %%
%% OUTPUT :                                                                              %%
%%%%%%%%% jobs - SPM.mat model files                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                        %%
%% EXEMPLES : %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                                                               %
%fspm            = {'/home/anna.skrzatek/data/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/SPM.mat'}
%fspm            = addsuffixtofilenames(gpath(e),'SPM.mat');                                            %
%fmask           = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/mask.nii,1'}      %
%roi_struct.name = 'Putamen_Right'                                                                                                                                              %
%roi_struct.froi = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_RestingState/Putamen_Right.nii,1'}                                              %
%par.jobname     = 'spm_voi_ts_extract_loop_test'                                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~exist ('par','var')
        par = '';
    end
    
    %% defpar
    defpar.jobname                            = 'spm_voi_ts_extract';
    defpar.walltime                           = '04:00:00';
    defpar.run                                = 0;
    
    par = complet_struct(par, defpar);
    
    %% SPM : util.voi
        %% subject loop
        if iscell(fspm(1))
            nrSubject = length(fspm)
        for i = 1 : nrSubject
        %     job{1}.spm.util.voi.spmmat                = {'/home/anna.skrzatek/data/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/SPM.mat'};
            job{i}.spm.util.voi.spmmat                = fspm(i);
            job{i}.spm.util.voi.adjust                = NaN;
            job{i}.spm.util.voi.session               = 1;
        %     job{1}.spm.util.voi.name                  = 'Putamen_Right';
            job{i}.spm.util.voi.name                  = roi_struct.name;
        %    job{1}.spm.util.voi.roi{1}.mask.image     = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S05_RS/model/model_1/mask.nii,1'};
            job{i}.spm.util.voi.roi{1}.mask.image     = fmask;
            job{i}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        %     job{1}.spm.util.voi.roi{2}.mask.image     = {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_RestingState/Putamen_Right.nii,1'};
            job{i}.spm.util.voi.roi{2}.mask.image     = roi_struct.froi;
            job{i}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            job{i}.spm.util.voi.expression            = 'i1&i2';

        end
            job = do_cmd_sge(job,par);

%end