clear
clc

par.run     = 0;
par.display = 0;
par.sge     = 1;
par.jobname = 'spm_firstlevel_VOI_ROI';


%% get ROIs & default values

% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/asya.ekmen/AMEDYST/rawIRMf/RSall/ROI_RestingState'; %%ATTENTION : nouvelles analyses avec la session 2 -> nouveau dossier de ROI
% fileROI = cellstr(char(gfile(dirROI,'^r'))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)

% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/studies/AMEDYST/RS/ROI_RS'; % change
% fileROI = cellstr(char(gfile(dirROI,{'^rROI_CbVI_L','^rROI_CbVI_R','^rROI_CbVIII_R','rROI_CingAntL','rROI_CingAntR','^rROI_FrontSupR'}))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)

%% case 1 .mat

dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/secondlevel_sts_tapas_PARK/PARKGAME_all/rois_all_S1_p001_k10';
fileROI = cellstr(char(gfile(dirROI,{'^k10'}))); %fileROI = remove_regex(fileROI,'T1');
char(fileROI)

% ? .nii

dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/secondlevel';
fileROI = cellstr(char(gfile(dirROI,{'^k10'}))); %fileROI = remove_regex(fileROI,'T1');
char(fileROI)


%% case 2 .mat

dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/rois_atlas';
fileROI = cellstr(char(gfile(dirROI,{'roi.mat$'}))); %fileROI = remove_regex(fileROI,'T1');
char(fileROI)

% ? .nii

%% case 3 .nii

dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_SMA_RS'; % change
fileROI = cellstr(char(gfile(dirROI,'^r'))); %fileROI = remove_regex(fileROI,'T1');
char(fileROI)

%%
nRun   = 1;

contrast_PPI = {
    'Effect of Interest'   1
    };


%% prepare job

data_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test';
list_subj = step00_subject_list();
subj_dir  = gdir(data_dir,cellstr2regex(list_subj));
subj_dir  = gdir(data_dir,{'072','074','075','076'}); % change

model_dir = cell(length(subj_dir),1);
for iSubj = 1 : length(subj_dir)
    model_dir{iSubj}(1,1) = gdir(subj_dir{iSubj},'^Presham','RS$','model','model_1$'); % change
    model_dir{iSubj}(2,1) = gdir(subj_dir{iSubj},'^Prestim','RS$','model','model_1$'); % change
    model_dir{iSubj}(3,1) = gdir(subj_dir{iSubj},'^Postsham','RS$','model','model_1$'); % change
    model_dir{iSubj}(4,1) = gdir(subj_dir{iSubj},'^Poststim','RS$','model','model_1$'); % change
end


nRun = 4;

% model_dir = gdir(subj_dir,'^P','RS$','LFF_BOX_glm');
models    = gfile(model_dir,'SPM.mat');


%list_subj  = {'_C01$','_C02$','_C03$','_C04$','_C05$','_C08$','_C09$','_C10$','_C12$','_C13$','_C14$','_C15$','_C16$','_C18$','_C19$','_C20$','_C21$','_C22$','_C23$','_C24$'};
nSubj = length(subj_dir);
nROI  = length(fileROI);
nCon  = size(contrast_PPI,1);

%jobs = cell(nSubj,nROI,nCon);
jobs = cell(nSubj,nROI,nRun);

% Fetch fmri files
%volume_dir = get_subdir_regex_multi(subj_dir,'prep[123]$');
% volume_dir = get_subdir_regex(subj_dir,'^P','RS$','^tedana');
% volume_file = cell(size(volume_dir));
% for iSubj = 1 : nSubj
%     volume_file{iSubj} = char(get_subdir_regex_files(volume_dir{iSubj},'s8wdn.*nii$'));
% end % iSubj

volume_dir  = cell(length(subj_dir),1);
volume_file = cell(length(subj_dir),1);
for iSubj = 1 : length(subj_dir)
    volume_dir{iSubj}(1,1) = gdir(subj_dir{iSubj},'^Presham','RS$','^tedana'); % change
    volume_dir{iSubj}(2,1) = gdir(subj_dir{iSubj},'^Prestim','RS$','^tedana'); % change
    volume_dir{iSubj}(3,1) = gdir(subj_dir{iSubj},'^Postsham','RS$','^tedana'); % change
    volume_dir{iSubj}(4,1) = gdir(subj_dir{iSubj},'^Poststim','RS$','^tedana'); % change
   
    for i = 1 : 4
        volume_file{iSubj}(i,1) = gfile(volume_dir{iSubj}(i,1),'^s6wts');
    end
   
end

% Verif VOI


%% Prepare jobs

% Onsets identical for everybody (resting state)
TR=1.6;
nbvol = 300;

for iSubj = 1 : nSubj
   
   
    for iROI = 1 : nROI
       
        [~,roi_name] = fileparts(fileROI{iROI});
       
        for iRun = 1 : nRun
           
            stat_dir = get_parent_path( model_dir{iSubj}(iRun) );
           
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.dir = fullfile( stat_dir, sprintf('Modele_VOI__%s',roi_name) );%VOI__rCerebellum_6_Left__run1_1
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.units = 'secs';
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.RT = 1.6;
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.fmri_t = 16;
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
           
            % VOI
            voi_file = fullfile(model_dir{iSubj}{iRun},sprintf('VOI__%s__run1_1.mat',roi_name));
            voi = load(voi_file);
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.sess.regress(1).name = 'Y';
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.sess.regress(1).val = voi.Y;
           
            % Volumes
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.sess.scans = spm_select('expand',cellstr(volume_file{iSubj}{iRun}));

            % Other
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.sess.multi = {''}; % we don't use .mat file for the onsets
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.sess.hpf = 128;
           
        end % iRun
       
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.volt = 1;
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.global = 'None';
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.mthresh = 0; % because already skullstriped before TEDANA
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.mask = {''};
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.cvi = 'AR(1)';
       
       
       
    end % iROI
   
end % iSubj

%% verification creation VOIs

for iSubj = 1 : nSubj
   
   
    for iROI = 1 : nROI
       
        [~,roi_name] = fileparts(fileROI{iROI});
       
        for iRun = 1 : nRun
            voi_file = gfile(model_dir{iSubj}{iRun},sprintf('VOI__%s__run1_1.mat',roi_name));
        end
    end
end

%% Run
% return
jobs = jobs(:);
% clear par
% par.display=1;
%jobs = jobs(1:10);
job_ending_rountines(jobs,[],par)

%%  Estimate


% model_dir = gdir(subj_dir,'^P','RS$','LFF_BOX_glm');


Modele_roi_dir = gdir(subj_dir,'^P','RS$','model','^Modele_VOI'); % change ?
SPM_roi_file = gfile(Modele_roi_dir,'SPM.mat');

clear par
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname  = 'spm_firstlevel_ROI_est';
job_first_level_estimate(SPM_roi_file,par)

%%  Define contrasts

clear par
par.sge=1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'spm_firstlevel_ROI_con';


PositiveEffectROI= [1 0];


contrast.names={'PositiveEffectROI'};
contrast.values={PositiveEffectROI};
contrast.types={'T'};

par.delete_previous=1

j=job_first_level_contrast(SPM_roi_file,contrast,par)