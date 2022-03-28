clear
clc

%addpath '/network/lustre/iss01/cenir/analyse/irm/users/cecile.gallea/ASYA/asyasuit/'

par.run     = 0;
par.display = 0;
par.sge     = 1;
par.jobname = 'spm_firstlevel_VOI_resliced_double_wbet';


%% get ROIs & default values

% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/asya.ekmen/AMEDYST/rawIRMf/RSall/ROI_RestingState'; %%ATTENTION : nouvelles analyses avec la session 2 -> nouveau dossier de ROI
% fileROI = cellstr(char(gfile(dirROI,'^r'))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)

% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/studies/AMEDYST/RS/ROI_RS'; % change
% fileROI = cellstr(char(gfile(dirROI,{'^rROI_CbVI_L','^rROI_CbVI_R','^rROI_CbVIII_R','rROI_CingAntL','rROI_CingAntR','^rROI_FrontSupR'}))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)

% %% case 1 .mat
% 
% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/secondlevel_sts_tapas_PARK/PARKGAME_all/rois_all_S1_p001_k10';
% fileROI = cellstr(char(gfile(dirROI,{'^k10'}))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)
% 
% % ? .nii
% 
% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/secondlevel';
% fileROI = cellstr(char(gfile(dirROI,{'^k10'}))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)
% 
% 
% %% case 2 .mat
% 
% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/rois_atlas';
% fileROI = cellstr(char(gfile(dirROI,{'roi.mat$'}))); %fileROI = remove_regex(fileROI,'T1');
% char(fileROI)
% 
% % ? .nii

%% case 3 .nii

% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_aal_pariet_mot_premot_cereb_BG';
% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_SMA_RS';
% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_RestingState';
% dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_Cereb_RS';
dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_pariet_mot_premot_cereb_BG_PPN';

fileROI = cellstr(char(gfile(dirROI,'.*'))); %fileROI = remove_regex(fileROI,'T1');

%pccmodeldir = get_subdir_regex(subj_dir,'.*RS','model','model_2');
%fileROI = cellstr(char(gfile(pccmodeldir,'.*PCC.*.mat')));

char(fileROI)

%% reslice fileROI to have the same voxel size (2.5) as the func img

%%
nRun   = 1;

contrast_PPI = {
    'Effect of Interest'   1
    };


%% prepare job
% % verifs
% main_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test';
% % subj_dir = gdir(main_dir,'^Subj|^___S')
% subj_dir = gdir(main_dir,'PARKGAME.*[a,c]$')
% SessDir = gdir(subj_dir,'.*RS$');
% StatDir = gdir(SessDir,'model','model_1');

data_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test';

%list_subj = step00_subject_list();
list_subj ={'PARKGAMEII_001_NB.*a$','PARKGAMEII_002_BM.*a$','PARKGAMEII_00[1,3]_SM.*c$','PARKGAMEII_007_SD.*a$','PARKGAMEII_008_JR.*a$','PARKGAMEII_023_LJ.*c$','PARKGAMEII_025_CA.*a$','PARKGAMEII_028_PC.*c$','PARKGAMEII_033_DD.*c$','PARKGAMEII_039_KM.*a$','PARKGAMEII_040_RE.*a$','PARKGAMEII_042_RS.*a$','PARKGAMEII_043_PD.*a$','PARKGAMEII_044_CK.*c$','PARKGAMEII_046_HJ.*c$','PARKGAMEII_047_BF.*c$','PARKGAMEII_048_SB.*a$','PARKGAMEII_052_HJ.*c$','PARKGAMEII_053_LM.*a$'}; 
%list_subj ={'PARKGAMEII_052_HJ.*c$','PARKGAMEII_053_LM.*a$'}
subj_dir  = gdir(data_dir,cellstr2regex(list_subj));
%subj_dir  = gdir(data_dir,{'072','074','075','076'}); % change

model_dir = cell(length(subj_dir),1);
for iSubj = 1 : length(subj_dir)
    model_dir{iSubj}(1,1) = gdir(subj_dir{iSubj},'RS$','model','^model_2$'); % change
end


%nRun = 2; % when analysis done in double-run condition
nRun  = 1;

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
    volume_dir{iSubj}(1,1) = gdir(subj_dir{iSubj},'RS$','tedana009a1_vtd'); % change   
    volume_file{iSubj} = gfile(volume_dir{iSubj},'^s6wts');   
end

% Verif VOI
%% verification creation VOIs
for iSubj = 1 : nSubj
% fileROI_id = cell(nROI,1);   
    for iROI = 1 : nROI
       
        [~,roi_name] = fileparts(fileROI{iROI});
       
        for iRun = 1 : nRun
            voi_file = gfile(model_dir{iSubj}{iRun},sprintf('VOI_%s_1.mat',roi_name));
%             voi_files = gfile(model_dir{iSubj}{iRun},sprintf('VOI_%s_1.mat','.*'));
%             for iROI = 1 : nROI
%                 [~,roi_name] = fileparts(voi_files{1}(iROI,:));
%                 fileROI_id{iROI} = roi_name(5:length(roi_name)-2);
%                 clear roi_name
%             end
%             fileROI_id
%             clear fileROI_id
        end
    end
end


%% Prepare jobs

% Onsets identical for everybody (resting state)
TR    = 1.6;
nbvol = 300;
par.redo = 1;

for iSubj = 1 : nSubj
    
    skip = [];
        
    for iROI = 1 : nROI
       
        [~,roi_name] = fileparts(fileROI{iROI}); % not aplied since we create one VOI individually from a sphere at coordinates (PCC)

%         roi_name = fileROI_id{iROI};
        for iRun = 1 : nRun
           
            stat_dir = get_parent_path( model_dir{iSubj}(iRun) );
            
%             beta_file = fullfile(stat_dir,'beta_0001.nii');
%             if ~par.redo   &&  exist(beta_file,'file')
%                 skip = [skip iSubj];
%                 fprintf('[%s]: skiping subj %d because %s exist \n',mfilename,iSubj,beta_file);
%             end
            
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.dir = fullfile( stat_dir, sprintf('Modele_VOI__%s', roi_name) ); %'PCC')); %VOI__rCerebellum_6_Left__run1_1
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.units = 'secs';
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.RT = 1.6;
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.fmri_t = 16;
            jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
           
            % VOI
            voi_file = fullfile(model_dir{iSubj}{iRun},sprintf('VOI_%s_1.mat', roi_name)); %'PCC')); %
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
%         jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.mthresh = 0; % because already skullstriped before TEDANA
%         jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.mask = {''};
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.mthresh = 0.1; % consistent with previous mask parameters -> according to Benoit
        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.mask = fullfile(model_dir{iSubj}{iRun},'mask.nii'); % masks from the firstlevel of model_1 or model_2 

        jobs{iSubj,iROI,iRun}.spm.stats.fmri_spec.cvi = 'AR(1)';
       
       
       
    end % iROI
   
end % iSubj


%% Run
% return
jobs = jobs(:);
% clear par
% par.display=1;
%jobs = jobs(1:10);
par.sge_queu = 'normal,bigmem';
par.walltime = '01:00:00';
par.run = 0;
par.sge = 1;

job_ending_rountines(jobs,[],par);

%%  Estimate


% model_dir = gdir(subj_dir,'^P','RS$','LFF_BOX_glm');
clear par
clear Modele_roi_dir
clear SPM_roi_file

par.sge = 1;
par.run = 0;

par.sge_queu = 'normal,bigmem';
par.walltime = '01:00:00';

for ir = 1 : nROI

    for subj = 1 : nSubj

        [~,roi_name] = fileparts(fileROI{ir});
        roi_name = roi_name(1:end-10);

        Modele_roi_dir{subj} = gdir(subj_dir{subj},'RS$','model',sprintf('Modele_VOI__%s.*',roi_name));
        SPM_roi_file{subj} = char(gfile(Modele_roi_dir{subj},'SPM.mat'));

    end
    
    cd (main_dir)
    
    par.jobname  = sprintf('spm_firstlevel_VOI_est_%s',roi_name);
    job_first_level_estimate(SPM_roi_file,par);

end

%Modele_roi_dir = gdir(subj_dir,'RS$','model','^Modele_VOI_'); % change ?

%SPM_roi_file = gfile(Modele_roi_dir,'SPM.mat');

%clear par
%par.sge = 1;
%par.sge_queu = 'normal,bigmem';
%par.walltime = '01:00:00';

%par.jobname  = 'spm_firstlevel_VOI_est_double_wbet';
%job_first_level_estimate(SPM_roi_file,par);

%%  Define contrasts

clear par
par.sge=0;
par.run = 1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'spm_firstlevel_ROI_con_double_wbet_pcc';


PositiveEffectROI= [1 0];


contrast.names={'PositiveEffectROI'};
contrast.values={PositiveEffectROI};
contrast.types={'T'};

par.delete_previous = 1;
for ir = 1 : nROI

    for subj = 1 : nSubj

        [~,roi_name] = fileparts(fileROI{ir});
        roi_name = roi_name(1:end-10);

        Modele_roi_dir{subj} = gdir(subj_dir{subj},'RS$','model',sprintf('Modele_VOI__%s.*',roi_name));
        SPM_roi_file{subj} = char(gfile(Modele_roi_dir{subj},'SPM.mat'));

    end

    j=job_first_level_contrast(SPM_roi_file,contrast,par);
    
end

%% Create the directories in subject firstlevel_RS/
firstlevel_dir = fullfile(main_dir,'firstlevel_RS');
list_subj ={'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_040_RE_a','PARKGAMEII_042_RS_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_046_HJ_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_053_LM_a'};

for ir = 1 : nROI
    for subj = 1 : length(list_subj)
        outdir = char(gdir(firstlevel_dir,list_subj{subj}));
        
        [~,roi_name] = fileparts(fileROI{ir});
        roi_name = roi_name(1:end-10);
        
        mkdir(outdir,sprintf('%s_V1',roi_name));
        mkdir(outdir,sprintf('%s_V2',roi_name));
    end
end
