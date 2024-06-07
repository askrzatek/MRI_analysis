%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script adapté par A.SKRZATEK le 09/04/2024/ de Salim OUARAB, 
%%% d'après l'article et le conseil d'Ottavia Dipasquale (NeuroImage, 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% To be executed after the Preproc & Analyse_RSFC_toolbox : as we will need inputs:
%%% e.mat
%%% bp_clean.nii 
clc
clear

addpath('/network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/dev_matvol/')
addpath('/network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/MNI/pet_atlas/')
addpath('/network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/toolbox/')

main_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';

cd (main_dir)
load('e_nonchir.mat');

suj   = gpath(e.getSerie('run_RS'));

dfunc = gdir(suj,'tedana009a1_vt');
dclean = gdir(dfunc,'rsfc');
rp_file = gfile(dfunc,'multiple_reg.*txt');

fclean = gfile(dclean,'^bp_clean.nii');

%% Alternative version from Cecile : yet completely wrong
% dfunc = gdir(suj,'ALFF');
% dclean = dfunc;
% 
% fclean = gfile(dclean,'^ess.*');

% for is = 1: length(dfunc)
%     mkdir(dfunc{is},'React_Cecile_model')
% end

dreact = get_subdir_regex(dfunc(:),'React$'); %.*Cecile');

r_movefile(fclean, dreact,'copy')
r_movefile(rp_file,dreact,'copy')

frclean = gfile(dreact,'^bp_clean.nii');
% frclean = gfile(dreact,'^ess');
% rprclean = gfile(dreact,'multiple');

%% Noradrenaline map

datlas = {'/network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/MNI/pet_atlas'}
ref_NA = gfile(datlas,'Atlas_NAT_IntensityNorm_1.nii');

% %% original directory with masks
% datlas = {'/network/lustre/iss02/cenir/analyse/irm/users/sandy.mournet/ToM_2021/IRM/ToM_BOLD/Scripts_Roxane/Réseaux/'};

% %% Acetylcholine map
% dace   = gdir(datlas,'_Acétylcholine')
% ref_ach = gfile(dace,'Atlas_VACh_IntensityNorm.nii');

% %% Dopamine map
% ddop   = gdir(datlas,'_Dopamine')
% ref_da = gfile(ddop,'Atlas_DAT_IntensityNorm.nii')

%% Reslicing grey matter mask to NT maps

% cmd = 'python3
% /network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/dev_matvol/python/react_masks.py
% /network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/data/subject_list.txt
% /network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/MNI/pet_atlas/pet_atlas_4d.nii
% /network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/MNI/pet_atlas/rgm_mask.nii.gz
% out_masks' % unused command ??

mask = gfile({'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/pet_atlas'},'^gm.*')
% mask = unzip_volume(mask)
%%%%%%%%%
clear par
par.sge = 0;
par.run = 1;
par.prefix = 'rNA_';
par.jobname = 'rNA'

job_reslice(mask,ref_NA,par)

%% NORADRENALINE reference map bp_clean reslice

clear par
par.sge = 1;
par.run = 0;
par.prefix = 'rNA_';
par.jobname = 'rNA';

job_reslice(frclean,repmat(ref_NA,50,1),par)

% %% DOPAMINE reference map grey matter mask reslice
% 
% clear par
% par.sge = 0;
% par.run = 1;
% par.prefix = 'rda_';
% par.jobname = 'rDA'
% 
% job_reslice(frclean,refd,par)

frna  = gfile(dreact,'^rNA_.*nii');
% frna  = gfile(dreact,'^rNA_bp_clean.*nii');

% frda   = gfile(dreact,'^rda_bp_clean.*nii')

%% attention : it's important to use module load fsl before executing

% JOB PREPARE
for i=1:length(suj)
    cmd1 = sprintf('fsl_glm -i %s -d %s -o %sts_NA_clean_4D.txt'...
        ,frna{i},ref_NA{1},dreact{i}) %NA
    jobs{i,1}=cmd1
end

% JOB EXECUTE
clear par
par.sge = 1
par.jobname = 'MaskReact_NA_job'
par.walltime      = '05';
par.mem      = '100G';
do_cmd_sge(jobs,par);

% %% DOPAMINE
% 
% clear jobs
% for i=1:length(suj)
%     cmd1 = sprintf('fsl_glm -i %s -d %s -o %sts_da_clean_4D.txt'...
%         ,frda{i},ref_da{1},dreact{i}) 
%     jobs{i,1}=cmd1
% end
% 
% clear par
% par.sge = 1
% par.jobname = 'SpaRegDA'
% par.walltime      = '05';
% par.mem      = '100G';
% do_cmd_sge(jobs,par);

%% Step2 : second individual glm

ftxt_na   = gfile(dreact,'ts_NA.*txt')
dirStats   = r_mkdir(dreact,'GLM_NA')

clear par
par.file_reg = '^bp_clean.*nii';

% par.file_reg = '^ess';
% par.rp = 1;
% par.rp_file = rp_file;
%par.rp_regex = 'multiple.*txt';
% par.mask_thr = -inf; 
par.user_regressor = ftxt_na; 
par.TR  = 1.6;

par.redo = 0;
par.run = 1;
par.sge = 0;
par.jobname = 'GLM_NA'
par.walltime = '05';
par.mem      = '16G';
par.output = dirStats; 

addpath('/home/anna.skrzatek/MRI_analysis/')
addpath('/home/anna.skrzatek/MRI_analysis/AUDICOG')
jobs = job_RS_model_specify(dreact,dirStats,par)
% jobs = job_first_level_specify(dreact,dirStats,[],par)

% job_rsfmri_model_specify(dreact,par)

%% ESTIMATION

fspm = gfile(dirStats,'SPM.mat');

clear par
par.write_residuals = 0;

par.run = 0;
par.sge = 1;
par.jobname = 'EstGlmNA'
par.mem      = '32G';
job_first_level_estimate(fspm,par)

% %% DA
% 
% ftxt_da    = gfile(dreact,'da.*txt')
% dirStats   = r_mkdir(dreact,'GLM_DA')
% 
% 
% 
% clear par
% par.file_reg = '^bp_clean.*nii';
% par.mask_thr = -inf; 
% par.file_regressor = ftxt_da; 
% par.TR  = 2.044;
% 
% par.run = 0;
% par.sge = 1;
% par.jobname = 'GLM_DA1'
% par.walltime = '05';
% par.mem      = '16G';
% par.output = dirStats; 
% 
% %jobs = job_first_level_specify(dreact,dirStats,onsets,par)
% 
% job_rsfmri_model_specify(dreact,par)
% fspm = gfile(dirStats,'SPM.mat')
% 
% clear par
% par.write_residuals = 0;  % pas besoin  ---> 0
% 
% par.run = 0;
% par.sge = 1;
% par.jobname = 'EstGlmDa1'
% par.mem      = '32G';
% 
% 
% job_first_level_estimate(fspm,par)
% 


%% Stats :

%% A.SKRZATEK
project_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';

d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;
tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv'));
covars.name = {'Age','Sex','Hearing_loss','Emotion','Alert_score_ANT'};

groups.name = {'Control' 'Tinnitus'};
pathway_contrasts = {'/GLM_NA/'};
% contrast_names = {'seed2voxel_pearson_ParaHipp_BA_31_Orb_PFC_lAudio_rAudio__Cingulate'};
contrast_names = {'beta.*01.nii$'};

%% getting scans & covariates organised per group
scans = cell(1,length(groups.name));
for igroup = 1:length(groups.name)
    j = 0 ;
    for iSubj = 1:length(e)
        
        ifile = e(iSubj).name;
        id = str2double(ifile(25:end)) ;
        subj_group = d.Groupe( find (d.Num_IRM == id) )  ;
        
        if subj_group == igroup
            j = j + 1 ;
            for icov = 1:length(covars.name)
                scans{igroup}.cov{j,1} = tab.Age(tab.code_IRM == id);
                scans{igroup}.cov{j,2} = tab.Genre(tab.code_IRM == id);
                scans{igroup}.cov{j,3} = tab.pca_audio1(tab.code_IRM == id);
                scans{igroup}.cov{j,4} = str2double(tab.pca_emotionnel1(tab.code_IRM == id));
                scans{igroup}.cov{j,5} = tab.log_ANT_RT_Alerting(tab.code_IRM == id);
                scans{igroup}.cov{j,6} = tab.ANT_STD_mean(tab.code_IRM == id);
%                 scans{igroup}.cov{j,5} = tab.ANT_STD_mean(tab.code_IRM == id);
            end
            scans{igroup}.contrast{j} = char(get_subdir_regex_files(fullfile(dreact{iSubj}, pathway_contrasts{1}), contrast_names{1})) ; % if multiple pathway_contrasts then change 1 to icontr
        end
    end
end

%% NA contrast

groups.val = {scans{1}.contrast(1,:), scans{2}.contrast(1,:)}; % NA

%% 2-sample t-test model specification

jobnames = {'2sample_ttest_RS_pca','2sample_ttest_RS_Alert_RT_STD_inter_pca', '2sample_ttest_RS_Alert_pca_model', '2sample_ttest_RS_RT_STD_pca_model','2sample_ttest_RS_Alert_wo_PCA_RT_STD_covar'};
outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/REACT/2sample_ttest_NA_pca', '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/REACT/2sample_ttest_NA_mixed_interaction_pca', '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/REACT/2sample_ttest_NA_Alert_pca', '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/REACT/2sample_ttest_NA_RT_STD_pca', '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/REACT/2sample_ttest_NA_Alert_covar_RT_STD_wo_pca'};

clear par
par.sge = 0;
par.run = 1;

for ijob = 1: length(jobnames)
    par.jobname = jobnames{ijob};
    
    % covars according to the model tested : Alert or RT_STD (or mixed in
    % the future)
    if ijob == 1
        par.covars = 1;
        par.intercov = [1,1,1,1,1,1];
        covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting', 'STD_RT_ANT'};
        covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5});vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6})};
        
%         par.intercov = [1,1,1,1];
%         covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1'};
%         covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4})};    
    elseif ijob == 2
        par.covars = 1;
        par.intercov = [1,1,1,1,2,2];
        covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting', 'STD_RT_ANT'};
        covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5});vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6})};
    elseif ijob == 3
        par.covars = 1;
        par.intercov = [1,1,1,1,2];
        covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting'};
        covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5})};
    elseif ijob == 4
        par.covars = 1;
        par.intercov = [1,1,1,1,2];
        covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','STD_RT_ANT'};
        covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6})};
    else
        par.covars = 1;
        par.intercov = [1,1,1,2];
        covars.name = {'Age','Genre','STD_RT_ANT','log_ANT_RT_Alerting'};
        covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5})};
    end
    
    varcov_2nd_level_2sample_model_spec(groups.val,outdirs(ijob),covars,par);
end

%% Estimation

for ir = 1 : length(jobnames)
    fspm{ir} = addsuffixtofilenames(outdirs(ir),'/SPM.mat');
    
    %% Model estimate
    clear par
    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = sprintf('spm_est_%s',jobnames{ir});
    job_first_level_estimate(fspm{ir},par)

end

%% Contrast definition

%% Single variable of interest
% F-stat
Main_NA_effect = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
Main_group_effect_Alert = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
Differential_group_Alert_effect = [0 0 0 0 0 0 1 -1];
Main_group_effect_RT_STD = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
Differential_group_RT_STD_effect = [0 0 0 0 0 0 1 -1];
Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0];
Main_effect_Emotion = [0 0 0 0 0 1 0 0];
Main_effect_Age = [0 0 1 0 0 0 0 0];
Main_effect_Sex = [0 0 0 1 0 0 0 0];

% t-stat
Control = [1 0];
Tinnitus =[0 1];

%% Contrast : write
clear par

par.sge = 0;
par.run = 1;
par.display = 0;
par.jobname = 'spm_write_con';

par.sessrep = 'none';

par.delete_previous = 1;
par.report          = 0;

for ir = 1 : length(jobnames)
    if ir == 1
        contrast_F.names = {
            'Main_Group_Effect_on_NA'
            'Main_effect_of_Alert_Score'
            'Main_effect_of_RT_STD'
            'Main_effect_Hearing_Loss'
            'Main_effect_Emotion'
            'Main_effect_Age'
            'Main_effect_Sex'
        }';
        contrast_F.values = {
            [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 1]
            [0 0 0 0 1 0 0 0]
            [0 0 0 0 0 1 0 0]
            [0 0 1 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]
        }';
    elseif ir == 2
        contrast_F.names = {
            'Main_Group_NA_Effect'
            'Main_Group_NA_Alert_Score_Interaction'
            'Differential_Group_NA_effect_on_Alert_Score'
            'Main_Group_NA_RT_STD_Interaction'
            'Differential_Group_NA_effect_on_RT_STD'
            'Main_NA_effect_Hearing_Loss'
            'Main_NA_effect_Emotion'
            'Main_NA_effect_Age'
            'Main_NA_effect_Sex'
            }';

        contrast_F.values = {
            [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 1 -1 0 0]
            [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 1 -1]
            [0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 0 0]
            [0 0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 0]
        }';
    elseif ir == 3
        contrast_F.names = {
            'Main_NA_Effect'
            'Main_Group_effect_on_Alert_Score'
            'Differential_Group_effect_on_Alert_Score'
            'Main_effect_Hearing_Loss'
            'Main_effect_Emotion'
            'Main_effect_Age'
            'Main_effect_Sex'
            }';

        contrast_F.values = {
            Main_NA_effect
            Main_group_effect_Alert
            Differential_group_Alert_effect
            Main_effect_Hearing_Loss
            Main_effect_Emotion
            Main_effect_Age
            Main_effect_Sex
        }';
    elseif ir == 4
        contrast_F.names = {
            'Main_NA_Effect'
            'Main_Group_effect_on_RT_STD'
            'Differential_Group_effect_on_RT_STD'
            'Main_effect_Hearing_Loss'
            'Main_effect_Emotion'
            'Main_effect_Age'
            'Main_effect_Sex'
        }';
        contrast_F.values = {
            Main_NA_effect
            Main_group_effect_RT_STD
            Differential_group_RT_STD_effect
            Main_effect_Hearing_Loss
            Main_effect_Emotion
            Main_effect_Age
            Main_effect_Sex
        }';
    else
        contrast_F.names = {
            'Main_Group_Effect_on_NA'
            'Main_effect_of_Alert_Score_on_NA'
            'Main_differential_Alert_Score_effect_on_NA'
            'Main_effect_of_RT_STD_on_NA'
            'Main_effect_Age'
            'Main_effect_Sex'
        }';
        contrast_F.values = {
            [1 0 0 0 0 0 0; 0 1 0 0 0 0 0]
            [0 0 0 0 0 1 0; 0 0 0 0 0 0 1]
            [0 0 0 0 0 1 -1]
            [0 0 0 0 1 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 1 0 0 0]
        }';
    end
    
%% t-contrast definition

    contrast_T.names = {
    'Control-Tinnitus'
    'Tinnitus-Control'
    }';

    contrast_T.values = {
        Control-Tinnitus
        Tinnitus-Control
    }';

    contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));
    contrast_T.types = cat(1,repmat({'T'},[1,length(contrast_T.names)]));

    contrast.names  = [contrast_T.names contrast_F.names ];
    contrast.values = [contrast_T.values contrast_F.values];
    contrast.types  = [contrast_T.types contrast_F.types];

%% job execute
     
    job_first_level_contrast(fspm{ir}, contrast, par);

end


%% S.OUARAB for Covid study
% out = load('out.mat')
% out = out.out;

% suj(contains(suj,out)) = []
% [pat, hv] = separate(suj);
% suj  = [pat; hv]
% dfunc = gdir(suj,'fMRI_RS_MNI');
% dreact = gdir(dfunc,'React');
% 
% 
% fcsv = gfile(dir,'info_covid_suj.csv$');
% info = readtable(fcsv{1});    
% sujname = info.suj;
% gender_str  = info.Gender;
% age     = str2num(cell2mat(info.age));   % Convert str to number 
% gender  = contains(gender_str, 'M');   % sum(gender) is the number of male and sum(~gender) is the number of female
% 
% 
% 
% % Ach
% glm_Ach = gdir(dreact,'GLM_ACH');
% 
% img1 = gfile(glm_Ach(1:length(pat)),'^beta.*1.*nii');   
% img2 = gfile(glm_Ach(length(pat)+1:end),'^beta.*1.*nii');  
% 
% 
% 
% 
% 
% 
% 
% dirStat = {'/network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/data/crc_covid/fMRI_stats/react/Ach_out7suj'}
% 
% clear par
% 
% %par.cov_name = {'age','gender'};
% %par.age = age;
% %par.gender = gender;
% 
% par.th_masking = 0;      % not use absolute threshold masking 
% par.use_imask  = 1;     % not use implicit mask
% 
% 
% par.use_gms  = 0;
% %par.normalisation = 2;   % Using proportional scaling
% par.rcontrast     = 0;   % don't result contrast for VBM
% 
% job_do_two_sample_ttest(dirStat,img1,img2,par)
% 
% 
% 
% 
% 
% % DA
% 
% 
% glm_Da = gdir(dreact,'GLM_DA');
% 
% img1 = gfile(glm_Da(1:length(pat)),'^beta.*1.*nii');   
% img2 = gfile(glm_Da(length(pat)+1:end),'^beta.*1.*nii');  
% 
% 
% 
% 
% 
% 
% 
% dirStat = {'/network/lustre/iss02/cenir/analyse/irm/users/salim.ouarab/data/crc_covid/fMRI_stats/react/Da_out7suj'}
% 
% clear par
% 
% %par.cov_name = {'age','gender'};
% %par.age = age;
% %par.gender = gender;
% 
% par.th_masking = 0;      % not use absolute threshold masking 
% par.use_imask  = 1;     % not use implicit mask
% 
% 
% par.use_gms  = 0;
% %par.normalisation = 2;   % Using proportional scaling
% par.rcontrast     = 0;   % don't result contrast for VBM
% 
% job_do_two_sample_ttest(dirStat,img1,img2,par)
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% info = gfile(dir,'info_maladie.csv')
% 
% out = load('out.mat')
% out = out.out;
% 
% suj(contains(suj,out)) = []
% 
% [pat, hv] = separate(suj);
% suj  = [pat; hv]
% dfunc = gdir(suj,'fMRI_RS_MNI');
% dreact = gdir(dfunc,'React');
% 
% 
% 
% 
% 
% data = readtable(info{1}) % supp FL05 et AL18
% sind = data.suj
% [~,name] = get_parent_path(pat,1);
% for i = 1:length(name)
%     for j = 1:length(sind)
%      if contains(name(i),sind(j))
%        delai(i) = data.delai(j)
%      end
%     end
% end
% 
% 
% 
% 
% 
% 
% cov2 = r_mkdir(vbm,'cov2')
% 
% clear par
% 
% par.name = 'score_ini_total'
% par.cov= Score_ini_total
% par.dir = cov2;
% par.age =AGE';
% par.sex= sex';
% %par.tiv = volumes';
% par.run=1;
% 
% 
% job=  multi_regression(fgr1,par)