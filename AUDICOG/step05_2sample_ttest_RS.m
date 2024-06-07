%% paramètres / variables
clc
clear
addpath /home/anna.skrzatek/MRI_analysis/AUDICOG/

main_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
project_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';
rsfc = 1;
% ROIs = {'ParaHipp','BA_31', 'Orb_PFC', 'lAudio', 'rAudio', 'Cingulate'};
ROIs = {'ParaHipp','lAudio', 'rAudio', 'Cingulate', 'OFC'};
cd (project_dir)

% outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/fALFF'} ;
% outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_BA31','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_ParaHipp'} ;
% outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_BA31','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_ParaHipp'} ;

for ir = 1:length(ROIs)
    outdirs{ir} = fullfile(project_dir,'/Results/RS_2sample_ttest/rsfc_ANT_RT_STD_pca_covariates',ROIs{ir});
    outdirs{ir} = fullfile(project_dir,'/Results/RS_2sample_ttest/rsfc_Alert_pca_covariates',ROIs{ir});
%     outdirs{ir} = fullfile(project_dir,'/Results/RS_2sample_ttest/rsfc_double_pca_covariates',ROIs{ir});
    mkdir(outdirs{ir});
end

load('e_nonchir.mat');
% fichier de correspondance numero IRM - comportement - groupe - age
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;

%% Regressors definition
% importing the table with multiple columns with patients characteristics
% from CSV: filtered by IRM==1
tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv'));
%tab_sel = ismember(tab.IRM, 1);
%group = tab.Group(tab_sel);
%age   = tab.Age(tab_sel);
%sex   = tab.Genre(tab_sel);
%ANT   = tab.log_ANT_RT_Alerting(tab_sel);
%ANT_std = tab.ANT_STD_mean(tab_sel);
%hearing_loss = tab.pca_audio1(tab_sel);
%emotion = str2double(tab.pca_emotionnel1(tab_sel));

%% Inputs

% pathway_contrasts = '/tedana009a1_vt/rsfc/';
% contrast_names   = { 'ALFF_clean.nii'  'fALFF_clean.nii' } ;

% pathway_contrasts = {'/model/ALFF_VOIreg_BA31','/model/ALFF_VOIreg_ParaHipp'};
% contrast_names = {'ess_0001.nii','ess_0001.nii'};

pathway_contrasts = {'/tedana009a1_vt/rsfc/'};
% contrast_names = {'seed2voxel_pearson_ParaHipp_BA_31_Orb_PFC_lAudio_rAudio__Cingulate'};
contrast_names = {'seed2voxel_pearson_ParaHipp_lAudio_rAudio__Cingulate_OFC'};
% rsfc = 1;
ncon = length(ROIs);

groups.name = {'Control' 'Tinnitus'};
groups.val = cell(size(ROIs));
% covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting'};
% covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','STD_RT_ANT'};
% covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting','STD_RT_ANT'};

% getting scans & covariates organised per group
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
%                 scans{igroup}.cov{j,6} = tab.ANT_STD_mean(tab.code_IRM == id);
                scans{igroup}.cov{j,5} = tab.ANT_STD_mean(tab.code_IRM == id);
                
            end
            if rsfc == 1
                pearson_map_rsfc = get_subdir_regex_files( fullfile(e(iSubj).getSerie('run_RS').path, pathway_contrasts{1}), contrast_names{1}) ; % if multiple pathway_contrasts then change 1 to icontr
                for icon = 1: ncon
                    scans{igroup}.contrast{j, icon } = [char(pearson_map_rsfc), ',' num2str(icon) ] ;
                end
            else
                for icontr = 1:length(contrast_names) % missing RS folder path
                    scans{igroup}.contrast{j, icontr } = fullfile( e(iSubj).getSerie('run_RS').path, pathway_contrasts{icontr}, contrast_names{icontr}) ; % if multiple pathway_contrasts then change 1 to icontr
                end
            end
        end
    end
end

covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5})};
% covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5});vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6})};

%% ANALYSIS PER CONTRAST

if rsfc == 1
    for ir = 1 : length(ROIs)
        groups.val{ir} = {scans{1}.contrast(:,ir), scans{2}.contrast(:,ir)};
%         groups.val = {scans{1}.contrast(:,1), scans{2}.contrast(:,1)}; % ParaHipp
%         groups.val = {scans{1}.contrast(:,2), scans{2}.contrast(:,2)}; % BA31
%         groups.val = {scans{1}.contrast(:,1), scans{2}.contrast(:,3)}; % LC
%         groups.val = {scans{1}.contrast(:,2), scans{2}.contrast(:,3)}; % Frontal Orb
    end
else
    groups.val = {scans{1}.contrast(:,1), scans{2}.contrast(:,1)}; % BA31
    groups.val = {scans{1}.contrast(:,2), scans{2}.contrast(:,2)}; % ParaHipp
end


%% Model specify
clear par
par.sge = 1;
par.run = 0;
par.covars = 1;
par.intercov = [1,1,1,1,2];
% par.intercov = [1,1,1,1,2,2];

%% ANALYSIS PER CONTRAST
if rsfc == 1
    for ir = 1:length(ROIs)
%         par.jobname = sprintf('2sample_ttest_%s',ROIs{ir});
        par.jobname = '2sample_ttest_RS_double_model';
        par.jobname = '2sample_ttest_RS_Alert_pca_model';
        par.jobname = '2sample_ttest_RS_RT_STD_pca_model';
        varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(ir),covars,par);
    end
%     varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(1),covars,par) % ParaHipp
%     varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(2),covars,par) % BA31
%     varcov_2nd_level_2sample_model_spec(groups.val,outdirs(3),covars,par) % LC
%     varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(3),covars,par) % Frontal Orb
else
    varcov_2nd_level_2sample_model_spec(groups.val,outdirs(1),covars,par) % BA31
    varcov_2nd_level_2sample_model_spec(groups.val,outdirs(2),covars,par) % ParaHipp
end

%% ANALYSIS PER CONTRAST
% fspm = addsuffixtofilenames(outdirs(1),'/SPM.mat'); % BA31
% fspm = addsuffixtofilenames(outdirs(2),'/SPM.mat'); % ParaHipp
for ir = 1 : length(ROIs)
    fspm{ir} = addsuffixtofilenames(outdirs(ir),'/SPM.mat');
    
    %% Model estimate
    clear par
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = sprintf('spm_est_Alert');
    par.jobname  = sprintf('spm_est_RT_STD');
    job_first_level_estimate(fspm{ir},par)

end


%% Contrast definition

%% Single variable of interest
% F-stat
Main_seed_effect = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
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

%% Double variable of interest
% Main_seed_effect = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
% Main_group_effect_Alert = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
% Differential_group_Alert_effect = [0 0 0 0 0 0 1 -1 0 0];
% Main_group_effect_RT_STD = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
% Differential_group_RT_STD_effect = [0 0 0 0 0 0 0 0 1 -1];
% Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0];
% Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];


% t-stat
Control = [1 0];
Tinnitus =[0 1];

% contrast_F.names = {
%     'Main_ANT_effect'
%     'Main_ANT_effect_on_Group_ALFF'
%     }';
% contrast_F.names = {
%     'Main_ANT_effect'
%     'Main_ANT_effect_on_Group_fALFF'
%     }';

contrast_F.names = {
    'Main_Seed_Effect'
    'Main_Group_effect_on_Alert_Score'
    'Differential_Group_effect_on_Alert_Score'
    'Main_effect_Hearing_Loss'
    'Main_effect_Emotion'
    'Main_effect_Age'
    'Main_effect_Sex'
%     'Main_Group_effect_on_RT_STD'
%     'Differential_Group_effect_on_RT_STD'
%     'Control-Tinnitus'
%     'Tinnitus-Control'
    }';
contrast_T.names = {
    'Control-Tinnitus'
    'Tinnitus-Control'
    }';

contrast_F.values = {
    Main_seed_effect
    Main_group_effect_Alert
    Differential_group_Alert_effect
    Main_effect_Hearing_Loss
    Main_effect_Emotion
    Main_effect_Age
    Main_effect_Sex
%     Main_group_effect_RT_STD
%     Differential_group_RT_STD_effect
%     Control_Tinnitus
%     Tinnitus_Control
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

%% Contrast : write
for ir = 1 : length(ROIs)
    clear par

    par.sge = 0;
    par.run = 1;
    par.display = 0;
    par.jobname = 'spm_write_con';

    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;
     
    job_first_level_contrast(fspm{ir}, contrast, par);

end

