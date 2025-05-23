%% paramètres / variables
clc
clear
addpath /home/anna.skrzatek/MRI_analysis/AUDICOG/

main_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
project_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG';
rsfc = 1;
% ROIs = {'ParaHipp','BA_31', 'Orb_PFC', 'lAudio', 'rAudio', 'Cingulate'};
% ROIs = {'ParaHipp', 'ACC', 'lAudio', 'rAudio', 'Cingulate', 'OFC'};
ROIs = {'ParaHipp_l','ParaHipp_r','ACC_l','ACC_r','A1_l','A1_r','Cingulate'}; %bilateral separately for more precision
% ROIs = {'ParaHipp', 'ACC', 'A1', 'Cingulate'}; %joint hemispheres for smaller number of multiple comparisons

%  1: rsfc_pca_audio2_pca_emo2_att1_att2 2: rsfc_pca_audio2 3: rsfc_pca_emo2_att1_att2 4: rsfc_pca_audio1_audio2 5: rsfc_pca_emo1_emo2_att1_att2 6:rsfc_pca_att1_att2_emo1_emo2 7:'rsfc_plsda_att1_att2_emo' 8:'rsfc_plsda_emo_att1_att2' 9:'rsfc_plsda_hearloss' 10:'rsfc_plsda_audio_emo_double_att_inter' 11:'rsfc_plsda_audio1_emo1_double_att' 12:'rsfc_plsda_audio_emo_double_cogni' 13:'rsfc_plsda_audio_emo_att_mem_inter' 14:'rsfc_plsda_audio_emo_att_mem' 15:'rsfc_plsda_audio_emo_cogni1_inter' 16:'rsfc_plsda_audio_emo_cogni1' 17:'rsfc_plsda_audio_emo_inter' 18:'rsfc_plsda_audio_emo' 19:'rsfc_verif_ANT_RT_STD_pca_audio1_emo1' 20:'rsfc_verif_Alert_pca_audio1_emo1' 21: 'rsfc_verif_ANT_RT_STD_wo_pca' 22: 'rsfc_verif_Alert_wo_pca' 23: 'rsfc_verif_double_Alert_STD_RT_pca_covariates'  24: 'rsfc_verif_double_Alert_STD_RT_wo_pca'rsfc_verif_Alert_pca_covariates' 21:'rsfc_verif_ANT_RT_STD_wo_pca' 22:'rsfc_verif_Alert_wo_pca' 23:'rsfc_verif_double_pca_covariates' 24:'rsfc_verif_double_wo_pca'
                                %1                                                          %2                          %3                                                      %4                                        %5                                                    %6                                                              %7                                                  %8                                                      %9                              %10                                                                         %11                                                                                 %12                                                                      %13                                                            %14                                                                         %15                                                        %16                                                          %17                                                     %18                                         %19                                                             %20                                                               %21                           %22                                     %23                                                                                     %24
models.names = {    'rsfc_pca_audio2_pca_emo2_att1_att2_no_interaction',           'rsfc_pca_audio2_pca_emo2_att1_att2',                           'rsfc_pca_audio2',          'rsfc_pca_emo2_att1_att2',                       'rsfc_pca_audio1_audio2',                  'rsfc_pca_emo1_emo2_att1_att2',                             'rsfc_pca_att1_att2_emo1_emo2',                             'rsfc_plsda_att1_att2_emo',                             'rsfc_plsda_emo_att1_att2',                             'rsfc_plsda_hearloss',           'rsfc_plsda__audio_emo_double_att_inter',                                           'rsfc_plsda_audio1_emo1_double_att',                                       'rsfc_plsda_audio_emo_double_cogni',                                         'rsfc_plsda_audio_emo_att_mem_inter',                             'rsfc_plsda_audio_emo_att_mem',                                                  'rsfc_plsda_audio_emo_cogni1_inter',                      'rsfc_plsda_audio_emo_cogni1',                  'rsfc_plsda_audio_emo_inter',                    'rsfc_plsda_audio_emo',                              'rsfc_verif_ANT_RT_STD_pca_audio1_emo1',                     'rsfc_verif_Alert_pca_audio1_emo1',                              'rsfc_verif_ANT_RT_STD_wo_pca',   'rsfc_verif_Alert_wo_pca',          'rsfc_verif_double_Alert_STD_RT_pca_covariates',                                     'rsfc_verif_double_Alert_STD_RT_wo_pca'};
models.covarnames = {{'Age','Genre','pca_audio2','pca_emo2','pca_att1','pca_att2'},{'Age','Genre','pca_audio2','pca_emo2','pca_att1','pca_att2'},{'Age','Genre','pca_audio2'},{'Age','Genre','pca_emo2','pca_att1','pca_att2'},{'Age','Genre','pca_audio1','pca_audio2'},{'Age','Genre','pca_emo1','pca_emo2','pca_att1','pca_att2'},{'Age','Genre','pca_att1','pca_att2','pca_emo1','pca_emo2'},{'Age','Genre','plsda_att1','plsda_att2','plsda_emo1'},{'Age','Genre','plsda_emo1','plsda_att1','plsda_att2'},{'Age','Genre','plsda_HFaudio1'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_att1','plsda_att2'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_att1','plsda_att2'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_cogni1','plsda_cogni2'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_att1','plsda_mem1'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_att1','plsda_mem1'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_cogni1'},{'Age','Genre','plsda_HFaudio1','plsda_emo1','plsda_cogni1'},{'Age','Genre','plsda_HFaudio1','plsda_emotionnel1'},{'Age','Genre','plsda_HFaudio1','plsda_emotionnel1'},{'Age','Genre','pca_audio1','pca_emotionnel1','STD_RT_ANT'},{'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting'},{'Age','Genre','STD_RT_ANT'},{'Age','Genre','log_ANT_RT_Alerting'},{'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting','STD_RT_ANT'},{'Age','Genre','log_ANT_RT_Alerting','STD_RT_ANT'}};
    
cd (project_dir)

% outdirs = {'/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF','/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/fALFF'} ;
% outdirs = {'/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_BA31','/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_ParaHipp'} ;
% outdirs = {'/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_BA31','/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_ParaHipp'} ;
for imodel = 1:length(models.names)
    model_outdir = fullfile(project_dir,'/Results/RS_2sample_ttest/RSFC_cogni_interpretable_models',models.names{imodel});
      mkdir(model_outdir)
    for ir = 1:length(ROIs)
        outdirs{ir} = fullfile(model_outdir,ROIs{ir});
      mkdir(outdirs{ir});
    end
    models.outdirs{imodel} = outdirs;
end
clear imodel

% load('e_nonchir.mat');
e = exam(main_dir, 'AUDICOG_Suj'); % all subjects with multi-echo
% e.addSerie('RS$','tedana', 'run_RS', 1 );
e.addSerie('RS$', 'run_RS', 1 );
% e.getSerie('run_RS').addVolume('^s5wts_OC.nii$','s5wts_OC',1);
% e.getSerie('run_RS').addRP('multiple_regressors','multiple_regressors',1)
e.explore


% fichier de correspondance numero IRM - comportement - groupe - age
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;

%% Regressors definition
% importing the table with multiple columns with patients characteristics

tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv')); % from CSV: filtered by IRM==1 - definietely some errors in it
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2.csv')); % behavioral data verified and up to date for group 1 & 2
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2-28-02-2025.csv')); % behavioral data verified and up to date for group 1 & 2
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2-27-03-2025.csv')); % behavioral data plsda + pca verified and up to date for group 1 & 2
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2-01-04-2025-scaled.csv')); % behavioral data plsda + pca verified, scaled and up to date for group 1 & 2
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2-18-04-2025-scaled.csv')); % behavioral data plsda + pca verified, scaled and up to date for group 1 & 2 + cognition variables reduced for interpretability

%tab_sel = ismember(tab.IRM, 1);
%group = tab.Group(tab_sel);
%age   = tab.Age(tab_sel);
%sex   = tab.Genre(tab_sel);
%ANT   = tab.log_ANT_RT_Alerting(tab_sel);
%ANT_std = tab.ANT_STD_mean(tab_sel);
%hearing_loss = tab.pca_audio1(tab_sel);
%emotion = str2double(tab.pca_emotionnel1(tab_sel));
%hearing_loss = tab.plsda_audio1(tab_sel);
%emotion = str2double(tab.plsda_emo1(tab_sel));

%% Inputs

% pathway_contrasts = '/tedana009a1_vt/rsfc/';
% contrast_names   = { 'ALFF_clean.nii'  'fALFF_clean.nii' } ;

% pathway_contrasts = {'/model/ALFF_VOIreg_BA31','/model/ALFF_VOIreg_ParaHipp'};
% contrast_names = {'ess_0001.nii','ess_0001.nii'};

pathway_contrasts = {'/tedana009a1_vt/rsfc/'};
% contrast_names = {'seed2voxel_pearson_ParaHipp_BA_31_Orb_PFC_lAudio_rAudio__Cingulate'};
% contrast_names = {'seed2voxel_pearson_ParaHipp_ACC_lAudio_rAudio__Cingulate_OFC'};
contrast_names = {'seed2voxel_pearson_ParaHipp_l_ParaHipp_r_ACC_l_ACC_r_A1_l_A1_r__Cingulate'};
% contrast_names = {'seed2voxel_pearson_ParaHipp_ACC_A1__Cingulate'};

% rsfc = 1;
ncon = length(ROIs);

groups.name = {'Control' 'Tinnitus'};
groups.val = cell(size(ROIs));

%% Choose  the model you want to apply : accordingly to the desired number of covariates
imodel = 1; %  %  %  1:rsfc_pca_audio2_pca_emo2_att1_att2_no_interaction 1: rsfc_pca_audio2_pca_emo2_att1_att2 2: rsfc_pca_audio2 3: rsfc_pca_emo2_att1_att2 4: rsfc_pca_audio1_audio2 5: rsfc_pca_emo1_emo2_att1_att2 6:rsfc_pca_att1_att2_emo1_emo2 7:'rsfc_plsda_att1_att2_emo' 8:'rsfc_plsda_emo_att1_att2' 9:'rsfc_plsda_hearloss' 10:'rsfc_plsda_audio_emo_double_att_inter' 11:'rsfc_plsda_audio1_emo1_double_att' 12:'rsfc_plsda_audio_emo_double_cogni' 13:'rsfc_plsda_audio_emo_att_mem_inter' 14:'rsfc_plsda_audio_emo_att_mem' 15:'rsfc_plsda_audio_emo_cogni1_inter' 16:'rsfc_plsda_audio_emo_cogni1' 17:'rsfc_plsda_audio_emo_inter' 18:'rsfc_plsda_audio_emo' 19:'rsfc_verif_ANT_RT_STD_pca_audio1_emo1' 20:'rsfc_verif_Alert_pca_audio1_emo1' 21: 'rsfc_verif_ANT_RT_STD_wo_pca' 22: 'rsfc_verif_Alert_wo_pca' 23: 'rsfc_verif_double_Alert_STD_RT_pca_covariates'  24: 'rsfc_verif_double_Alert_STD_RT_wo_pca'
for imodel = 1:24
    sprintf('Model %s chosen',models.names{imodel})
    sprintf('Covariates to be used are %s %s %s %s %s %s', models.covarnames{imodel}{:})


    % covars.name = {'Age','Genre','pca_audio1','pca_emotion1','ANT_S_Alerting'};
    % covars.name = {'Age','Genre','pca_audio1','pca_emotionnel1','STD_RT_ANT'};

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

            %% Adding Covariates and defining the model
    %         for icov = 1:length(models.covarnames{imodel})
                scans{igroup}.cov{j,1} = tab.Age(tab.code_IRM == id);
    %             scans{igroup}.cov{j,2} = str2double(tab.Genre(tab.code_IRM == id));
                scans{igroup}.cov{j,2} = tab.Genre(tab.code_IRM == id);

                clear par
                par.covars = 1;

                if imodel == 1 || imodel ==2 || imodel ==3 || imodel == 4 || imodel == 5 || imodel == 6 || imodel == 7 || imodel == 8 || imodel == 9 || imodel == 10 || imodel == 11 || imodel == 12 || imodel == 13 || imodel == 14 || imodel == 15 || imodel == 16 || imodel == 17 || imodel == 18 || imodel == 19 || imodel == 20 || imodel == 23 % PCA covars included
    %                 scans{igroup}.cov{j,3} = tab.audio_pca1(tab.code_IRM == id); % only useful for the oldest models (15-20)
    %                 scans{igroup}.cov{j,4} = tab.emo_pca1(tab.code_IRM == id); % only useful for the oldest models (15-20)
                    if imodel == 1
                        scans{igroup}.cov{j,3} = tab.audio_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.emo_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.cogni_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.cogni_pca2(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,1,1];
                    elseif imodel == 2
                        scans{igroup}.cov{j,3} = tab.audio_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.emo_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.cogni_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.cogni_pca2(tab.code_IRM == id);
                        par.intercov = [1,1,2,2,2,2];
                    elseif imodel == 3
                        scans{igroup}.cov{j,3} = tab.audio_pca2(tab.code_IRM == id);
                        par.intercov = [1,1,2];
                    elseif imodel == 4
                        scans{igroup}.cov{j,3} = tab.emo_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.cogni_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.cogni_pca2(tab.code_IRM == id);
                        par.intercov = [1,1,2,2,2];
                    elseif imodel == 5
                        scans{igroup}.cov{j,3} = tab.audio_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.audio_pca2(tab.code_IRM == id);
                        par.intercov = [1,1,2,2];
                    elseif imodel == 6
                        scans{igroup}.cov{j,3} = tab.emo_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.emo_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.cogni_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.cogni_pca2(tab.code_IRM == id);                    
                        par.intercov = [1,1,2,2,2,2];
                    elseif imodel == 7
                        scans{igroup}.cov{j,3} = tab.cogni_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.cogni_pca2(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.emo_pca1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.emo_pca2(tab.code_IRM == id);                    
                        par.intercov = [1,1,2,2,2,2];
                    elseif imodel == 8
                        scans{igroup}.cov{j,3} = tab.plsda_att1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_att2(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_emo1(tab.code_IRM == id);                    
                        par.intercov = [1,1,2,2,2];
                    elseif imodel == 9
                        scans{igroup}.cov{j,3} = tab.plsda_emo1(tab.code_IRM == id);                    
                        scans{igroup}.cov{j,4} = tab.plsda_att1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_att2(tab.code_IRM == id);
                        par.intercov = [1,1,2,2,2];
                    elseif imodel == 10
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);                
                        par.intercov = [1,1,2];
                    elseif imodel == 11
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_att1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.plsda_att2(tab.code_IRM == id);
                        par.intercov = [1,1,2,2,2,2];
                    elseif imodel == 12
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_att1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.plsda_att2(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,1,1];
                    elseif imodel == 13
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_cogni1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.plsda_cogni2(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,1,1];
                    elseif imodel == 14
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_att1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.plsda_mem1(tab.code_IRM == id);
                        par.intercov = [2,2,1,2,2,1];
                    elseif imodel == 15
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_att1(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.plsda_mem1(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,1,1];
                    elseif imodel == 16
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_cogni1(tab.code_IRM == id);
                        par.intercov = [1,1,2,2,2];
                    elseif imodel == 17
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.plsda_cogni1(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,1];
                    elseif imodel == 18
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        par.intercov = [1,1,2,2];
                    elseif imodel == 19
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        par.intercov = [1,1,1,1];
                    elseif imodel == 20
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.ANT_STD_mean(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,2];
                    elseif imodel == 21
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.ANT_S_Alerting(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,2];
                    elseif imodel == 24
                        scans{igroup}.cov{j,3} = tab.plsda_HFaudio1(tab.code_IRM == id);
                        scans{igroup}.cov{j,4} = tab.plsda_emo1(tab.code_IRM == id);
                        scans{igroup}.cov{j,5} = tab.ANT_S_Alerting(tab.code_IRM == id);
                        scans{igroup}.cov{j,6} = tab.ANT_STD_mean(tab.code_IRM == id);
                        par.intercov = [1,1,1,1,2,2];
                    end
                elseif imodel == 22 % No PCA covars included
                    scans{igroup}.cov{j,3} = tab.ANT_STD_mean(tab.code_IRM == id);
                    par.intercov = [1,1,2];
                elseif imodel == 23
                    scans{igroup}.cov{j,3} = tab.ANT_S_Alerting(tab.code_IRM == id);
                    par.intercov = [1,1,2];
                else
                    scans{igroup}.cov{j,3} = tab.ANT_S_Alerting(tab.code_IRM == id);
                    scans{igroup}.cov{j,4} = tab.ANT_STD_mean(tab.code_IRM == id);
                    par.intercov = [1,1,2,2];
                end
    %         end

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

    for icov = 1:length(models.covarnames{imodel})
        models.covarvals{imodel}(icov) = {vertcat(scans{1}.cov{:,icov},scans{2}.cov{:,icov})};
    end

    % covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5})};
    % covars.val = {vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5});vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6})};

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


    %% Model specify per contrast
    % clear par
    par.sge = 0;
    par.run = 1;

    covars.name = models.covarnames{imodel};
    covars.val   = models.covarvals{imodel};

    % par.covars = 1; %% done before
    % par.intercov = [1,1,1,1,2]; %% done before in the loops corresponding to
    % each model
    % par.intercov = [1,1,1,1,2,2];

    if rsfc == 1
        for ir = 1:length(ROIs)
            par.jobname = sprintf('2sample_ttest_RS_%s_%s',models.names{imodel},ROIs{ir});
    %         par.jobname = sprintf('2sample_ttest_%s',ROIs{ir});
    %         par.jobname = '2sample_ttest_RS_double_model';
    %         par.jobname = '2sample_ttest_RS_Alert_pca_model';
    %         par.jobname = '2sample_ttest_RS_RT_STD_pca_model';
            varcov_2nd_level_2sample_model_spec(groups.val{ir},models.outdirs{imodel}(ir),covars,par);
        end
    %     varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(1),covars,par) % ParaHipp
    %     varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(2),covars,par) % BA31
    %     varcov_2nd_level_2sample_model_spec(groups.val,outdirs(3),covars,par) % LC
    %     varcov_2nd_level_2sample_model_spec(groups.val{ir},outdirs(3),covars,par) % Frontal Orb
    else
        varcov_2nd_level_2sample_model_spec(groups.val,outdirs(1),covars,par) % BA31
        varcov_2nd_level_2sample_model_spec(groups.val,outdirs(2),covars,par) % ParaHipp
    end
end
%% MODEL ESTIMATION PER CONTRAST
% fspm = addsuffixtofilenames(outdirs(1),'/SPM.mat'); % BA31
% fspm = addsuffixtofilenames(outdirs(2),'/SPM.mat'); % ParaHipp
for imodel = 1 : length(models.names)
    for ir = 1 : length(ROIs)
        fspm{ir} = addsuffixtofilenames(models.outdirs{imodel}(ir),'/SPM.mat');
    
        %% Model estimate
        clear par
        par.run = 1;
        par.sge = 0;
        par.sge_queu = 'normal,bigmem';
        par.jobname  = sprintf('spm_est_%s_%s_%s_%s_%s_%s',models.covarnames{imodel}{:});
        job_first_level_estimate(fspm{ir},par)
    end
end


%% Contrast definition
%  imodel = 1;

% t-stat
Control = [1 0];
Tinnitus =[0 1];

contrast_T.names = {
    'Control>Tinnitus'
    'Control<Tinnitus'
    }';

contrast_T.values = {
    Control-Tinnitus
    Tinnitus-Control
    }';

contrast_T.types = cat(1,repmat({'T'},[1,length(contrast_T.names)]));

for imodel = 1:24
    % Single variable of interest
    % F-stat
    switch(imodel) % %  %  1: rsfc_pca_audio2_pca_emo2_att1_att2 2: rsfc_pca_audio2 3: rsfc_pca_emo2_att1_att2 4: rsfc_pca_audio1_audio2 5: rsfc_pca_emo1_emo2_att1_att2 6:rsfc_pca_att1_att2_emo1_emo2 7:'rsfc_plsda_att1_att2_emo' 8:'rsfc_plsda_emo_att1_att2' 9:'rsfc_plsda_hearloss' 10:'rsfc_plsda_audio_emo_double_att_inter' 11:'rsfc_plsda_audio1_emo1_double_att' 12:'rsfc_plsda_audio_emo_double_cogni' 13:'rsfc_plsda_audio_emo_att_mem_inter' 14:'rsfc_plsda_audio_emo_att_mem' 15:'rsfc_plsda_audio_emo_cogni1_inter' 16:'rsfc_plsda_audio_emo_cogni1' 17:'rsfc_plsda_audio_emo_inter' 18:'rsfc_plsda_audio_emo' 19:'rsfc_verif_ANT_RT_STD_pca_audio1_emo1' 20:'rsfc_verif_Alert_pca_audio1_emo1' 21: 'rsfc_verif_ANT_RT_STD_wo_pca' 22: 'rsfc_verif_Alert_wo_pca' 23: 'rsfc_verif_double_Alert_STD_RT_pca_covariates'  24: 'rsfc_verif_double_Alert_STD_RT_wo_pca'
        case 1 % with scaled PCAs Audio2 & Emo2 & Att1 & Att2 wo-Interaction
            clear contrast_F
            clear contrast
            disp('Interest PCAs variables model : audio2, emo 2, double cogni')
            Main_seed_effect = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Audio2     = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion2   = [0 0 0 0 0 1 0 0];
            Main_effect_Attention1 = [0 0 0 0 0 0 1 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Audio2'
                'Main_effect_Emotion2'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Audio2
                Main_effect_Emotion2
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];
            
            
        case 2 % with scaled PCAs Audio2 & Emo2 & Att1 & Att2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest PCAs variables model : audio2, emo 2, double cogni')
            Main_seed_effect = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Audio2     = [0 0 0 0 1 0 0 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0 0 0];
            Main_effect_Emotion2   = [0 0 0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Attention1 = [0 0 0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0 0 0 0 0];
            
            Inter_group_Audio2     = [0 0 0 0 1 -1 0 0 0 0 0 0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Emotion2   = [0 0 0 0 0 0 1 -1 0 0 0 0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention1 = [0 0 0 0 0 0 0 0 1 -1 0 0 ]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Attention2 = [0 0 0 0 0 0 0 0 0 0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Audio2'
                'Main_effect_Emotion2'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Audio2'
                'Inter_Group_Emotion2'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Audio2
                Main_effect_Emotion2
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Audio2
                Inter_group_Emotion2
                Inter_group_Attention1
                Inter_group_Attention2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];
            
            
        case 3 % with PCA Audio2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest scaled PCAs variables model : double audio')
            Main_seed_effect = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Audio2 = [0 0 0 0 1 0 ; 0 0 0 0 0 1];
            Main_effect_Age = [0 0 1 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0];
            
            Inter_group_Audio2 = [0 0 0 0 1 -1]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            
            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Audio2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Audio2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Audio2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Audio2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];
            
        case 4 % with scaled PCAs Emo2 & Att1 & Att2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : double cogni')
            Main_seed_effect       = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Emotion2   = [0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Attention1 = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0 0 0];
            
            Inter_group_Emotion2   = [0 0 0 0 1 -1 0  0 0  0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention1 = [0 0 0 0 0  0 1 -1 0  0]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Attention2 = [0 0 0 0 0  0 0  0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Emotion2'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Emotion2'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Emotion2
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Emotion2
                Inter_group_Attention1
                Inter_group_Attention2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

        case 5 % with PCAs Audio1 & Audio2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest scaled PCAs variables model : double audio')
            Main_seed_effect   = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Audio1 = [0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0];
            Main_effect_Audio2 = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Main_effect_Age    = [0 0 1 0 0 0 0 0];
            Main_effect_Sex    = [0 0 0 1 0 0 0 0];
            
            Inter_group_Audio1 = [0 0 0 0 1 -1 0 0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Audio2 = [0 0 0 0 0 0 1 -1]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            
            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Audio1'
                'Main_effect_Audio2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Audio1'
                'Inter_Group_Audio2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Audio1
                Main_effect_Audio2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Audio1
                Inter_group_Audio2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

        case 6 % with scaled PCAs Emo1 & Emo2 & Att1 & Att2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest PCAs variables model : double emo, double cogni')
            Main_seed_effect       = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Emotion1   = [0 0 0 0 1 0 0 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0 0 0];
            Main_effect_Emotion2   = [0 0 0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Attention1 = [0 0 0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0 0 0 0 0];
            
            Inter_group_Emotion1 = [0 0 0 0 1 -1 0 0 0 0 0 0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Emotion2 = [0 0 0 0 0 0 1 -1 0 0 0 0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention1 = [0 0 0 0 0 0 0 0 1 -1 0 0 ]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Attention2 = [0 0 0 0 0 0 0 0 0 0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Emotion1'
                'Main_effect_Emotion2'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Emotion1'
                'Inter_Group_Emotion2'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Emotion1
                Main_effect_Emotion2
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Emotion1
                Inter_group_Emotion2
                Inter_group_Attention1
                Inter_group_Attention2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

        case 7 % with scaled PCAs Att1 & Att2 & Emo1 & Emo2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest PCAs variables model : double cogni, double emo')
            Main_seed_effect       = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Attention1 = [0 0 0 0 1 0 0 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0 0 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Emotion1   = [0 0 0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Emotion2   = [0 0 0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0 0 0 0 0];
            
            Inter_group_Attention1 = [0 0 0 0 1 -1 0  0 0  0 0  0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention2 = [0 0 0 0 0  0 1 -1 0  0 0  0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Emotion1   = [0 0 0 0 0  0 0  0 1 -1 0  0]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Emotion2   = [0 0 0 0 0  0 0  0 0  0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Emotion1'
                'Main_effect_Emotion2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                'Inter_Group_Emotion1'
                'Inter_Group_Emotion2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Emotion1
                Main_effect_Emotion2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Attention1
                Inter_group_Attention2
                Inter_group_Emotion1
                Inter_group_Emotion2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

        % 7:'rsfc_plsda_att1_att2_emo'
        case 8 % with PLSDAs Att1 & Att2 & Emo Interaction
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : double cogni, emo')
            Main_seed_effect       = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Attention1 = [0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Emotion    = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0 0 0];
            
            Inter_group_Attention1 = [0 0 0 0 1 -1 0  0 0  0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention2 = [0 0 0 0 0  0 1 -1 0  0]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Emotion    = [0 0 0 0 0  0 0  0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Emotion'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                'Inter_Group_Emotion'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Emotion
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Attention1
                Inter_group_Attention2
                Inter_group_Emotion
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %8:'rsfc_plsda_emo_att1_att2' 
        case 9 % with PLSDAs Emo & Att1 & Att2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : double cogni, emo')
            Main_seed_effect       = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Emotion    = [0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Attention1 = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Attention2 = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age        = [0 0 1 0 0 0 0 0 0 0];
            Main_effect_Sex        = [0 0 0 1 0 0 0 0 0 0];
            
            Inter_group_Emotion    = [0 0 0 0 1 -1 0  0 0  0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention1 = [0 0 0 0 0  0 1 -1 0  0]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Attention2 = [0 0 0 0 0  0 0  0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Emotion'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Emotion'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Emotion
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Emotion
                Inter_group_Attention1
                Inter_group_Attention2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];
        
        %9:'rsfc_plsda_hearloss' 
        case 10 % with PLSDAs Hearing Interaction
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : hearing')
            Main_seed_effect         = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 ; 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0];
            
            Inter_group_Hearing_Loss = [0 0 0 0 1 -1 ] ; % 0 0 0 0 0 -1 0 0 0 0 0 0];
            
            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Hearing_Loss'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Hearing_Loss
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];
        
        
        %10:'rsfc_plsda_audio_emo_double_att_inter' 
        case 11 % with PLSDAs Hearing & Emo & Att1 & Att2 Interaction
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : audio, emo, double cogni')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Attention1   = [0 0 0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Attention2   = [0 0 0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0 0 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0 0 0 0 0 0];
            
            Inter_group_Hearing_Loss = [0 0 0 0 1 -1 0  0 0  0 0  0] ; % 0 0 0 0 0 -1 0 0 0 0 0 0];
            Inter_group_Emotion      = [0 0 0 0 0  0 1 -1 0  0 0  0]; %0 0 0 0 0 0 0 -1 0 0 0 0];
            Inter_group_Attention1   = [0 0 0 0 0  0 0  0 1 -1 0  0]; % 0 0 0 0 0 0 0 0 0 -1 0 0];
            Inter_group_Attention2   = [0 0 0 0 0  0 0  0 0  0 1 -1];% ; 0 0 0 0 0 0 0 0 0 0 0 -1];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Hearing_Loss'
                'Inter_Group_Emotion'
                'Inter_Group_Attention1'
                'Inter_Group_Attention2'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Hearing_Loss
                Inter_group_Emotion
                Inter_group_Attention1
                Inter_group_Attention2
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            % 11:'rsfc_plsda_audio1_emo1_double_att' 
        case 12 % with PLSDAs Hearing & Emo & Att1 & Att2
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : audio, emo, double cogni no interaction')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 1 0 0];
            Main_effect_Attention1   = [0 0 0 0 0 0 1 0];
            Main_effect_Attention2   = [0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Attention1'
                'Main_effect_Attention2'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Attention1
                Main_effect_Attention2
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %12:'rsfc_plsda_audio_emo_double_cogni' 
        case 13 % with PLSDAs Hearing & Emo & Cogni1 & Cogni2
            clear contrast_F
            clear contrast
            disp('Interest PLS-DAs variables model : double cogni')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 1 0 0];
            Main_effect_Cognition1   = [0 0 0 0 0 0 1 0];
            Main_effect_Cognition2   = [0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Cognition1'
                'Main_effect_Cognition2'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Cognition1
                Main_effect_Cognition2
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %13:'rsfc_plsda_audio_emo_att_mem_inter' 
        case 14 % with PLSDAs Hearing & Emo & Att & Mem Interaction
            clear contrast_F
            clear contrast
            disp('PLSD-A HL & Emo & Att & Mem Interaction model')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 0 0 1 0 0 0 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 0 0 1 0 0 0 0 ; 0 0 0 0 0 0 0 0 1 0 0 0];
            Main_effect_Attention    = [0 0 0 0 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 0 0 0 0 1 0];
            Main_effect_Memory       = [0 0 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 0 0 0 0];
            Main_effect_Sex          = [0 0 0 0 1 0 0 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0 0 0];
            
            Inter_group_Emotion      = [0 0 0  0 0  0 0 1 -1 0  0 0] ; %0 0 0 0 0 0 0 0 -1 0 0 0];
            Inter_group_Attention    = [0 0 0  0 0  0 0 0  0 1 -1 0] ; %0 0 0 0 0 0 0 0 0 0 -1 0];
            Inter_group_Age          = [0 0 1 -1 0  0 0 0  0 0  0 0] ; %0 0 0 -1 0 0 0 0 0 0 0 0];
            Inter_group_Sex          = [0 0 0  0 1 -1 0 0  0 0  0 0] ; %0 0 0 0 0 -1 0 0 0 0 0 0];
            

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Attention'
                'Main_effect_Memory'
                'Main_effect_Age'
                'Main_effect_Sex'
                'Inter_Group_Emotion'
                'Inter_Group_Attention'
                'Inter_Group_Age'
                'Inter_Group_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Attention
                Main_effect_Memory
                Main_effect_Age
                Main_effect_Sex
                Inter_group_Emotion
                Inter_group_Attention
                Inter_group_Age
                Inter_group_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %14:'rsfc_plsda_audio_emo_att_mem' 
        case 15 % with PLSDAs Hearing & Emo & Att & Mem
            clear contrast_F
            clear contrast
            disp('No interest variable model with PLS-DAs')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 1 0 0];
            Main_effect_Attention    = [0 0 0 0 0 0 1 0];
            Main_effect_Memory       = [0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Attention'
                'Main_effect_Memory'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Attention
                Main_effect_Memory
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            % 15:'rsfc_plsda_audio_emo_cogni1_inter' 
        case 16 % with PLSDAs Hearing & Emo & Cogni Interaction
            clear contrast_F
            clear contrast
            disp('No interest variable model with PLS-DAs')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            Control_Tinnitus         = [1 0 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0 0];
            Tinnitus_Control         = [-1 0 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 0 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
            Main_effect_Cognition    = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0 0 0 0];

            Inter_Group_Audio        = [0 0 0 0 1 -1 0  0 0  0] ; %0 0 0 0 0 -1 0 0 0 0];
            Inter_Group_Emotion      = [0 0 0 0 0  0 1 -1 0  0] ; %0 0 0 0 0 0 0 -1 0 0];
            Inter_Group_Cognition    = [0 0 0 0 0  0 0  0 1 -1] ; %0 0 0 0 0 0 0 0 0 -1];
            
            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Cognition'
                'Inter_Group_Hearing_Loss'
                'Inter_Group_Emotion'
                'Inter_Group_Cognition'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Cognition
                Inter_Group_Audio
                Inter_Group_Emotion
                Inter_Group_Cognition
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %16:'rsfc_plsda_audio_emo_cogni1' 
        case 17 % with PLSDAs Hearing & Emo & Cogni
            clear contrast_F
            clear contrast
            disp('No interest variable model with PLS-DAs')
            Main_seed_effect         = [1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 1 0];
            Main_effect_Cognition    = [0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Cognition'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Cognition
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %17:'rsfc_plsda_audio_emo_inter' 
        case 18 % with PLSDAs & Interaction
            clear contrast_F
            clear contrast
            disp('Emo & Hearing as interest variables')
            Main_seed_effect         = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0];
            Main_effect_Emotion      = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0 0 0];
            
            Inter_Group_Hearing_Loss = [0 0 0 0 1 -1 0  0];
            Inter_Group_Emotion      = [0 0 0 0 0  0 1 -1];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Interaction_Group_x_Hearing_Loss'
                'Main_effect_Emotion'
                'Interaction_Group_x_Emotion'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Inter_Group_Hearing_Loss
                Main_effect_Emotion
                Inter_Group_Emotion
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %18:'rsfc_plsda_audio_emo' 
        case 19 % with PCAs
            clear contrast_F
            clear contrast
            disp('No interest variable model with PLS-DAs')
            Main_seed_effect         = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_effect_Hearing_Loss = [0 0 0 0 1 0];
            Main_effect_Emotion      = [0 0 0 0 0 1];
            Main_effect_Age          = [0 0 1 0 0 0];
            Main_effect_Sex          = [0 0 0 1 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %19:'rsfc_verif_ANT_RT_STD_pca_audio1_emo1' 
        case 20 % with PCAs
            clear contrast_F
            clear contrast
            disp('Single variable models with PCAs')
            Main_seed_effect                 = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_RT_STD_effect               = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 0 0 1 -1];
            Main_effect_Hearing_Loss         = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion              = [0 0 0 0 0 1 0 0];
            Main_effect_Age                  = [0 0 1 0 0 0 0 0];
            Main_effect_Sex                  = [0 0 0 1 0 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_RT_STD_Effect'
                'Main_Group_x_RT_STD_Interaction'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_RT_STD_effect
                Differential_group_RT_STD_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %20:'rsfc_verif_Alert_pca_audio1_emo1' 
        case 21
            clear contrast_F
            clear contrast
            Main_seed_effect                = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            Main_Alert_effect               = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Differential_group_Alert_effect = [0 0 0 0 0 0 1 -1];
            Main_effect_Hearing_Loss        = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion             = [0 0 0 0 0 1 0 0];
            Main_effect_Age                 = [0 0 1 0 0 0 0 0];
            Main_effect_Sex                 = [0 0 0 1 0 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_Alert_Score_Effect'
                'Main_Group_x_Alert_Score_Interaction'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_Alert_effect
                Differential_group_Alert_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %21: 'rsfc_verif_ANT_RT_STD_wo_pca' 
        case 22 % Without PCAs
            clear contrast_F
            clear contrast
            disp('Single variable models without PCAs')
            Main_seed_effect                 = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            Main_RT_STD_effect               = [0 0 0 0 1 0 ; 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 1 -1];        
            Main_effect_Age                  = [0 0 1 0 0 0];
            Main_effect_Sex                  = [0 0 0 1 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_RT_STD_Effect'
                'Main_Group_x_RT_STD_Interaction'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_RT_STD_effect
                Differential_group_RT_STD_effect
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %22: 'rsfc_verif_Alert_wo_pca' 
        case 23
            clear contrast_F
            clear contrast
            Main_seed_effect                = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            Main_Alert_effect               = [0 0 0 0 1 0 ; 0 0 0 0 0 1];
            Differential_group_Alert_effect = [0 0 0 0 1 -1];        
            Main_effect_Age                 = [0 0 1 0 0 0];
            Main_effect_Sex                 = [0 0 0 1 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_Alert_Score_Effect'
                'Main_Group_x_Alert_Score_Interaction'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_Alert_effect
                Differential_group_Alert_effect
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %23: 'rsfc_verif_double_Alert_STD_RT_pca_covariates'  
        case 24
            clear contrast_F
            clear contrast
            disp('Double variable models') % with PCA

            % Double variable of interest
            Main_seed_effect                 = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            Main_Alert_effect                = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
            Differential_group_Alert_effect  = [0 0 0 0 0 0 1 -1 0 0];
            Main_RT_STD_effect               = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 0 0 0 0 1 -1];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            Main_effect_Hearing_Loss         = [0 0 0 0 1 0 0 0 0 0];
            Main_effect_Emotion              = [0 0 0 0 0 1 0 0 0 0];
            Main_effect_Age                  = [0 0 1 0 0 0 0 0 0 0];
            Main_effect_Sex                  = [0 0 0 1 0 0 0 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_Alert_Score_Effect'
                'Main_Group_x_Alert_Score_Interaction'
                'Main_RT_STD_Effect'
                'Main_Group_x_RT_STD_Interaction'
                'Main_effect_Hearing_Loss'
                'Main_effect_Emotion'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_Alert_effect
                Differential_group_Alert_effect
                Main_RT_STD_effect
                Differential_group_RT_STD_effect
                Main_effect_Hearing_Loss
                Main_effect_Emotion
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];

            %24: 'rsfc_verif_double_Alert_STD_RT_wo_pca'
        case 25
            clear contrast_F
            clear contrast
            %without PCA
            Main_seed_effect                 = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            Main_Alert_effect                = [0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0];
            Differential_group_Alert_effect  = [0 0 0 0 1 -1 0 0];
            Main_RT_STD_effect               = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 0 0 1 -1];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            Main_effect_Age                  = [0 0 1 0 0 0 0 0];
            Main_effect_Sex                  = [0 0 0 1 0 0 0 0];

            contrast_F.names = {
                'Main_Seed_Effect'
                'Main_Alert_Score_Effect'
                'Main_Group_x_Alert_Score_Interaction'
                'Main_RT_STD_Effect'
                'Main_Group_x_RT_STD_Interaction'
                'Main_effect_Age'
                'Main_effect_Sex'
                }';

            contrast_F.values = {
                Main_seed_effect
                Main_Alert_effect
                Differential_group_Alert_effect
                Main_RT_STD_effect
                Differential_group_RT_STD_effect
                Main_effect_Age
                Main_effect_Sex
                }';
            contrast_F.types = cat(1,repmat({'F'},[1,length(contrast_F.names)]));

            contrast.names  = [contrast_T.names contrast_F.names ];
            contrast.values = [contrast_T.values contrast_F.values];
            contrast.types  = [contrast_T.types contrast_F.types];


    end


    % Contrast : write
    for ir = 1 : length(ROIs)
        clear par
        fspm{ir} = addsuffixtofilenames(models.outdirs{imodel}(ir),'/SPM.mat');

        par.sge = 0;
        par.run = 1;
        par.display = 0;
        par.jobname = 'spm_write_con';

        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(fspm{ir}, contrast, par);

    end
end
