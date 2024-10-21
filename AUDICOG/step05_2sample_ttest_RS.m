%% paramètres / variables
clc
clear
addpath /home/anna.skrzatek/MRI_analysis/AUDICOG/

main_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/DATA/Non_chirurgicaux';
project_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';
rsfc = 1;
% ROIs = {'ParaHipp','BA_31', 'Orb_PFC', 'lAudio', 'rAudio', 'Cingulate'};
ROIs = {'ParaHipp','lAudio', 'rAudio', 'Cingulate', 'OFC'};
models.names = {    'rsfc_verif_ANT_RT_STD_pca_covariates',                     'rsfc_verif_Alert_pca_covariates',                                   'rsfc_verif_ANT_RT_STD_wo_pca','rsfc_verif_Alert_wo_pca',          'rsfc_verif_double_pca_covariates',                                               'rsfc_verif_double_wo_pca'};
models.covarnames = {{'Age','Genre','pca_audio1','pca_emotionnel1','STD_RT_ANT'},{'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting'},{'Age','Genre','STD_RT_ANT'},{'Age','Genre','log_ANT_RT_Alerting'},{'Age','Genre','pca_audio1','pca_emotionnel1','log_ANT_RT_Alerting','STD_RT_ANT'},{'Age','Genre','log_ANT_RT_Alerting','STD_RT_ANT'}};
    
cd (project_dir)

% outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/fALFF'} ;
% outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_BA31','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_ParaHipp'} ;
% outdirs = {'/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_BA31','/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/RS_2sample_ttest/ALFF_ParaHipp'} ;
for imodel = 1:length(models.names)
    model_outdir = fullfile(project_dir,'/Results/RS_2sample_ttest/',models.names{imodel});
%     mkdir(model_outdir)
    for ir = 1:length(ROIs)
        outdirs{ir} = fullfile(model_outdir,ROIs{ir});
%         mkdir(outdirs{ir});
    end
    models.outdirs{imodel} = outdirs;
end
clear imodel

load('e_nonchir.mat');
% fichier de correspondance numero IRM - comportement - groupe - age
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;

%% Regressors definition
% importing the table with multiple columns with patients characteristics

tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv')); % from CSV: filtered by IRM==1 - definietely some errors in it
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2.csv')); % behavioral data verified and up to date for group 1 & 2

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

%% Choose  the model you want to apply : accordingly to the desired number of covariates
imodel = 5; % 1:rsfc_verif_ANT_RT_STD_pca_covariates 2:rsfc_verif_Alert_pca_covariates 3:rsfc_verif_ANT_RT_STD_wo_pca 4:rsfc_verif_Alert_wo_pca 5:rsfc_verif_double_pca_covariates 6:rsfc_verif_double_wo_pca
sprintf('Model %s chosen',models.names{imodel})
sprintf('Covariates to be used are %s %s %s %s %s %s', models.covarnames{imodel}{:})


% covars.name = {'Age','Genre','pca_audio1','pca_emotion1','ANT_Alert_Score'};
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
            scans{igroup}.cov{j,2} = str2double(tab.Genre(tab.code_IRM == id));
            
            clear par
            par.covars = 1;
            if imodel == 1 || imodel ==2 || imodel == 5 % PCA covars included
                scans{igroup}.cov{j,3} = tab.pca_audio1(tab.code_IRM == id);
                scans{igroup}.cov{j,4} = tab.pca_emotion1(tab.code_IRM == id);
                if imodel == 1
                    scans{igroup}.cov{j,5} = tab.RT_STD(tab.code_IRM == id);
                    par.intercov = [1,1,1,1,2];
                elseif imodel == 2
                    scans{igroup}.cov{j,5} = tab.ANT_Alert_Score(tab.code_IRM == id);
                    par.intercov = [1,1,1,1,2];
                elseif imodel == 5
                    scans{igroup}.cov{j,5} = tab.ANT_Alert_Score(tab.code_IRM == id);
                    scans{igroup}.cov{j,6} = tab.RT_STD(tab.code_IRM == id);
                    par.intercov = [1,1,1,1,2,2];
                end
            elseif imodel == 3 % No PCA covars included
                scans{igroup}.cov{j,3} = tab.RT_STD(tab.code_IRM == id);
                par.intercov = [1,1,2];
            elseif imodel == 4
                scans{igroup}.cov{j,3} = tab.ANT_Alert_Score(tab.code_IRM == id);
                par.intercov = [1,1,2];
            else
                scans{igroup}.cov{j,3} = tab.ANT_Alert_Score(tab.code_IRM == id);
                scans{igroup}.cov{j,4} = tab.RT_STD(tab.code_IRM == id);
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

%% MODEL ESTIMATION PER CONTRAST
% fspm = addsuffixtofilenames(outdirs(1),'/SPM.mat'); % BA31
% fspm = addsuffixtofilenames(outdirs(2),'/SPM.mat'); % ParaHipp
for imodel = 1: length(models.names)
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
%imodel = 2;

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
for imodel = 3:6
    % Single variable of interest
    % F-stat
    switch(imodel) % 1:rsfc_verif_ANT_RT_STD_pca_covariates 2:rsfc_verif_Alert_pca_covariates 3:rsfc_verif_ANT_RT_STD_wo_pca 4:rsfc_verif_Alert_wo_pca 5:rsfc_verif_double_pca_covariates 6:rsfc_verif_double_wo_pca
        case 1 % with PCAs
            clear contrast_F
            clear contrast
            disp('Single variable models with PCAs')
            Main_seed_effect = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];

            Main_RT_STD_effect = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 0 0 1 -1];
            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion = [0 0 0 0 0 1 0 0];
            Main_effect_Age = [0 0 1 0 0 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0 0 0];

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

        case 2
            clear contrast_F
            clear contrast
            Main_seed_effect = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            Main_Alert_effect = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Differential_group_Alert_effect = [0 0 0 0 0 0 1 -1];
            Main_effect_Hearing_Loss = [0 0 0 0 1 0 0 0];
            Main_effect_Emotion = [0 0 0 0 0 1 0 0];
            Main_effect_Age = [0 0 1 0 0 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0 0 0];

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

        case 3 % Without PCAs
            clear contrast_F
            clear contrast
            disp('Single variable models without PCAs')
            Main_seed_effect = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            Main_RT_STD_effect = [0 0 0 0 1 0 ; 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 1 -1];        
            Main_effect_Age = [0 0 1 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0];

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

        case 4
            clear contrast_F
            clear contrast
            Main_seed_effect = [1 0 0 0 0 0 ; 0 1 0 0 0 0];
            Main_Alert_effect = [0 0 0 0 1 0 ; 0 0 0 0 0 1];
            Differential_group_Alert_effect = [0 0 0 0 1 -1];        
            Main_effect_Age = [0 0 1 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0];

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

        case 5
            clear contrast_F
            clear contrast
            disp('Double variable models') % with PCA

            % Double variable of interest
            Main_seed_effect = [1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            Main_Alert_effect = [0 0 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 0 0];
            Differential_group_Alert_effect = [0 0 0 0 0 0 1 -1 0 0];
            Main_RT_STD_effect = [0 0 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 0 0 0 0 1 -1];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 0 0];
            Main_effect_Age = [0 0 1 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0];

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


        case 6
            clear contrast_F
            clear contrast
            %without PCA
            Main_seed_effect = [1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            Main_Alert_effect = [0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0];
            Differential_group_Alert_effect = [0 0 0 0 1 -1 0 0];
            Main_RT_STD_effect = [0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1];
            Differential_group_RT_STD_effect = [0 0 0 0 0 0 1 -1];
            %Control_Tinnitus = [1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0];
            %Tinnitus_Control = [-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0];
            Main_effect_Age = [0 0 1 0 0 0];
            Main_effect_Sex = [0 0 0 1 0 0];

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
