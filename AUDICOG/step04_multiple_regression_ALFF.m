%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script for multiple regression analysis : AUDICOG + behav                                                  %%
%% Anna SKRZATEK January 2024                                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

%% Initialise
addpath /home/anna.skrzatek/MRI_analysis
addpath /home/anna.skrzatek/MRI_analysis/AUDICOG

project_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG'
main_dir = fullfile(project_dir,'/DATA','/Non_chirurgicaux')

cd (main_dir)
load e

sujDir = e.getPath();
Input_RS = e.getSerie('run_RS').getVolume('s5wts_OC').getPath();
e.addSerie('RS','tedana','tedana_RS')

%% Regressors definition
% importing the table with multiple columns with patients characteristics
% from CSV: filtered by IRM==1
clear par
par.rsfc = 0;
par.ALFF = 1;

if par.rsfc ==1
    seeds = {'ParaHipp','lAudio','rAudio','Cingulate','OFC'};
    pathway_contrasts = {'/tedana009a1_vt/rsfc/'};
    contrast_names = {'seed2voxel_pearson_ParaHipp_lAudio_rAudio__Cingulate_OFC'};
    ncon = length(seeds);
    % contrast_names = {'seed2voxel_pearson_ParaHipp_BA_31_Orb_PFC_lAudio_rAudio__Cingulate'};
end

cd (project_dir)
tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv'));
tab = readtable(fullfile(project_dir,'DATA/AUDICOG_behavioral_data_Groups1_2.csv')); % behavioral data verified and up to date for group 1 & 2
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;


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
            scans{igroup}.cov{j,2} = tab.Genre(tab.code_IRM == id);
            scans{igroup}.cov{j,3} = str2double(tab.Group(tab.code_IRM == id));
            scans{igroup}.cov{j,4} = tab.pca_audio1(tab.code_IRM == id);
            scans{igroup}.cov{j,5} = tab.pca_emotion1(tab.code_IRM == id);
            scans{igroup}.cov{j,6} = tab.ANT_Alert_Score(tab.code_IRM == id);
            scans{igroup}.cov{j,7} = tab.RT_STD(tab.code_IRM == id);
            if par.rsfc == 1
                pearson_map_rsfc = get_subdir_regex_files( fullfile(e(iSubj).getSerie('run_RS').path, pathway_contrasts{1}), contrast_names{1}) ; % if multiple pathway_contrasts then change 1 to icontr
                for icon = 1: length(seeds)
                    scans{igroup}.contrast{j, icon } = [char(pearson_map_rsfc), ',' num2str(icon) ] ;
                end
            elseif par.ALFF ==1
                ALFF_InputDir = get_subdir_regex(e(iSubj).getSerie('tedana_RS').gpath(),'rsfc');
                scans{igroup}.contrast{j, 1 } = get_subdir_regex_files(ALFF_InputDir,'^ALFF_clean.nii');
            end
        end
    end
end

% tab_sel = ismember(tab.IRM, 1);

age   = vertcat(scans{1}.cov{:,1},scans{2}.cov{:,1});
sex   = vertcat(scans{1}.cov{:,2},scans{2}.cov{:,2});
group = vertcat(scans{1}.cov{:,3},scans{2}.cov{:,3});
Audio1   = vertcat(scans{1}.cov{:,4},scans{2}.cov{:,4});
Emo1   = vertcat(scans{1}.cov{:,5},scans{2}.cov{:,5});
Alert_ANT   = vertcat(scans{1}.cov{:,6},scans{2}.cov{:,6});
RT_STD_ANT   = vertcat(scans{1}.cov{:,7},scans{2}.cov{:,7});

%% Creating the par structure for multiple regression (both groups combined)
clear par

target_regressors.name  = {'Alerting_Score_ANT','RT_STD_ANT','Hearing_Loss','Emotion'};
target_regressors.value = {Alert_ANT,RT_STD_ANT,Audio1,Emo1};
covars.name = {'Age';'Sex';'Group'};
covars.val = {age; sex; group};

par.nb_cons = length(target_regressors.name);

par.ALFF = 1;
par.rsfc = 0;
if par.ALFF == 1
    par.nb_cond = 1;
    par.jobname = 'ANT_ALFF_reg_model_spec';
    input = {scans{1}.contrast{:,1},scans{2}.contrast{:,1}};
elseif par.rsfc == 1
    par.nb_cond = length(seeds);
    par.jobname = 'rsfc_reg_model_spec';
    input = {scans{1}.contrast,scans{2}.contrast};
end
par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.mem = 10000;


mkdir(project_dir,'Results/Multireg')
parent_dir = fullfile(project_dir,'Results/Multireg');
for ireg = 1: par.nb_cons
    mkdir(parent_dir,target_regressors.name{ireg})
    for ir = 1: par.nb_cond
        if par.rsfc == 1
            mkdir(parent_dir,sprintf('%s/%s',target_regressors.name{ireg},seeds{ir}))
            outputdir{ireg}{ir} = fullfile(parent_dir,sprintf('%s/%s',target_regressors.name{ireg},seeds{ir}));
    %     mkdir(ALFF_OutputDir,sprintf('Multireg_%s',target_regressors.name{ireg}))
        elseif par.ALFF == 1
            mkdir(parent_dir,sprintf('%s/ALFF',target_regressors.name{ireg}))
            outputdir{ireg}{ir} = fullfile(parent_dir,sprintf('%s/ALFF',target_regressors.name{ireg}));
        end
    end
end

varcov_multiregression_model_spec(outputdir, input , covars, target_regressors, par);

%% Model estimation
fspm = addsuffixtofilenames(outputdir,'/SPM.mat');

clear par
par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.jobname = 'ALFF_ANT_multireg_est';
par.mem = 10000;

par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{1});

job_first_level_estimate({fspm},par)

%% Contrast creation for each SPM.mat
% t-statistics : vérifier les
    PosCorrelation = [0 0 0 0 1] ;
    NegCorrelation = [0 0 0 0 -1] ;

% fspm doesn't change: no need to load
%% Contrast names
    contrast_T.names = {
    sprintf('Pos correlation_%s_on_fALFF',target_regressors.name{1})
    sprintf('Neg correlation_%s_on_fALFF',target_regressors.name{1})}';

%% Contrast values
contrast_T.values = {
    PosCorrelation
    NegCorrelation}';

%% Contrast type
contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast.names  = [contrast_T.names];
contrast.values = [contrast_T.values];
contrast.types  = [contrast_T.types];

%% Contrast : write
clear par

par.sge = 0;
par.run = 1;
par.display = 0;
par.jobname = sprintf('spm_write_%s_con',target_regressors.name{1});

% par.sessrep = 'both';
par.sessrep = 'none';

par.delete_previous = 1;
par.report          = 0;

job_first_level_contrast({fspm},contrast,par);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF connectivity has been done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF connectivity hasn't been done yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list of network ROIs according to Fan et al. 2007 (MEG analyse)
AlertN = {'Temp_Sup_R', 'Pariet_Lob_Inf_L', 'Fusi_L, Front_Inf_L', 'Pariet_Lob_Sup_L'};
OrientN = {'Fusi_L','Fusi_R','Precent_L','Pariet_Lob_Sup_R','Front_Sup_L','Pariet_Lob_Sup_L','Postcent_R'};
ExecN = {'Front_Sup_L','Front_Inf_R','Fusi_L','Front_Inf_L','Front_Mid_R','Fusi_R','Cingul_Ant_R'};


%% Faire tourner les calculs de connnectivité

% define input volumes and confounds
clear par
par.run = 1;
par.redo = 0;
par.volume   =  e.getSerie('run_RS').getVolume('s5wts_OC') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound =  e.getSerie('run_RS').getRP('multiple_regressors');
par.mask_threshold = 0.001;

% define ROI, using several methods
par.roi_type.atlas_cat12 = 'aal3';

%path_masks_nback= '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/NBack/Masks_positive_effect_2back' ;
%path_masks_audio= '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/Audio/2023_01_27 - Audio' ;
%path_masks_yeo_aal3 = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/Yeo_networks_with_aal3_masks' ;
masks_CAREN_SN = get_subdir_regex_files('/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN','RSN01') ;
masks_CAREN_SN = cellstr(masks_CAREN_SN{:,1});

masks_CAREN_CEN = get_subdir_regex_files('/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN','RSN02') ;
masks_CAREN_CEN = cellstr(masks_CAREN_CEN{:,1});

masks_CAREN_DMN = get_subdir_regex_files('/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN','RSN04') ;
masks_CAREN_DMN = cellstr(masks_CAREN_DMN{:,1});

%% descriptions inside AAL correspondance table, once it is adapted for use and create abbrevs

[ids,labels,aal_caren_tab] = xlsread('/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/aal_id_labels_abbrevs_CAREN.xlsx');
abbrev      = {aal_caren_tab{2:length(aal_caren_tab),3}};
description = {aal_caren_tab{2:length(aal_caren_tab),2}};
%SN_labels_id = ids(1:length(masks_CAREN_SN),4);
%CEN_labels_id = ids(1:length(masks_CAREN_CEN),5);
%DMN_labels_id = ids(1:length(masks_CAREN_DMN),7);

%% attention à l'ordre des fichiers par rapport à l'ordre des numéros des labels
%% quand les paths sont sorted, alors Reg150 se trouve juste après le Reg15
%%

par.roi_type.mask_global = {
%     % path                                                 abbrev        description
        fullfile(masks_CAREN_SN(1:length(masks_CAREN_SN))), 'l_precentral_SN', 'Left_Precentral_Cortex_Salience' 
        % et ainsi de suite
        } ;


%% perform the extraction

TS = job_extract_timeseries(par);

%%
% % define some networks : not mandatory

par.network.salience = par.roi_type.mask_global{1:lenght(masks_CAREN_SN),2} ;
par.network.cen      = par.roi_type.mask_global{length(masks_CAREN_SN)+1:length(masks_CAREN_SN)+length(masks_CAREN_CEN),2} ;
par.network.dmn      = par.roi_type.mask_global{length(masks_CAREN_SN)+length(masks_CAREN_CEN)+1:length(par.roi_type.mask_global),2} ;


% perform connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

%% plot
guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% perfrom seed-to-voxel correlation (seed == ROI)
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 
