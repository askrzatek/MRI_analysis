%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script for multiple regression analysis : AUDICOG + behav                                                  %%
%% Anna SKRZATEK January 2024                                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
load e

%% Initialise
addpath /home/anna.skrzatek/MRI_analysis
addpath /home/anna.skrzatek/MRI_analysis/AUDICOG

project_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG'
main_dir = fullfile(project_dir,'/DATA','/Non_chirurgicaux')

cd (main_dir)

sujDir = e.getPath();
Input_RS = e.getSerie('run_RS').getVolume('s5wts_OC').getPath();
e.addSerie('RS','tedana','tedana_RS')
ALFF_InputDir = get_subdir_regex(e.getSerie('tedana_RS').gpath(),'rsfc');
ALFF_Input = get_subdir_regex_files(ALFF_InputDir,'^ALFF_clean.nii');
%mkdir(project_dir,'/Results/fALFF');
ALFF_OutputDir = fullfile(project_dir,'Results/ALFF');


%% Regressors definition
% importing the table with multiple columns with patients characteristics
% from CSV: filtered by IRM==1
tab = readtable(fullfile(project_dir,'DATA/ANT_Alerting_RT_multiregression.csv'));
tab_sel = ismember(tab.IRM, 1);
group = tab.Group(tab_sel);
age   = tab.Age(tab_sel);
sex   = tab.Genre(tab_sel);
ANT   = tab.log_ANT_RT_Alerting(tab_sel);

%% Creating the par structure for multiple regression (both groups combined)
clear par
par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.jobname = 'ANT_ALFF_reg_model_spec';
par.mem = 10000;
par.nb_cond = 1;
par.nb_cons = length(ALFF_Input);
target_regressors.name  = {'Alerting_ANT_RT_log'};
%target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1,targetUPDRSIII_Sup1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)
target_regressors.value = {ANT};
covars.name = {'Age';'Sex';'Group'};
covars.val = {age; sex; group};
%mkdir(ALFF_OutputDir,'Multireg_ANT');
outputdir = fullfile(ALFF_OutputDir,'Multireg_ANT');

varcov_multiregression_model_spec({cellstr(outputdir)}, {ALFF_Input(:)}, covars, target_regressors, par);

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

%path_masks_nback= '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/NBack/Masks_positive_effect_2back' ;
%path_masks_audio= '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/Audio/2023_01_27 - Audio' ;
%path_masks_yeo_aal3 = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/Yeo_networks_with_aal3_masks' ;
masks_CAREN_SN = get_subdir_regex_files('/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN','RSN01') ;
masks_CAREN_SN = cellstr(masks_CAREN_SN{:,1});

masks_CAREN_CEN = get_subdir_regex_files('/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN','RSN02') ;
masks_CAREN_CEN = cellstr(masks_CAREN_CEN{:,1});

masks_CAREN_DMN = get_subdir_regex_files('/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN','RSN04') ;
masks_CAREN_DMN = cellstr(masks_CAREN_DMN{:,1});

%% descriptions inside AAL correspondance table, once it is adapted for use and create abbrevs

[ids,labels,aal_caren_tab] = xlsread('/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/aal_id_labels_abbrevs_CAREN.xlsx');
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
