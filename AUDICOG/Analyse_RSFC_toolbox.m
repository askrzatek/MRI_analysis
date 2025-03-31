%% connectivity
clear
clc
main_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG';
ROI_dir = fullfile(main_dir,'Networks_Masks');


cd (main_dir)
% load('e_nonchir.mat');
data_dir = fullfile(main_dir,'/DATA/Non_chirurgicaux')

e = exam(data_dir, 'AUDICOG_Suj'); % all subjects with multi-echo
e.addSerie('RS$','tedana', 'run_RS', 1 );
e.getSerie('run_RS').addVolume('^s5wts_OC.nii$','s5wts_OC',1);
e.getSerie('run_RS').addRP('multiple_regressors','multiple_regressors',1)
e.explore
%% Faire tourner les calculs de connnectivité

% define input volumes and confounds
clear par
par.run = 1;
par.sge = 0;
par.jobname = 'timeseries_extract_AUDICOG_with_LC';
par.display = 0;
par.redo = 0;
par.volume   =  e.getSerie('run_RS').getVolume('s5wts_OC') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound =  e.getSerie('run_RS').getRP('multiple_regressors');
par.mask_threshold = 0.001;

% define ROI, using several methods
%par.roi_type.atlas_cat12 = 'aal3';

% path_masks_nback= '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/NBack/Masks_positive_effect_2back' ;
path_masks_audio= '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/Audio/2023_01_27 - Audio' ;
% path_masks_yeo_aal3 = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/Yeo_networks_with_aal3_masks' ; 
path_masks_Tinnitus = get_subdir_regex(ROI_dir, 'Tinnitus_Meta');
path_masks_Alert    = get_subdir_regex(ROI_dir, 'Alerting_effect');

par.roi_type.mask_global = {
%     % path                                                           abbrev          description
       char(get_subdir_regex_files(path_masks_Tinnitus, 'ParaHipp')), 'ParaHipp',  'Bilateral_ParaHippocampus'
%        char(get_subdir_regex_files(path_masks_Tinnitus, 'BA_31')),    'BA_31',     'Brodmann_Area_31'
% NOT ENOUGH SIGNAL IN 3T BOLD ACQUISITIONS SO AFTER RESLICING THE MASK, THE
% MASK BECOMES EMPTY FOR MIDBRAIN STRUCTURES LIKE LC
%        char(get_subdir_regex_files(path_masks_Alert, 'LC_l')),          'LC_l',        'Left_Locus_Coereleus'
%        char(get_subdir_regex_files(path_masks_Alert, 'LC_r')),          'LC_r',        'Right_Locus_Coereleus'
%        char(get_subdir_regex_files(path_masks_Alert, 'LC.nii')),          'LC',        'Bilateral_Locus_Coereleus'
%        char(get_subdir_regex_files(path_masks_Alert, 'Frontal_Orb')),          'Orb_PFC',        'Combined_Orbital_Prefrontal_Cortex'

       fullfile( path_masks_audio,  'Left_Auditory_activation_fwe05_0v_50subj.nii'), 'lAudio'   ,  'Left_Loca_Audio_activation'   
       fullfile( path_masks_audio,  'Right_Auditory_activation_fwe05_0v_50subj.nii'), 'rAudio'   ,  'Right_Loca_Audio_activation'   
%      fullfile(path_masks_nback, 'L_Parietal_fwe05.nii') , 'l_parietal_nback', 'Left_activation_nback_parietal'
%      fullfile(path_masks_nback, 'R_Parietal_fwe05.nii') , 'r_parietal_nback', 'Right_activation_nback_parietal'
%      fullfile(path_masks_nback, 'L_MFG_fwe05.nii')      , 'l_MedialFG_nback', 'Left_activation_nback_frontal'
%      fullfile(path_masks_nback, 'R_MFG_fwe05.nii')      , 'r_MedialFG_nback', 'Right_activation_nback_frontal'
%      fullfile(path_masks_yeo_aal3, 'lAG.nii')           , 'lAG_yeo'             , 'Left_AG_Yeo'
%      fullfile(path_masks_yeo_aal3, 'rAG.nii')           , 'rAG_yeo'             , 'Right_AG_Yeo'
%      fullfile(path_masks_yeo_aal3, 'lIFG.nii')          , 'lIFG_yeo'            , 'Left_IFG_Yeo'
%      fullfile(path_masks_yeo_aal3, 'rIFG.nii')          , 'rIFG_yeo'            , 'Rigth_IFG_Yeo'
%      fullfile(path_masks_yeo_aal3, 'lMFG.nii')          , 'lMFG_yeo'            , 'Left_MFG_Yeo'
%      fullfile(path_masks_yeo_aal3, 'rMFG.nii')          , 'rMFG_yeo'            , 'Right_IFG_Yeo'
%      fullfile(path_masks_yeo_aal3, 'lMPFC.nii')         , 'lMPFC_yeo'           , 'Left_MPFC_Yeo'
%      fullfile(path_masks_yeo_aal3, 'rMPFC.nii')         , 'rMPFC_yeo'           , 'Right_MPFC_Yeo'
%      fullfile(path_masks_yeo_aal3, 'lMTC.nii')          , 'lMTC_yeo'            , 'Left_MTC_Yeo'
%      fullfile(path_masks_yeo_aal3, 'rMTC.nii')          , 'rMTC_yeo'            , 'Right_MTC_Yeo'
%      fullfile(path_masks_yeo_aal3, 'lPMC.nii')          , 'lPMC_yeo'            , 'Left_PMC_Yeo'
%      fullfile(path_masks_yeo_aal3, 'rPMC.nii')          , 'rPMC_yeo'            , 'Right_PMC_Yeo'
   } ;

par.roi_type.sphere_global = {
%     % [x y z]mm      radius(mm)   abbrev    fullname
        [+2,-42,+32],  15,          'Cingulate', 'Cingulate_Cluster1_Moring2022'
        [-26,+46,-2],  10,          'OFC', 'Orbitofrontal_Cortex_Cluster_Haupt2019'
%     [ -14,+4,-3 ],   3,          's3lGPe', 's3 left GPe'
%     [ -15,-2,-4 ],   2,          's2lGPi', 's2 left GPi'
    };

% perform the extraction

TS = job_extract_timeseries(par);

%%
% % define some networks : not mandatory


% par.network.DMN = {'lSFGmedial', 'rSFGmedial','lPFCventmed', 'rPFCventmed','lREC', 'rREC','lPCC', 'rPCC', 'lTPOsup', 'rTPOsup', 'lTPOmid', 'rTPOmid', ...
%     'lMTG', 'rMTG','lSMG',  'rSMG','lANG', 'rANG', 'lSTG', 'rSTG', 'lIPG', 'rIPG'   };

% par.network.DMN = {'lAG_yeo', 'rAG_yeo','lIFG_yeo', 'rIFG_yeo','lMFG_yeo', 'rMFG_yeo','lMPFC_yeo', 'rMPFC_yeo', 'lMTC_yeo', 'rMTC_yeo', 'lPMC_yeo', 'rPMC_yeo' };


% par.network.salience = { 'lINS','rINS', 'lACCsub',  'rACCsub','lOLF', 'rOLF', 'lVTA', 'rVTA'   } ;

% par.network.auditory = { 'lAudio',  'rAudio', 'ltMDl',  'rtMDl' } ;

% par.network.limbic = {'lAMYG', 'rAMYG', 'lHIP', 'rHIP', 'lPHG', 'rPHG', 'lNacc', 'rNacc' } ;  % Parahippocampe est inclue ici
% par.network.limbic = {'lAMYG', 'rAMYG', 'lHIP', 'rHIP', 'lNacc', 'rNacc' } ;

% par.network.attention = {'r_parietal_nback' , 'l_parietal_nback', 'r_MedialFG_nback', 'l_MedialFG_nback' } ;



% perform connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

%% plot
guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% perfrom seed-to-voxel correlation (seed == ROI)
par.jobname = 'seed2brain_analysis_AUDICOG';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 

%% Sauvegarder les matrices sur owncloud pour pouvoir les récupérer sur mac

output_dir = '/home/lise.hobeika/ownCloud/Postdoc_acouphènes/Manips/Imagerie/Data/Matrices_correlations'; 

for isubj = 1:length(e)

    ifile =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/connectivity__DMN_salience_auditory_limbic_attention__timeseries__aal3__lrlrlrlrlrlrlrlrlr.mat') ;
    load(ifile)   ;

    ifile_out = ['S' , num2str(isubj) ,'_connectivity_RS.mat' ]; 

    
copyfile( ifile,  fullfile( output_dir, ifile_out ) ) ; 

%     ifile_out = ['S' , num2str(isubj) ,'_connectivity_timeseries_aal3.xlsx' ]; 
%     writematrix(  connectivity_matrix, fullfile( output_dir, ifile_out )  ) ; % , 'conn_network',  )  ;
% %   
%     ifile_out = ['S' , num2str(isubj) ,'_connectivity_timeseries_Salience.xlsx' ]; 
%     writematrix(  network(2).mx  , fullfile( output_dir, ifile_out )  ) ; % , 'conn_network',  )  ;
%   
end



%% Save le tableau de correspondance de l'atlas

% %     save(fullfile( output_dir, 'info_networks.mat' ) , 'connectivity_matrix' , 'conn_network',  )  ;
%     writetable(  ts_table, fullfile( output_dir, 'Correspondances_aal3.xlsx' )  ) ; % , 'conn_network',  )  ;

        writetable(  network(2).table , fullfile( output_dir, 'Correspondances_Salience.xlsx' )  ) ; % , 'conn_network',  )  ;

