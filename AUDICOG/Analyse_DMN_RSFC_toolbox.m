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
par.jobname = 'timeseries_extract_AUDICOG_DMN_CAREN_Robust_hubs';
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
path_masks_CAREN    = get_subdir_regex(ROI_dir, 'CAREN');
path_masks_DMN    = get_subdir_regex(path_masks_CAREN, 'RSN04');

par.roi_type.mask_global = {
%     % path                                                                        abbrev          description
       char(get_subdir_regex_files(path_masks_DMN, 'Reg43')),                       'ParaHipp_l',   'Left_ParaHippocampus'
       char(get_subdir_regex_files(path_masks_DMN, 'Reg44')),                       'ParaHipp_r',   'Right_ParaHippocampus'

       char(get_subdir_regex_files(path_masks_CAREN, 'ACC_L_all_DMN_CAREN')),       'ACC_l',        'Left_Anterior_Cingulate'
       char(get_subdir_regex_files(path_masks_CAREN, 'ACC_R_all_DMN_CAREN')),       'ACC_r',        'Right_Anterior_Cingulate'

       char(get_subdir_regex_files(path_masks_DMN, 'Reg39')),                       'PCC_l',   'Left_ParaHippocampus'
       char(get_subdir_regex_files(path_masks_DMN, 'Reg40')),                       'PCC_r',   'Right_ParaHippocampus'

       char(get_subdir_regex_files(path_masks_DMN, 'Reg65')),                       'Pariet_Inf_l', 'Left_Inferior_Parietal'
       char(get_subdir_regex_files(path_masks_DMN, 'Reg66')),                       'Pariet_Inf_r', 'Right_Inferior_Parietal'

       char(get_subdir_regex_files(path_masks_DMN, 'Reg89')),                       'Temp_Mid_l',   'Left_Middle_Temporal'
       char(get_subdir_regex_files(path_masks_DMN, 'Reg90')),                       'Temp_Mid_r',   'Right_Middle_Temporal'
       
       char(get_subdir_regex_files(path_masks_CAREN, 'Front_Inf_L_all_DMN_CAREN')), 'Front_Inf_l',  'Left_Inferior_Frontal'
       char(get_subdir_regex_files(path_masks_CAREN, 'Front_Inf_R_all_DMN_CAREN')), 'Front_Inf_r',  'Right_Inferior_Frontal'

       char(get_subdir_regex_files(path_masks_DMN, 'Reg05')),                       'Front_Mid_l',  'Left_Middle_Frontal'
       char(get_subdir_regex_files(path_masks_DMN, 'Reg06')),                       'Front_Mid_r',  'Right_Middle_Frontal'
       
       char(get_subdir_regex_files(path_masks_DMN, 'Reg71')),                       'Precuneus_l',  'Left_Precuneus'
       char(get_subdir_regex_files(path_masks_DMN, 'Reg72')),                       'Precuneus_r',  'Right_Precuneus'

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

%par.roi_type.sphere_global = {
%     % [x y z]mm      radius(mm)   abbrev    fullname
%        [+2,-42,+32],  15,          'Cingulate', 'Cingulate_Cluster1_Moring2022'

%         [-26,+46,-2],  10,          'OFC', 'Orbitofrontal_Cortex_Cluster_Haupt2019'

%     [ -14,+4,-3 ],   3,          's3lGPe', 's3 left GPe'
%     [ -15,-2,-4 ],   2,          's2lGPi', 's2 left GPi'
%    };

% perform the extraction
par.outname = 'AUDICOG_PH_ACC_PCC_ParInf_TempMid_FrontInf_FrontMid_Precun_L_R';
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

% rename and save files in a common directory
output_dir = '/home/anna.skrzatek/AUDICOG_backup/Correl_Matrix'; 

for isubj = 1:length(e)
%% part used to rename the files in the subject directories if one forgot to used the par.output_name in the rsfc - therefore for further copying we only need dst files
%    src =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst_DMN = fullfile( e(isubj).getSerie('run_RS').path, 'rsfc/timeseries__AUDICOG_PH_ACC_PCC_ParInf_TempMid_FrontInf_FrontMid_Precun_L_R.mat') ;
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst2_DMN = fullfile( e(isubj).getSerie('run_RS').path, 'rsfc/static_conn__timeseries__AUDICOG_PH_ACC_PCC_ParInf_TempMid_FrontInf_FrontMid_Precun_L_R.mat') ;
        
%     if exist(src, 'file')
%         movefile(src, dst);
%         movefile(src2, dst2);
%     else
          if ~exist(dst_DMN,'file')
              fprintf('expected file "%s" not found\n', dst_DMN);
          end
          
          if ~exist(dst2_DMN,'file')
              fprintf('expected file "%s" not found\n', dst2_DMN);
          end
%         fprintf('expected file "%s" not found\n', src);
%         fprintf('expected file "%s" not found\n', src2);
%    end

    %% timeseries matrix
    ifile = dst_DMN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_AUDICOG_PH_ACC_PCC_ParietInf_TempMid_FrontInf_FrontMid_Precun_L_R_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_DMN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_AUDICOG_PH_ACC_PCC_ParietInf_TempMid_FrontInf_FrontMid_Precun_L_R_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_DMN_united_CAREN2_labels.csv';
    %ifile_out = 'Connectivity_RS_SM_Visual_Net_CAREN_Tinnitus_Alert_AudioACT_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'));

end

%% plot
guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% perfrom seed-to-voxel correlation (seed == ROI)
par.jobname = 'seed2brain_analysis_AUDICOG';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 


%% Save le tableau de correspondance de l'atlas

% %     save(fullfile( output_dir, 'info_networks.mat' ) , 'connectivity_matrix' , 'conn_network',  )  ;
%     writetable(  ts_table, fullfile( output_dir, 'Correspondances_aal3.xlsx' )  ) ; % , 'conn_network',  )  ;

%        writetable(  network(1).table , fullfile( output_dir, 'Correspondances_DMN.xlsx' )  ) ; % , 'conn_network',  )  ;

