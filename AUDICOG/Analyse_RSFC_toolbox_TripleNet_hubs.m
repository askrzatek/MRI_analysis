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
par.jobname = 'timeseries_extract_AUDICOG_left_right_sep';
par.display = 0;
par.redo = 0;
par.volume   =  e.getSerie('run_RS').getVolume('s5wts_OC') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound =  e.getSerie('run_RS').getRP('multiple_regressors');
par.mask_threshold = 0.001;

% define ROI, using several methods
%par.roi_type.atlas_cat12 = 'aal3';

path_masks_DMN    = get_subdir_regex(path_masks_CAREN, 'RSN04');
path_masks_SN    = get_subdir_regex(path_masks_CAREN, 'RSN01');
path_masks_CEN    = get_subdir_regex(path_masks_CAREN, 'RSN02');


par.roi_type.mask_global = {
%     % path                                                                        abbrev          description
% DMN
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
% SN
       char(get_subdir_regex_files(path_masks_SN, 'Reg19')),                        'Front_Sup_Med_l', 'Left_Superior_Frontal_Medial'
       
       char(get_subdir_regex_files(path_masks_SN, 'Reg33')),                        'Insula_l',  'Left_Insula'
       char(get_subdir_regex_files(path_masks_SN, 'Reg34')),                        'Insula_r',  'Right_Insula'

       char(get_subdir_regex_files(path_masks_SN, 'Reg67')),                        'SupraMarginal_l',     'Left_Supramarginal'
       char(get_subdir_regex_files(path_masks_SN, 'Reg68')),                        'Supramarginal_r',     'Right_Supramarginal'

       char(get_subdir_regex_files(path_masks_SN, 'Reg155')),                       'ACC_Sup_l',     'Left_Superior_Anterior_Cingulate'
       char(get_subdir_regex_files(path_masks_SN, 'Reg156')),                       'ACC_Sup_r',     'Right_Superior_Anterior_Cingulate'
% CEN
       char(get_subdir_regex_files(path_masks_CEN, 'Reg03')),                        'Front_Sup_l', 'Left_Superior_Frontal'
       char(get_subdir_regex_files(path_masks_CEN, 'Reg04')),                        'Front_Sup_r', 'Right_Superior_Frontal'
       
       char(get_subdir_regex_files(path_masks_CEN, 'Reg05')),                        'Front_Mid_l', 'Left_Middle_Frontal'
       char(get_subdir_regex_files(path_masks_CEN, 'Reg06')),                        'Front_Mid_r', 'Right_Middle_Frontal'
       
       char(get_subdir_regex_files(path_masks_CEN, 'Reg15')),                         'SMA_l',     'Left_Supplementary_Motor_Area'
       char(get_subdir_regex_files(path_masks_CEN, 'Reg16')),                         'SMA_r',     'Right_Supplementary_Motor_Area'

       char(get_subdir_regex_files(path_masks_CEN, 'Reg19')),                        'Front_Sup_Med_l', 'Left_Superior_Frontal_Medial'
       char(get_subdir_regex_files(path_masks_CEN, 'Reg20')),                        'Front_Sup_Med_r', 'Right_Superior_Frontal_Medial'
       
       char(get_subdir_regex_files(path_masks_CEN, 'Reg63')),                         'Pariet_Sup_l',     'Left_Superior_Parietal'
       char(get_subdir_regex_files(path_masks_CEN, 'Reg64')),                         'Pariet_Sup_r',     'Right_Superior_Parietal'
       
       char(get_subdir_regex_files(path_masks_CEN, 'Reg65')),                        'Pariet_Inf_l',     'Left_Inferior_Parietal'
       char(get_subdir_regex_files(path_masks_CEN, 'Reg66')),                        'Pariet_Inf_r',     'Right_Inferior_Parietal'

%        char(get_subdir_regex_files(ROI_dir, 'Heschl_R_L')),           'A1',        'Bilateral_Heschl_Gyrus'
%        char(get_subdir_regex_files(ROI_dir, 'A1_L')),                      'A1_l',        'Left_Heschl_Gyrus'
%        char(get_subdir_regex_files(ROI_dir, 'A1_R')),                      'A1_r',        'Right_Heschl_Gyrus'       
   } ;

%par.roi_type.sphere_global = {
%     % [x y z]mm      radius(mm)   abbrev    fullname
%        [+2,-42,+32],  15,          'Cingulate', 'Cingulate_Cluster1_Moring2022'

%         [-26,+46,-2],  10,          'OFC', 'Orbitofrontal_Cortex_Cluster_Haupt2019'

%     [ -14,+4,-3 ],   3,          's3lGPe', 's3 left GPe'
%     [ -15,-2,-4 ],   2,          's2lGPi', 's2 left GPi'
%    };

% perform the extraction
par.outname = 'AUDICOG_DMN_SN_CEN_Hubs_CAREN';
TS = job_extract_timeseries(par);

%%
% % define some networks : not mandatory

par.network.DMN = {'ParaHipp_l','ParaHipp_r','ACC_l','ACC_r','PCC_l','PCC_r','Pariet_Inf_l','Pariet_Inf_r','Temp_Mid_l','Temp_Mid_r','Front_Inf_l','Front_Inf_r','Front_Mid_l','Front_Mid_r','Precuneus_l','Precuneus_r'};
par.network.SN  = {'Front_Sup_Med_l','Insula_l','Insula_r','SupraMarginal_l','Supramarginal_r','ACC_Sup_l','ACC_Sup_r'} ;
par.network.CEN = {'Front_Sup_r','Front_Sup_l','Front_Mid_l','Front_Mid_r','SMA_l','SMA_r','Front_Sup_Med_l','Front_Sup_Med_r','Pariet_Sup_l','Pariet_Sup_r','Pariet_Inf_l','Pariet_Inf_r'}

% perform connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

% rename and save files in a common directory
output_dir = '/home/anna.skrzatek/AUDICOG_backup/Correl_Matrix'; 

for isubj = 1:length(e)
%% part used to rename the files in the subject directories if one forgot to used the par.output_name in the rsfc - therefore for further copying we only need dst files
%    src =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst_3Net = fullfile( e(isubj).getSerie('run_RS').path, 'rsfc/timeseries__AUDICOG_DMN_SN_CEN_Hubs_CAREN.mat') ;
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst2_3Net = fullfile( e(isubj).getSerie('run_RS').path, 'rsfc/static_conn__timeseries__AUDICOG_DMN_SN_CEN_Hubs_CAREN.mat') ;
        
%     if exist(src, 'file')
%         movefile(src, dst);
%         movefile(src2, dst2);
%     else
          if ~exist(dst_3Net,'file')
              fprintf('expected file "%s" not found\n', dst_3Net);
          end
          
          if ~exist(dst2_3Net,'file')
              fprintf('expected file "%s" not found\n', dst2_3Net);
          end
%         fprintf('expected file "%s" not found\n', src);
%         fprintf('expected file "%s" not found\n', src2);
%    end

    %% timeseries matrix
    ifile = dst_3Net;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_AUDICOG_DMN_SN_CEN_Hubs_CAREN_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_3Net;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_AUDICOG_DMN_SN_CEN_Hubs_CAREN_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_DMN_SN_CEN_Hubs_CAREN_labels.csv';
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

