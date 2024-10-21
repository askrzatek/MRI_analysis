%% connectivity
clear
clc
main_dir = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG';
ROI_dir = fullfile(main_dir,'Networks_Masks');


cd (main_dir)
load('e_nonchir.mat'); 


%% Step 1: faire tourner les calculs de connnectivité

% define input volumes and confounds
clear par
par.run             = 0;
par.sge             = 1;
par.mem             = '16G';
par.jobname         = 'timeseries_extract_AUDICOG_salience_tinnitus';
par.display         = 0;
par.redo            = 0;
par.volume          =  e.getSerie('run_RS').getVolume('s5wts_OC') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound        =  e.getSerie('run_RS').getRP('multiple_regressors');
par.mask_threshold  = 0.001;

%% define ROI, using several methods

% par.roi_type.atlas_cat12 = 'aal3';

path_masks_audio          = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Results/Audio/2023_01_27 - Audio' ;
path_masks_Tinnitus       = get_subdir_regex(ROI_dir, 'Tinnitus_Meta');
path_masks_Alert          = get_subdir_regex(ROI_dir, 'Alerting_effect');
path_masks_Salience_CAREN = '/network/lustre/iss02/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01';

%% get labels for CAREN ROIs from table according to RSN we are interested in : in this case Salience = RSN01

tab          = readtable(fullfile(ROI_dir,'aal_id_labels_abbrevs_CAREN.csv'));
salience_net.labels  = tab.label_AAL(tab.salience(~isnan(tab.salience)));
salience_net.abbrevs = tab.abbrev_AAL(tab.salience(~isnan(tab.salience)));
atlas_rois_list = [];
% salience_net.rois    = get_subdir_regex_files(path_masks_Salience_CAREN,'.*');

for iroi = 1 : length(salience_net.labels)
                                %     % path                                                                                                abbrev                      description
    atlas_rois_list = [atlas_rois_list; char(get_subdir_regex_files(path_masks_Salience_CAREN,sprintf('.*Reg%d.nii$',tab.salience(iroi)))), salience_net.abbrevs(iroi), salience_net.labels(iroi)];
end

%% we could do the same for other networks : just to have the timeseries extracted per network if we want to compare


%% define all ROIs independently from the atlas

par.roi_type.mask_global = vertcat(atlas_rois_list,{
%     % path                                                           abbrev          description
       char(get_subdir_regex_files(path_masks_Tinnitus, 'ParaHipp')), 'ParaHipp',  'Bilateral_ParaHippocampus'
       fullfile( path_masks_audio,  'Left_Auditory_activation_fwe05_0v_50subj.nii'), 'lAudio'   ,  'Left_Loca_Audio_activation'   
       fullfile( path_masks_audio,  'Right_Auditory_activation_fwe05_0v_50subj.nii'), 'rAudio'   ,  'Right_Loca_Audio_activation'   

%        char(get_subdir_regex_files(path_masks_Alert, 'Frontal_Orb')),          'Orb_PFC',        'Combined_Orbital_Prefrontal_Cortex'

% NOT ENOUGH SIGNAL IN 3T BOLD ACQUISITIONS SO AFTER RESLICING THE MASK, THE MASK BECOMES EMPTY FOR MIDBRAIN STRUCTURES LIKE LC
%        char(get_subdir_regex_files(path_masks_Alert, 'LC_l')),          'LC_l',        'Left_Locus_Coereleus'
%        char(get_subdir_regex_files(path_masks_Alert, 'LC_r')),          'LC_r',        'Right_Locus_Coereleus'
%        char(get_subdir_regex_files(path_masks_Alert, 'LC.nii')),          'LC',        'Bilateral_Locus_Coereleus'

   }) ;

%% define ROIs accordingly to the coordinates from the literature

par.roi_type.sphere_global = {
%     % [x y z]mm      radius(mm)   abbrev    fullname
        [+2,-42,+32],  15,          'Cingulate', 'Cingulate_Cluster1_Moring2022'
        [-26,+46,-2],  10,          'OFC', 'Orbitofrontal_Cortex_Cluster_Haupt2019'
%     [ -14,+4,-3 ],   3,          's3lGPe', 's3 left GPe'
%     [ -15,-2,-4 ],   2,          's2lGPi', 's2 left GPi'
    };

%% Perform the timeseries' extraction

TS = job_extract_timeseries(par);

%% Step 2: create correlation matrix

%% Define some networks : not mandatory


% par.network.DMN = {'lSFGmedial', 'rSFGmedial','lPFCventmed', 'rPFCventmed','lREC', 'rREC','lPCC', 'rPCC', 'lTPOsup', 'rTPOsup', 'lTPOmid', 'rTPOmid', ...
%     'lMTG', 'rMTG','lSMG',  'rSMG','lANG', 'rANG', 'lSTG', 'rSTG', 'lIPG', 'rIPG'   };

% par.network.DMN = {'lAG_yeo', 'rAG_yeo','lIFG_yeo', 'rIFG_yeo','lMFG_yeo', 'rMFG_yeo','lMPFC_yeo', 'rMPFC_yeo', 'lMTC_yeo', 'rMTC_yeo', 'lPMC_yeo', 'rPMC_yeo' };


% par.network.salience = { 'lINS','rINS', 'lACCsub',  'rACCsub','lOLF', 'rOLF', 'lVTA', 'rVTA'   } ;

% par.network.auditory = { 'lAudio',  'rAudio', 'ltMDl',  'rtMDl' } ;

% par.network.limbic = {'lAMYG', 'rAMYG', 'lHIP', 'rHIP', 'lPHG', 'rPHG', 'lNacc', 'rNacc' } ;  % Parahippocampe est inclue ici
% par.network.limbic = {'lAMYG', 'rAMYG', 'lHIP', 'rHIP', 'lNacc', 'rNacc' } ;

% par.network.attention = {'r_parietal_nback' , 'l_parietal_nback', 'r_MedialFG_nback', 'l_MedialFG_nback' } ;



%% Create connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

%% plot
guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% Step 3: Perfrom seed-to-voxel correlation (seed == ROI)
par.jobname = 'seed2brain_analysis_AUDICOG';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 

%% Step 4: Sauvegarder les matrices sur owncloud pour pouvoir les récupérer sur mac

output_dir = '/home/anna.skrzatek/AUDICOG_backup/Correl_Matrix'; 

for isubj = 1:length(e)

    src =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__ppfffffffffrrSSfiiMMpppssppppttttAAPlr__CO.mat') ;
    dst = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries_Salience_Tinnitus_Alert_Audio.mat') ;
    src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__ppfffffffffrrSSfiiMMpppssppppttttAAPlr__CO.mat') ;
    dst2 = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn_timeseries_Salience_Tinnitus_Alert_Audio.mat') ;
    
    if exist(src, 'file')
        movefile(src, dst);
        movefile(src2, dst2);
    else
        fprintf('expected file "%s" not found\n', src);
        fprintf('expected file "%s" not found\n', src2);
    end
    
    %% timeseries matrix
    ifile = dst;
    rs_cmat = load(ifile)   ;

    ifile_out = ['AUDICOG_Sujet' , num2str(isubj) ,'_connectivity_RS_Salience_CAREN_Tinnitus_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2;
    rs_cmat = load(ifile)   ;

    ifile_out = ['AUDICOG_Sujet' , num2str(isubj) ,'_connectivity_RS_Salience_CAREN_Tinnitus_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_Salience_CAREN_Tinnitus_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'));

    %writematrix(rs_cmat,addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix'))
    %%% only available on Matlab 2019b
%    copyfile( ifile,  fullfile( output_dir, ifile_out )) ;

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

