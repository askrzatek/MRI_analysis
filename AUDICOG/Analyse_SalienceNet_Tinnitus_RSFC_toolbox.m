%% connectivity
clear
clc
addpath('/home/anna.skrzatek/MRI_analysis/')
addpath('/network/iss/cenir/analyse/irm/users/salim.ouarab/')

main_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG';
ROI_dir = fullfile(main_dir,'Networks_Masks');


cd (fullfile(main_dir,'DATA/Non_chirurgicaux'))
%load('e_nonchir.mat'); 
e = exam(fullfile(main_dir,'DATA/Non_chirurgicaux'), '.*AUDICOG_Suj'); % all subjects with multi-echo


%% Get files paths #matvol

% % Anat
e.addSerie('T1w$', 'anat_T1', 1 );
e.getSerie('anat').addVolume('^v_.*nii','v',1);

% Func
run_name = 'RS';
e.addSerie([run_name           '$'], ['run_' run_name], 1);
e.addSerie('RS','tedana','tedana',1)

e.getSerie('run').addVolume('^v_.*nii$',   'v', 3);
e.getSerie('tedana').addVolume('^s5wts_OC','s5wts') %% autoadd job
e.explore

rp = fullfile(e.getSerie('tedana').getPath(),'multiple_regressors.txt');

%% Step 1: faire tourner les calculs de connnectivité

% define input volumes and confounds
clear par
par.run             = 1;
par.sge             = 0;
par.mem             = '16G';
par.jobname         = 'timeseries_extract_AUDICOG_salience_tinnitus';
par.jobname         = 'timeseries_extract_AUDICOG_cEx_tinnitus';
par.jobname         = 'timeseries_extract_AUDICOG_DMN_tinnitus';
par.jobname         = 'timeseries_extract_AUDICOG_SMN_VIS';
par.jobname         = 'timeseries_extract_AUDICOG_CAREN_all_v2';
par.display         = 0;
par.redo            = 0;
par.volume          = e.getSerie('tedana').getVolume('s5wts');
par.confound        = rp;
% par.volume          =  e.getSerie('run_RS').getVolume('s5wts_OC') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
% par.confound        =  e.getSerie('run_RS').getRP('multiple_regressors');
par.mask_threshold  = 0.001;

%% define ROI, using several methods

% par.roi_type.atlas_cat12 = 'aal3';

path_masks_audio          = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/Audio/2023_01_27 - Audio' ;
path_masks_Tinnitus       = get_subdir_regex(ROI_dir, 'Tinnitus_Meta');
path_masks_Alert          = get_subdir_regex(ROI_dir, 'Alerting_effect');
path_masks_Salience_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01';
path_masks_CEN_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN02';
path_masks_DMN_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN04';
path_masks_SMN_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN03';
path_masks_VIS_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN05';

%% get labels for CAREN ROIs from table according to RSN we are interested in : in this case Salience = RSN01

tab          = readtable(fullfile(ROI_dir,'aal_id_labels_abbrevs_CAREN.csv'));

salience_net.labels  = strcat(tab.label_AAL(tab.salience(~isnan(tab.salience))),'_SN');
salience_net.abbrevs = strcat(tab.abbrev_AAL(tab.salience(~isnan(tab.salience))),'_SN');
salience_net.idx = tab.id_AAL(tab.salience(~isnan(tab.salience)));
cEx_net.labels  = strcat(tab.label_AAL(tab.cen(~isnan(tab.cen))),'_CEN');
cEx_net.abbrevs = strcat(tab.abbrev_AAL(tab.cen(~isnan(tab.cen))),'_CEN');
cEx_net.idx = tab.id_AAL(tab.cen(~isnan(tab.cen)));
DM_net.labels  = strcat(tab.label_AAL(tab.dmn(~isnan(tab.dmn))),'_DMN');
DM_net.abbrevs = strcat(tab.abbrev_AAL(tab.dmn(~isnan(tab.dmn))),'_DMN');
DM_net.idx = tab.id_AAL(tab.dmn(~isnan(tab.dmn)));
SM_net.labels  = strcat(tab.label_AAL(tab.smn(~isnan(tab.smn))),'_SMN');
SM_net.abbrevs = strcat(tab.abbrev_AAL(tab.smn(~isnan(tab.smn))),'_SMN');
SM_net.idx = tab.id_AAL(tab.smn(~isnan(tab.smn)));
Vis_net.labels  = strcat(tab.label_AAL(tab.vis(~isnan(tab.vis))),'_VIS');
Vis_net.abbrevs = strcat(tab.abbrev_AAL(tab.vis(~isnan(tab.vis))),'_VIS');
Vis_net.idx = tab.id_AAL(tab.vis(~isnan(tab.vis)));


salience_net.skip = [];
cEx_net.skip = [57];
DM_net.skip = [51,93];
SM_net.skip = [6];
Vis_net.skip = [93];

atlas_rois_list_SN = [];
atlas_rois_list_CEN = [];
atlas_rois_list_DMN = [];
atlas_rois_list_SMN = [];
atlas_rois_list_Vis = [];

% salience_net.rois    = get_subdir_regex_files(path_masks_Salience_CAREN,'.*');

for iroi = 1 : length(salience_net.labels)
% nroi = dir(path_masks_Salience_CAREN);
% for iroi = 1: length(nroi)-2
    if ~ismember(salience_net.idx(iroi), salience_net.skip)
        if salience_net.idx(iroi) < 10
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_SN = [atlas_rois_list_SN; char(get_subdir_regex_files(path_masks_Salience_CAREN,sprintf('.*Reg0%d.nii$',salience_net.idx(iroi)))), salience_net.abbrevs(iroi), salience_net.labels(iroi)];
        else
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_SN = [atlas_rois_list_SN; char(get_subdir_regex_files(path_masks_Salience_CAREN,sprintf('.*Reg%d.nii$',salience_net.idx(iroi)))), salience_net.abbrevs(iroi), salience_net.labels(iroi)];
        end
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', salience_net.idx(iroi))
    end
end

for iroi = 1 : length(cEx_net.labels)
%nroi = dir(path_masks_CEN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(cEx_net.idx(iroi), cEx_net.skip)
        if cEx_net.idx(iroi) < 10
                                        %     % path                                                                                                abbrev                      description
            atlas_rois_list_CEN = [atlas_rois_list_CEN; char(get_subdir_regex_files(path_masks_CEN_CAREN,sprintf('.*Reg0%d.nii$',cEx_net.idx(iroi)))), cEx_net.abbrevs(iroi), cEx_net.labels(iroi)];

        else
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_CEN = [atlas_rois_list_CEN; char(get_subdir_regex_files(path_masks_CEN_CAREN,sprintf('.*Reg%d.nii$',cEx_net.idx(iroi)))), cEx_net.abbrevs(iroi), cEx_net.labels(iroi)];
    
        end
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', cEx_net.idx(iroi))
    end
end

for iroi = 1 : length(DM_net.labels)
%nroi = dir(path_masks_DMN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(DM_net.idx(iroi), DM_net.skip)
        if DM_net.idx(iroi) < 10
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_DMN = [atlas_rois_list_DMN; char(get_subdir_regex_files(path_masks_DMN_CAREN,sprintf('.*Reg0%d.nii$',DM_net.idx(iroi)))), DM_net.abbrevs(iroi), DM_net.labels(iroi)];

        else
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_DMN = [atlas_rois_list_DMN; char(get_subdir_regex_files(path_masks_DMN_CAREN,sprintf('.*Reg%d.nii$',DM_net.idx(iroi)))), DM_net.abbrevs(iroi), DM_net.labels(iroi)];
        end
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', DM_net.idx(iroi))
    end
end

for iroi = 1 : length(SM_net.labels)

    if ~ismember(SM_net.idx(iroi), SM_net.skip)
        if SM_net.idx(iroi) < 10
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_SMN = [atlas_rois_list_SMN; char(get_subdir_regex_files(path_masks_SMN_CAREN,sprintf('.*Reg0%d.nii$',SM_net.idx(iroi)))), SM_net.abbrevs(iroi), SM_net.labels(iroi)];

        else
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_SMN = [atlas_rois_list_SMN; char(get_subdir_regex_files(path_masks_SMN_CAREN,sprintf('.*Reg%d.nii$',SM_net.idx(iroi)))), SM_net.abbrevs(iroi), SM_net.labels(iroi)];
        end
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', SM_net.idx(iroi))
    end
end

for iroi = 1 : length(Vis_net.labels)
%nroi = dir(path_masks_DMN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(Vis_net.idx(iroi), Vis_net.skip)
        if Vis_net.idx(iroi) < 10
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_Vis = [atlas_rois_list_Vis; char(get_subdir_regex_files(path_masks_Vis_CAREN,sprintf('.*Reg0%d.nii$',Vis_net.idx(iroi)))), Vis_net.abbrevs(iroi), Vis_net.labels(iroi)];

        else
                                %     % path                                                                                                abbrev                      description
            atlas_rois_list_Vis = [atlas_rois_list_Vis; char(get_subdir_regex_files(path_masks_VIS_CAREN,sprintf('.*Reg%d.nii$',Vis_net.idx(iroi)))), Vis_net.abbrevs(iroi), Vis_net.labels(iroi)];
        end
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', Vis_net.idx(iroi))
    end
end

%% we could do the same for other networks : just to have the timeseries extracted per network if we want to compare


%% define all ROIs independently from the atlas

par.roi_type.mask_global = vertcat(atlas_rois_list_SN, atlas_rois_list_CEN, atlas_rois_list_DMN, atlas_rois_list_SMN, atlas_rois_list_Vis, {
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
par.outname = 'AUDICOG_CAREN_SN_CEN_DMN_AudioACT';
par.outname = 'Salience_Tinnitus_Alert_Audio';
par.outname = 'Executive_Tinnitus_Alert_Audio';
par.outname = 'DMN_Tinnitus_Alert_Audio';
par.outname = 'SN_CEN_DMN_SMN_VIS_CAREN';
par.outname = 'AUDICOG_CAREN_3Net';
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

par.network.SN = atlas_rois_list_SN(:,2);
par.network.CEN = atlas_rois_list_CEN(:,2);
par.network.DMN = atlas_rois_list_DMN(:,2);
par.network.SMN = atlas_rois_list_SMN(:,2);
par.network.Vis = atlas_rois_list_Vis(:,2);

%% Create connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

%% plot
% fichier de correspondance numero IRM - comportement - groupe - age
cd(main_dir)
d = readtable( [ './DATA/' , 'Correspondance_Numero_Comportement_IRM.csv' ])  ;

% séparation des paths en fonction des groupes expérimentaux
groups.name = {'Control','Tinnitus'};
groups.TS = cell((length(e)/2),length(groups.name));
groups.e  = cell((length(e)/2),length(groups.name));
for igroup = 1:length(groups.name)
    j = 0 ;
    for iSubj = 1:length(e)
        ifile = e(iSubj).name;
        id = str2double(ifile(25:end)) ;
        subj_group = d.Groupe(find(d.Num_IRM == id));

        if subj_group == igroup
            j = j + 1 ;
            groups.TS{j,igroup} = TS(iSubj);
            groups.e{j,igroup}  = e(iSubj);
        end
    end
end    
con_FC = vertcat(groups.TS{:,1});
tin_FC = vertcat(groups.TS{:,2});

control_guidata = plot_resting_state_connectivity_matrix(con_FC);
tinnitus_guidata = plot_resting_state_connectivity_matrix(tin_FC);

guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% Step 3: Perfrom seed-to-voxel correlation (seed == ROI)
par.jobname = 'seed2brain_analysis_AUDICOG';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 

%% Step 4: Sauvegarder les matrices sur owncloud pour pouvoir les récupérer sur mac

output_dir = '/home/anna.skrzatek/AUDICOG_backup/Correl_Matrix'; 

for isubj = 1:length(e)
%% part used to rename the files in the subject directories if one forgot to used the par.output_name in the rsfc - therefore for further copying we only need dst files
%    src =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst_allnet = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__SN_CEN_DMN_SMN_VIS_CAREN.mat');
    
    dst_tinnet = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat');
    dst_controlnet = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_SM_Visual_Net_Tinnitus_Alert_Audio.mat');
    dst_SN = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__Salience_Tinnitus_Alert_Audio.mat') ;
    dst_CEN = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__Executive_Tinnitus_Alert_Audio.mat') ;
    dst_DMN = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__DMN_Tinnitus_Alert_Audio.mat') ;

    dst_3Net = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_3Net.mat') ;

%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst2_allnet = fullfile(e(isubj).getSerie('run_RS').path,'tedana009a1_vt/rsfc/SN_CEN_DMN_SMN_Vis__static_conn__timeseries__SN_CEN_DMN_SMN_VIS_CAREN.mat');
    
    dst2_tinnet = fullfile(e(isubj).getSerie('run_RS').path,'tedana009a1_vt/rsfc/SN_CEN_DMN__static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat');
    dst2_controlnet = fullfile(e(isubj).getSerie('run_RS').path,'tedana009a1_vt/rsfc/SMN_Vis__static_conn__timeseries__AUDICOG_CAREN_SM_Visual_Net_Tinnitus_Alert_Audio.mat');
    dst2_SN = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/SN__static_conn__timeseries__Salience_Tinnitus_Alert_Audio.mat') ;
    dst2_CEN = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/CEN__static_conn__timeseries__Executive_Tinnitus_Alert_Audio.mat') ;
    dst2_DMN = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/DMN__static_conn__timeseries__DMN_Tinnitus_Alert_Audio.mat') ;
        
    dst2_3Net = fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/SN_CEN_DMN__static_conn__timeseries__AUDICOG_CAREN_3Net.mat') ;
%     if exist(src, 'file')
%         movefile(src, dst);
%         movefile(src2, dst2);
%     else
          if ~exist(dst_allnet,'file')
              fprintf('expected file "%s" not found\n', dst_allnet);
          end
          
          if ~exist(dst2_allnet,'file')
              fprintf('expected file "%s" not found\n', dst2_allnet);
          end
          
          if ~exist(dst_SN,'file')
              fprintf('expected file "%s" not found\n', dst_SN);
          end
          
          if ~exist(dst2_SN,'file')
              fprintf('expected file "%s" not found\n', dst2_SN);
          end
          
          if ~exist(dst_CEN,'file')
              fprintf('expected file "%s" not found\n', dst_CEN);
          end
          
          if ~exist(dst2_CEN,'file')
              fprintf('expected file "%s" not found\n', dst2_CEN);
          end
          
          if ~exist(dst_DMN,'file')
              fprintf('expected file "%s" not found\n', dst_DMN);
          end
          
          if ~exist(dst2_DMN,'file')
              fprintf('expected file "%s" not found\n', dst2_DMN);
          end
          
          if ~exist(dst_3Net,'file')
              fprintf('expected file "%s" not found\n', dst2_3Net);
          end
          if ~exist(dst2_3Net,'file')
              fprintf('expected file "%s" not found\n', dst2_3Net);
          end
%         
%         fprintf('expected file "%s" not found\n', src);
%         fprintf('expected file "%s" not found\n', src2);
%    end

%% ALL NETWORKS
    %% timeseries matrix
    ifile = dst_allnet;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_RS_SN_CEN_DMN_SMN_Visual_Net_CAREN_all_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_allnet;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_RS_SN_CEN_DMN_SMN_Visual_Net_CAREN_all_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_SN_CEN_DMN_SMN_Vis_CAREN_labels.csv';
    %ifile_out = 'Connectivity_RS_SM_Visual_Net_CAREN_Tinnitus_Alert_AudioACT_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'));


%% SALIENCE NETWORK
    %% timeseries matrix
    ifile = dst_SN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx' , num2str(isubj) ,'_connectivity_RS_SN_CAREN_Tinnitus_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_SN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx' , num2str(isubj) ,'_connectivity_RS_SN_CAREN_Tinnitus_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_Salience_CAREN_Tinnitus_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'));

%% CENTRAL EXECUTIVE NETWORK
    %% timeseries matrix
    ifile = dst_CEN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_RS_CEN_CAREN_Tinnitus_Alert_AudioACT_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_CEN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_RS_CEN_CAREN_Tinnitus_Alert_AudioACT__static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_CEN_CAREN_AudioACT_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'));

%% DEFAULT MODE NETWORK
    %% timeseries matrix
    ifile = dst_DMN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_RS_DMN_CAREN_Tinnitus_Alert_AudioACT_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_DMN;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_RS_DMN_CAREN_Tinnitus_Alert_AudioACT__static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_DMN_CAREN_AudioACT_labels.csv';
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

%% 3 Networks: SN, CEN, DMN
    %% timeseries matrix
    ifile = dst_3Net;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_RS_SN_CEN_DMN_CAREN_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst2_3Net;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_RS_SN_CEN_DMN_CAREN_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_SN_CEN_DMN_CAREN_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/home/anna.skrzatek/AUDICOG_backup/Analyses/Correl_Matrix/'));


end



%% Save le tableau de correspondance de l'atlas

% %     save(fullfile( output_dir, 'info_networks.mat' ) , 'connectivity_matrix' , 'conn_network',  )  ;
%     writetable(  ts_table, fullfile( output_dir, 'Correspondances_aal3.xlsx' )  ) ; % , 'conn_network',  )  ;

        writetable(  network(2).table , fullfile( output_dir, 'Correspondances_Salience.xlsx' )  ) ; % , 'conn_network',  )  ;

