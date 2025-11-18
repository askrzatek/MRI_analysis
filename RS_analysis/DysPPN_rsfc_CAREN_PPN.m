%% DysPPN connectivity : PPN x Cortical-SubCortical
clear
clc
%main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/cecile.gallea/AGENT10_Daniela';
main_dir = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN';
ROI_dir = fullfile(main_dir,'ROI_PPN_Dys/PPN_SMN');

data_dir = fullfile(main_dir,'AGENT10');

cd (main_dir)
addpath(main_dir)

%% Load the data & structure it

e = exam(data_dir, '.*AGENT10.*')
e.addSerie('.*MINEA.*','tedana.*','RS',1)
e.getSerie('RS').addVolume('^s5wts_OC.nii','s5wts',1)
%e.getSerie('RS').addVolume('rp_spm','rp',1)
e.getSerie('RS').addRP('rp_spm','rp',1)
e.explore

%% Step 1: faire tourner les calculs de connnectivité

% define input volumes and confounds
clear par
par.run             = 1;
par.sge             = 0;
par.mem             = '16G';
par.jobname         = 'timeseries_extract_DYS_PPN_CAREN_all';
par.display         = 0;
par.redo            = 0;
par.volume          =  e.getSerie('RS').getVolume('s5wts') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound        =  e.getSerie('RS').getRP('rp');
par.mask_threshold  = 0.001;

%% define ROI, using several methods

% par.roi_type.atlas_cat12 = 'aal3';
path_masks_PPN      = get_subdir_regex(ROI_dir,'PPN');
path_masks_SMN_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN03';
path_masks_Salience_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN01';
path_masks_CEN_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN02';
path_masks_DMN_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN04';
path_masks_VIS_CAREN = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Networks_Masks/CAREN/RSN05';

%% create file name table from the files list - search for FreeSurfer correspondance tables or complete manually


%% get labels for CAREN ROIs from table according to RSN we are interested in : in this case Salience = RSN01

tab          = readtable(fullfile(ROI_dir,'id_labels_abbrevs_AGEN10_PPN_SMN.csv'));

PPN.labels  = tab.label(tab.PPN(~isnan(tab.PPN)));
PPN.abbrevs = tab.abbrev(tab.PPN(~isnan(tab.PPN)));
PPN.idx = tab.id(tab.PPN(~isnan(tab.PPN)));

tab          = readtable(fullfile(ROI_dir,'aal_id_labels_abbrevs_CAREN.csv'));

SM_net.labels  = strcat(tab.label_AAL(tab.smn(~isnan(tab.smn))),'_SMN');
SM_net.abbrevs = strcat(tab.abbrev_AAL(tab.smn(~isnan(tab.smn))),'_SMN');
SM_net.idx = tab.id_AAL(tab.smn(~isnan(tab.smn)));
salience_net.labels  = strcat(tab.label_AAL(tab.salience(~isnan(tab.salience))),'_SN');
salience_net.abbrevs = strcat(tab.abbrev_AAL(tab.salience(~isnan(tab.salience))),'_SN');
salience_net.idx = tab.id_AAL(tab.salience(~isnan(tab.salience)));
cEx_net.labels  = strcat(tab.label_AAL(tab.cen(~isnan(tab.cen))),'_CEN');
cEx_net.abbrevs = strcat(tab.abbrev_AAL(tab.cen(~isnan(tab.cen))),'_CEN');
cEx_net.idx = tab.id_AAL(tab.cen(~isnan(tab.cen)));
DM_net.labels  = strcat(tab.label_AAL(tab.dmn(~isnan(tab.dmn))),'_DMN');
DM_net.abbrevs = strcat(tab.abbrev_AAL(tab.dmn(~isnan(tab.dmn))),'_DMN');
DM_net.idx = tab.id_AAL(tab.dmn(~isnan(tab.dmn)));
Vis_net.labels  = strcat(tab.label_AAL(tab.vis(~isnan(tab.vis))),'_VIS');
Vis_net.abbrevs = strcat(tab.abbrev_AAL(tab.vis(~isnan(tab.vis))),'_VIS');
Vis_net.idx = tab.id_AAL(tab.vis(~isnan(tab.vis)));

PPN.skip       = [];
SM_net.skip = [];
% SM_net.skip = [6];
salience_net.skip = [];
cEx_net.skip = [];
% cEx_net.skip = [57];
DM_net.skip = [];
% DM_net.skip = [51,93];
Vis_net.skip = [];
% Vis_net.skip = [93];

atlas_rois_list_PPN = [];
atlas_rois_list_SMN = [];
atlas_rois_list_SN = [];
atlas_rois_list_CEN = [];
atlas_rois_list_DMN = [];
atlas_rois_list_Vis = [];

for iroi = 1 : length(PPN.labels)
    if ~ismember(PPN.idx(iroi), PPN.skip)
                                %     % path                                                                                                abbrev             description
            atlas_rois_list_PPN = [atlas_rois_list_PPN; char(get_subdir_regex_files(path_masks_PPN,sprintf('.*Reg%d.nii$',PPN.idx(iroi)))), PPN.abbrevs(iroi), PPN.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', PPN.idx(iroi))
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

for iroi = 1 : length(salience_net.labels)
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

for iroi = 1 : length(Vis_net.labels)
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

par.roi_type.mask_global = vertcat(atlas_rois_list_PPN, atlas_rois_list_SMN, atlas_rois_list_SN, atlas_rois_list_CEN, atlas_rois_list_DMN, atlas_rois_list_Vis, {
%     % path                                                           abbrev          description
   }) ;

%% Perform the timeseries' extraction
% par.outname = 'Dystonic_PPN_SMN_AssBG_AssAnt_AssPost_Limbic_Primary';
par.outname = 'Dystonic_PPN_SMN_SN_CEN_DMN_VIS_CAREN';
TS = job_extract_timeseries(par);

%% Step 2: create correlation matrix

%% Define some networks : not mandatory

par.network.PPN       = PPN.abbrevs;
par.network.SMN       = SM_net.abbrevs;
par.network.SN        = salience_net.abbrevs;
par.network.CEN       = cEx_net.abbrevs;
par.network.DMN       = DM_net.abbrevs;
par.network.Vis       = Vis_net.abbrevs;

%% Create connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

%% plot
guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% Step 3: Perfrom seed-to-voxel correlation (seed == ROI)
par.jobname = 'seed2brain_analysis_DysPPN_CAREN';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 

%% Step 4: Sauvegarder les matrices sur owncloud pour pouvoir les récupérer sur mac

output_dir = '/home/anna.skrzatek/DYS_PPN/RSFC'; 

for isubj = 1:length(e)
%% part used to rename the files in the subject directories if one forgot to used the par.output_name in the rsfc - therefore for further copying we only need dst files
    dst_CAREN = fullfile(e(isubj).getSerie('RS').path, 'rsfc/PPN_SMN_SN_CEN_DMN_VIS_CAREN__static_conn__timeseries__Dystonic_PPN_SMN_SN_CEN_DMN_VIS_CAREN.mat');
    dst2_CAREN = fullfile(e(isubj).getSerie('RS').path,'rsfc/timeseries__Dystonic_PPN_SMN_SN_CEN_DMN_VIS_CAREN.mat');
    if ~exist(dst_CAREN,'file')
        fprintf('expected file "%s" not found\n', dst_CAREN);
    end
    if ~exist(dst2_CAREN,'file')
        fprintf('expected file "%s" not found\n', dst2_CAREN);
    end
          
    %% timeseries matrix
    ifile = dst2_cort;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_RS_DysPPN_SMN_SN_CEN_DMN_VIS_CAREN_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst_cort;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_RS_DysPPN_SMN_SN_CEN_DMN_VIS_CAREN_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_DysPPN_SMN_SN_CEN_DMN_VIS_CAREN_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/'));
end

%% Save le tableau de correspondance de l'atlas

% %     save(fullfile( output_dir, 'info_networks.mat' ) , 'connectivity_matrix' , 'conn_network',  )  ;
%     writetable(  ts_table, fullfile( output_dir, 'Correspondances_aal3.xlsx' )  ) ; % , 'conn_network',  )  ;

        writetable(  network(2).table , fullfile( output_dir, 'Correspondances_Salience.xlsx' )  ) ; % , 'conn_network',  )  ;

