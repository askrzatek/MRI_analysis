%% DysPPN connectivity : PPN x Cortical-SubCortical
clear
clc
%main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/cecile.gallea/AGENT10_Daniela';
main_dir = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN';
ROI_dir = fullfile(main_dir,'ROI_PPN_Dys/Cortical_Networks');
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
par.jobname         = 'timeseries_extract_DysFC';
par.display         = 0;
par.redo            = 0;
par.volume          =  e.getSerie('RS').getVolume('s5wts') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound        =  e.getSerie('RS').getRP('rp');
par.mask_threshold  = 0.001;

%% define ROI, using several methods

% par.roi_type.atlas_cat12 = 'aal3';
path_masks_PPN = get_subdir_regex(ROI_dir,'PPN');
path_masks_SMN       = get_subdir_regex(ROI_dir, 'SMN');
path_masks_Assoc_BG       = get_subdir_regex(ROI_dir, 'Assoc_BG');
path_masks_Assoc_Post          = get_subdir_regex(ROI_dir, 'Assoc_Post');
path_masks_Assoc_Ant              = get_subdir_regex(ROI_dir, 'Assoc_Ant');
path_masks_Limbic                   = get_subdir_regex(ROI_dir, 'Limbic');
path_masks_Primary                     = get_subdir_regex(ROI_dir, 'Primary');

%% create file name table from the files list - search for FreeSurfer correspondance tables or complete manually


%% get labels for CAREN ROIs from table according to RSN we are interested in : in this case Salience = RSN01

tab          = readtable(fullfile(ROI_dir,'id_labels_abbrevs_AGEN10.csv'));

PPN.labels  = tab.label(tab.PPN(~isnan(tab.PPN)));
PPN.abbrevs = tab.abbrev(tab.PPN(~isnan(tab.PPN)));
PPN.idx = tab.id(tab.PPN(~isnan(tab.PPN)));

SMN.labels  = tab.label(tab.SMN(~isnan(tab.SMN)));
SMN.abbrevs = tab.abbrev(tab.SMN(~isnan(tab.SMN)));
SMN.idx = tab.id(tab.SMN(~isnan(tab.SMN)));

AssBG.labels  = tab.label(tab.Assoc_BG(~isnan(tab.Assoc_BG)));
AssBG.abbrevs = tab.abbrev(tab.Assoc_BG(~isnan(tab.Assoc_BG)));
AssBG.idx = tab.id(tab.Assoc_BG(~isnan(tab.Assoc_BG)));

AssPost.labels  = tab.label(tab.Assoc_Post(~isnan(tab.Assoc_Post)));
AssPost.abbrevs = tab.abbrev(tab.Assoc_Post(~isnan(tab.Assoc_Post)));
AssPost.idx = tab.id(tab.Assoc_Post(~isnan(tab.Assoc_Post)));

AssAnt.labels  = tab.label(tab.Assoc_Ant(~isnan(tab.Assoc_Ant)));
AssAnt.abbrevs = tab.abbrev(tab.Assoc_Ant(~isnan(tab.Assoc_Ant)));
AssAnt.idx = tab.id(tab.Assoc_Ant(~isnan(tab.Assoc_Ant)));

Limbic.labels  = tab.label(tab.Limbic(~isnan(tab.Limbic)));
Limbic.abbrevs = tab.abbrev(tab.Limbic(~isnan(tab.Limbic)));
Limbic.idx = tab.id(tab.Limbic(~isnan(tab.Limbic)));

Primary.labels  = tab.label(tab.Primary(~isnan(tab.Primary)));
Primary.abbrevs = tab.abbrev(tab.Primary(~isnan(tab.Primary)));
Primary.idx = tab.id(tab.Primary(~isnan(tab.Primary)));

PPN.skip = [];
SMN.skip = [];
AssBG.skip = [];
AssAnt.skip = [];
AssPost.skip = [];
Limbic.skip = [];
Primary.skip = [];

atlas_rois_list_PPN = [];
atlas_rois_list_SMN = [];
atlas_rois_list_AssBG = [];
atlas_rois_list_AssAnt = [];
atlas_rois_list_AssPost = [];
atlas_rois_list_Limbic = [];
atlas_rois_list_Primary = [];

for iroi = 1 : length(PPN.labels)
% nroi = dir(path_masks_Salience_CAREN);
% for iroi = 1: length(nroi)-2
    if ~ismember(PPN.idx(iroi), PPN.skip)
                                %     % path                                                                                                abbrev             description
            atlas_rois_list_PPN = [atlas_rois_list_PPN; char(get_subdir_regex_files(path_masks_PPN,sprintf('.*Reg%d.nii$',PPN.idx(iroi)))), PPN.abbrevs(iroi), PPN.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', PPN.idx(iroi))
    end
end

for iroi = 1 : length(SMN.labels)
% nroi = dir(path_masks_Salience_CAREN);
% for iroi = 1: length(nroi)-2
    if ~ismember(SMN.idx(iroi), SMN.skip)
                                %     % path                                                                                                abbrev             description
            atlas_rois_list_SMN = [atlas_rois_list_SMN; char(get_subdir_regex_files(path_masks_SMN,sprintf('.*Reg%d.nii$',SMN.idx(iroi)))), SMN.abbrevs(iroi), SMN.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', SMN.idx(iroi))
    end
end

for iroi = 1 : length(AssBG.labels)
%nroi = dir(path_masks_CEN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(AssBG.idx(iroi), AssBG.skip)
                                %     % path                                                                                                           abbrev               description
            atlas_rois_list_AssBG = [atlas_rois_list_AssBG; char(get_subdir_regex_files(path_masks_Assoc_BG,sprintf('.*Reg%d.nii$',AssBG.idx(iroi)))), AssBG.abbrevs(iroi), AssBG.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', AssBG.idx(iroi))
    end
end

for iroi = 1 : length(AssAnt.labels)
%nroi = dir(path_masks_CEN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(AssAnt.idx(iroi), AssAnt.skip)
                                %     % path                                                                                                               abbrev                description
            atlas_rois_list_AssAnt = [atlas_rois_list_AssAnt; char(get_subdir_regex_files(path_masks_Assoc_Ant,sprintf('.*Reg%d.nii$',AssAnt.idx(iroi)))), AssAnt.abbrevs(iroi), AssAnt.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', AssAnt.idx(iroi))
    end
end

for iroi = 1 : length(AssPost.labels)
%nroi = dir(path_masks_CEN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(AssPost.idx(iroi), AssPost.skip)
                                %     % path                                                                                                                   abbrev                 description
            atlas_rois_list_AssPost = [atlas_rois_list_AssPost; char(get_subdir_regex_files(path_masks_Assoc_Post,sprintf('.*Reg%d.nii$',AssPost.idx(iroi)))), AssPost.abbrevs(iroi), AssPost.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', AssPost.idx(iroi))
    end
end

for iroi = 1 : length(Limbic.labels)
%nroi = dir(path_masks_CEN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(Limbic.idx(iroi), Limbic.skip)
                                %     % path                                                                                                            abbrev                description
            atlas_rois_list_Limbic = [atlas_rois_list_Limbic; char(get_subdir_regex_files(path_masks_Limbic,sprintf('.*Reg%d.nii$',Limbic.idx(iroi)))), Limbic.abbrevs(iroi), Limbic.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', Limbic.idx(iroi))
    end
end

for iroi = 1 : length(Primary.labels)
%nroi = dir(path_masks_CEN_CAREN);
%for iroi = 1: length(nroi)-2

    if ~ismember(Primary.idx(iroi), Primary.skip)
                                %     % path                                                                                                                abbrev                 description
            atlas_rois_list_Primary = [atlas_rois_list_Primary; char(get_subdir_regex_files(path_masks_Primary,sprintf('.*Reg%d.nii$',Primary.idx(iroi)))), Primary.abbrevs(iroi), Primary.labels(iroi)];
    else
        sprintf('Skipped region %d - not enough signal left after minimal denoising', Primary.idx(iroi))
    end
end

%% we could do the same for other networks : just to have the timeseries extracted per network if we want to compare


%% define all ROIs independently from the atlas

par.roi_type.mask_global = vertcat(atlas_rois_list_PPN, atlas_rois_list_SMN, atlas_rois_list_AssBG, atlas_rois_list_AssAnt, atlas_rois_list_AssPost, atlas_rois_list_Limbic, atlas_rois_list_Primary, {
%     % path                                                           abbrev          description
   }) ;

%% Perform the timeseries' extraction
par.outname = 'Dystonic_PPN_SMN_AssBG_AssAnt_AssPost_Limbic_Primary';
TS = job_extract_timeseries(par);

%% Step 2: create correlation matrix

%% Define some networks : not mandatory

par.network.PPN = PPN.abbrevs;
par.network.SMN = SMN.abbrevs;
par.network.AssBG = AssBG.abbrevs;
par.network.AssAnt = AssAnt.abbrevs;
par.network.AssPost = AssPost.abbrevs;
par.network.Limbic = Limbic.abbrevs;
par.network.Primary = Primary.abbrevs;

%% Create connectivity matrix
TS = job_timeseries_to_connectivity_matrix(TS,par);

%% plot
guidata = plot_resting_state_connectivity_matrix(TS, {e.getSerie('run_RS').getExam().name});

%% Step 3: Perfrom seed-to-voxel correlation (seed == ROI)
par.jobname = 'seed2brain_analysis_DysPPN';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 

%% Step 4: Sauvegarder les matrices sur owncloud pour pouvoir les récupérer sur mac

output_dir = '/home/anna.skrzatek/DYS_PPN/RSFC'; 

for isubj = 1:length(e)
%% part used to rename the files in the subject directories if one forgot to used the par.output_name in the rsfc - therefore for further copying we only need dst files
%    src =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst_cort = fullfile(e(isubj).getSerie('RS').path, 'rsfc/PPN_SMN_AssBG_AssAnt_AssPost_Limbic_Primary__static_conn__timeseries__Dystonic_PPN_SMN_AssBG_AssAnt_AssPost_Limbic_Primary.mat');
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
%     src2 =  fullfile( e(isubj).getSerie('run_RS').path, 'tedana009a1_vt/rsfc/static_conn__timeseries__AUDICOG_CAREN_SN_CEN_DMN_AudioACT.mat') ;
    dst2_cort = fullfile(e(isubj).getSerie('RS').path,'rsfc/timeseries__Dystonic_PPN_SMN_AssBG_AssAnt_AssPost_Limbic_Primary.mat');
        
%     if exist(src, 'file')
%         movefile(src, dst);
%         movefile(src2, dst2);
%     else
          if ~exist(dst_cort,'file')
              fprintf('expected file "%s" not found\n', dst_cort);
          end
          
          if ~exist(dst2_cort,'file')
              fprintf('expected file "%s" not found\n', dst2_cort);
          end
          
          if ~exist(dst_cereb,'file')
              fprintf('expected file "%s" not found\n', dst_cereb);
          end
          
          if ~exist(dst2_cereb,'file')
              fprintf('expected file "%s" not found\n', dst2_cereb);
          end
          
%         fprintf('expected file "%s" not found\n', src);
%         fprintf('expected file "%s" not found\n', src2);
%    end

%% CORTICAL NETWORKS
    %% timeseries matrix
    ifile = dst2_cort;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name, '_idx' , num2str(isubj) ,'_connectivity_RS_DysPPN_Cortical_timeseries.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/'), rs_cmat.timeseries)

    %% static conn    
    ifile = dst_cort;
    rs_cmat = load(ifile)   ;

    ifile_out = [e(isubj).name , '_idx', num2str(isubj) ,'_connectivity_RS_DysPPN_Cortical_static_conn.csv' ]; 
    dlmwrite(addprefixtofilenames(ifile_out,'/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/'), rs_cmat.static_connectivity_matrix)
    
    ifile_out = 'Connectivity_RS_DysPPN_labels.csv';
    writetable(rs_cmat.ts_table, addprefixtofilenames(ifile_out,'/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/'));
end

%% Save le tableau de correspondance de l'atlas

% %     save(fullfile( output_dir, 'info_networks.mat' ) , 'connectivity_matrix' , 'conn_network',  )  ;
%     writetable(  ts_table, fullfile( output_dir, 'Correspondances_aal3.xlsx' )  ) ; % , 'conn_network',  )  ;

        writetable(  network(2).table , fullfile( output_dir, 'Correspondances_Salience.xlsx' )  ) ; % , 'conn_network',  )  ;

