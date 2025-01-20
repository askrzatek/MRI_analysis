%% DysPPN connectivity : PPN x Cortical-SubCortical
clear
clc
%main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/cecile.gallea/AGENT10_Daniela';
main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN';
ROI_dir = fullfile(main_dir,'DYS_PPN_ROIs/Cortical_Networks');
data_dir = fullfile(main_dir,'AGENT_10');

cd (main_dir)

%% Load the data & structure it

e = exam(main_dir, ".*Subj.*")
e.addSerie(".*RS","RS",1)
e.getSerie("RS").addVolume('s5wts_OC','s5wts_OC',1)
%load('e_nonchir.mat'); 


%% Step 1: faire tourner les calculs de connnectivité

% define input volumes and confounds
clear par
par.run             = 1;
par.sge             = 0;
par.mem             = '16G';
par.jobname         = 'timeseries_extract_DysFC';
par.display         = 0;
par.redo            = 0;
par.volume          =  e.getSerie('run_RS').getVolume('s5wts_OC') ;  % Choisir le s5 ou s8 (en fonction de la taille des structures etudiees par exemple...)
par.confound        =  e.getSerie('run_RS').getRP('multiple_regressors');
par.mask_threshold  = 0.001;

%% define ROI, using several methods

% par.roi_type.atlas_cat12 = 'aal3';
path_masks_PPN = get_subdir_regex(ROI_dir,"PPN")
path_masks_SMN       = get_subdir_regex(ROI_dir, 'SMN');
path_masks_Assoc_BG       = get_subdir_regex(ROI_dir, 'Assoc_BG');
path_masks_Assoc_Post          = get_subdir_regex(ROI_dir, 'Assoc_Post');
path_masks_Assoc_Ant              = get_subdir_regex(ROI_dir, 'Assoc_Post');
path_masks_Limbic                   = get_subdir_regex(ROI_dir, 'Assoc_Ant');
path_masks_Primary                     = get_subdir_regex(ROI_dir, 'Primary');

%% create file name table from the files list - search for FreeSurfer correspondance tables or complete manually


%% get labels for CAREN ROIs from table according to RSN we are interested in : in this case Salience = RSN01

tab          = readtable(fullfile(ROI_dir,'id_labels_abbrevs_AGEN10.csv'));

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

SMN.skip = [];
AssBG.skip = [];
AssAnt.skip = [];
AssPost.skip = [];
Limbic.skip = [];
Primary.skip = [];

atlas_rois_list_SMN = [];
atlas_rois_list_AssBG = [];
atlas_rois_list_AssAnt = [];
atlas_rois_list_AssPost = [];
atlas_rois_list_Limbic = [];
atlas_rois_list_Primary = [];


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

par.roi_type.mask_global = vertcat(atlas_rois_list_SMN, atlas_rois_list_AssBG, atlas_rois_list_AssAnt, atlas_rois_list_AssPost, atlas_rois_list_Limbic, atlas_rois_list_Primary, {
%     % path                                                           abbrev          description
   }) ;

%% Perform the timeseries' extraction
par.outname = 'AUDICOG_Dystonic_PPN_SMN_AssBG_AssAnt_AssPost_Limbic_Primary';
TS = job_extract_timeseries(par);

%% Step 2: create correlation matrix

%% Define some networks : not mandatory


par.network.SMN = SMN.abbrevs;
par.network.AssBG = AssBG.abbrevs;
par.network.AssAnt = AssAnt.abbrevs;
par.network.AssPost = AssPost.abbrevs;
par.network.Limbic = Limbic.abbrevs;
par.network.Primary = Primary.abbrevs;

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
par.jobname = 'seed2brain_analysis_DysPPN';
TS = job_timeseries_to_connectivity_seedbased(TS,par);
% 

%% Step 4: Sauvegarder les matrices sur owncloud pour pouvoir les récupérer sur mac

output_dir = '/home/anna.skrzatek/'; 

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

