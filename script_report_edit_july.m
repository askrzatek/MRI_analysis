%function script_report_edit (fspm, fanat, SPM.xCon, SPM.xCon.name, list coordo, list noms coordo)
clc
clear

% %% Init
    if ~exist('par','var')
        par = ''; % for defpar
    end

    addpath /home/anna.skrzatek/

    %% defpar
    defpar.subj_regex = '.*';        %
    defpar.subdir     = '^model_meica';     % name of the working dir
    defpar.file_regex = 'SPM.mat$'; %
    defpar.verbose    = 1;           % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all

    par = complet_struct(par,defpar);
    %
    cd /home/anna.skrzatek/data/
    working_dir = '/home/anna.skrzatek/data/';
    
    load e

    dir_func = e.getSerie('run_ACTIVATION') .toJob;
    e.addSerie('^model','model',1);
    dir_con = e.getSerie('model') .toJob;

    if iscell(dir_con{1})
        nrContrast = length(dir_con);
    else
        nrContrast=1;
    end

    if iscell(dir_func{1})
        nrSubject = length(dir_func);
    else
        nrSubject=1;
    end

    %Coordlist_values = {[0; 0; 0]; [34; -24; 68]; [-34; -24; 68]; [-4; -48; -24]; [0; -24; 10]}; 
    Coordlist.values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); 
    Coordlist.names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'}; 

    %skip = [];
 
    %% LOOP BEGIN
    % 
    %for subj = 1:nrSubject
    %     %% Define the path 
    % % for each subject take the fspm file in model_meika and wms_t1 series,
    % % then take the SPM.xCon.Vcon or SPM.xCon.spmT and its names + ask for ROIs
    % % pointing if ROI coordinates == null (input), otherwise just load the setCoords
    % % for ROIs, check for names of each ROI's coord set
    % % the best would be for it to do as for the run_data : pass by the tags and
    % % the e-object structure for each subject
    %

    % main_dir = fullfile(pwd,'nifti'); % no need
    % stim_dir = fullfile(pwd,'behav'); % no need

    subj = 19;
        model_dir = char(addsuffixtofilenames(e(subj).path,'/model_meica'));
        %fspm_path = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti/2019_02_27_REMINARY_HM_005_V2/meica_model/':
        fspm_path = char(addsuffixtofilenames(model_dir,'/SPM.mat')); % fspm subject's run model

        mkdir(model_dir,'figures_meica'); % add a condition if exists ask to overwrite

        output_dir = char(addsuffixtofilenames(model_dir,'/figures_meica'));
        fspm = {fspm_path};
        t1 = e(subj).getSerie('anat').getVolume('^wms').path(1,:);
        for Ci = 1:nrContrast
        %Ci = 4;

            addpath(char(model_dir))
            load SPM.mat
            
    %% Display subject results for contrast 1

            job_spm_single_results_display(fspm, Ci, t1, Coordlist, output_dir, SPM, working_dir);
        end
        %% end of Contrast loop
    %% end of Subject loop
    %end
% THE END of the job
