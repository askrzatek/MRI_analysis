%function script_report_edit
%function script_report_edit (fspm, fanat, SPM.xCon, SPM.xCon.name, list coordo, list noms coordo)
%clc
%clear
% %% Init
    if ~exist('par','var')
        par = ''; % for defpar
    end
%%  INPUTS
%   FSPM : cell variable defining the path to the SPM.mat file
%   CI : integer variable with the number of the contrast to display
%   T1 : character variable defining the path to the T1 MRI volume to
%   display
%   COORDLIST : structure variable defining ROI's names and coordinates to
%   display
%   OUTPUT_DIR : directory path where you want your figures to be saved
%   SPMFILE : structure variable defining SPM model
%   WD : working directory path where to find this script
%
%%
    addpath(wd)
    addpath(output_dir)
    spm('defaults','FMRI');
        matlabbatch{1}.spm.stats.results.spmmat = fspm;
        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        
        matlabbatch{1}.spm.stats.results.conspec.contrasts = Ci; % learn how to choose multiple contrasts or do a Contrast loop
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
        matlabbatch{1}.spm.stats.results.conspec.extent = 5;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.print = false;
        matlabbatch{1}.spm.stats.results.write.none = 1;
        spm_jobman('run',matlabbatch);
        
        
        hReg = findobj('Tag','hReg');
        xSPM = evalin('base', 'xSPM');
        
    %% Load section == 3DT1 normalized


 %       try
    %   hReg = findobj('Tag','hReg');
            
            spm_sections(xSPM,hReg,t1);

            %% Go to the position : aka ROI choice
            for iROI = 1 : length(Coordlist.values)

                hX    = findobj('Tag','hX'); % pick on, because we need the upper USerData
                hFxyz = get(hX,'UserData');
                UD    = get(hFxyz,'UserData');


                nxyz = Coordlist.values(:,iROI);
            %nxyz = rand(3,1)*30 -10;

                spm_XYZreg('SetCoords',nxyz,UD.hReg,hFxyz);
                spm_results_ui('UpdateSPMval',UD)
                
            %% Save figure

                Fgraph = spm_figure('GetWin','Graphics');
                
                %Fgraph.Children(2).Children(7).Children(2) - display
                %coordinates display next to images maybe
                % dCoord = findobj('Tag','OVmenu_Coordinates');

                graphName = char(addsuffixtofilenames(addsuffixtofilenames(Coordlist.names(iROI),'_'),spmfile.xCon(Ci).name));
                
                global st
                X  = frame2im(getframe(st.fig));
                sz = size(X);
                sz = [sz(1) sz(1) sz(2) sz(2)];
                p1 = get(st.vols{1}.ax{1}.ax,'Position');
                p2 = get(st.vols{1}.ax{3}.ax,'Position');
                a  = [p1(1) p1(2)  p2(1)+p2(3)-p1(1) p2(2)+p2(4)-p1(2)] + 0.005*[-1 -1 2 2];
                a  = max(min(a,1),0);
                sz = ([1-a(2)-a(4),1-a(2),a(1),a(1)+a(3)] .* (sz-1)) + 1;
                sz = round(sz);
                X  = X(sz(1):sz(2),sz(3):sz(4),:);
                cd (output_dir);
                imwrite(X,graphName,'png');
                fprintf('Saving image as:\n');
                fprintf('  %s\n',spm_file(graphName,'link','web(''%s'')'));
                %saveas(Fgraph, graphName, 'jpeg') % choose the name according to the Subject & the Contrast & the ROI
            end
%        end
    
    
    addpath /home/anna.skrzatek/

    %% defpar
    defpar.subj_regex = '.*';        %
    defpar.subdir     = '^model_meica';     % name of the working dir
    defpar.file_regex = 'SPM.mat$'; %
    defpar.verbose    = 1;           % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all

    par = complet_struct(par,defpar);
    %
    %cd /home/anna.skrzatek/data/

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
    Coordlist_values = cat(2,[0.0; 0.0; 0.0], [34.0; -24.0; 68.0], [-34.0; -24.0; 68.0], [-4.0; -48.0; -24.0], [0.0; -24.0; 10.0]); 
    Coordlist_names = {'centre'; 'RIGHT SM'; 'LEFT SM'; 'CEREBELLUM'; 'MOTOR BI'}; 

    %skip = [];

    % % start de loop by defining the path of the directory to fetch your db -
    % % or just loading the e.mat could do the job
    % % define input & output folder (create the folder if needed but the best
    % % would be to put the img in the model_meica folder)
    % % needs files and directories defined in firstlevel.m : e + main_dir +
    % % dir_func + model_dir + dir_anat
    % % stim_dir + stim_files
    % 
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

    subj = 1;
        model_dir = char(addsuffixtofilenames(e(subj).path,'/model_meica'));
        %fspm_path = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti/2019_02_27_REMINARY_HM_005_V2/meica_model/':
        fspm_path = char(addsuffixtofilenames(model_dir,'/SPM.mat')); % fspm subject's run model

        mkdir(model_dir,'figures_meica'); % add a condition if exists ask to overwrite

        output_dir = char(addsuffixtofilenames(model_dir,'/figures_meica'));
        fspm = {fspm_path};
        
    % % How do I use it ? is it to access different files by the regex - tree architecture
    % 
    % % defpar.subj_regex = '.*';        %
    % % defpar.subdir     = '^meica_model';     % name of the working dir
    % % defpar.file_regex = 'SPM.mat$'; %
    % % defpar.verbose    = 1;           % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all
    % % 
    % % par = complet_struct(par,defpar);
    %     
        t1 = e(subj).getSerie('anat').getVolume('^wms').path;
    %     %t1 = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti/2019_02_27_REMINARY_HM_005_V2/S03_t1mpr_S256_0_8iso_p2/wms_S03_t1mpr_S256_0_8iso_p2.nii';
    % 
    %     % job_contrast2jpeg(fspm, t1, SPM.xCon, SPM.xCon.name, list coordo, list noms coordo )
        % if noContrast
        %for Ci = 1:nrContrast
        Ci = 4;
          
            %% end of ROI loop
%    end
        %% end of Contrast loop
    %% end of Subject loop
%end
% THE END of the job