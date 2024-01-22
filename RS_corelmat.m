%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS Correlation Matrix script
%% PARKGAME II
%% January 2022
%% A.SKRZATEK
%% inspired from S.OUARAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
%% fetch Input & Output dirs
main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/';
cd (main_dir)

%% either existing e-file from RS processing
load e
%% or new e-object

V1 = exam(main_dir,'PARKGAME.*.V1_[a,c]$');
V1 = V1(1:10) + V1(13:14) + V1(16:18);
V1.addSerie('S0.*RS$','model$','model_2$','model',1);


V2 = exam(main_dir,'PARKGAME.*.V2_[a,c]$');
V2.addSerie('S0.*RS$','model$','model_2$','model',1)

%%
patient_regex_V1 = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII_040.*RE.*_a','PARKGAMEII_042.*RS.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII_046.*HJ.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c','PARKGAMEII_053.*LM.*_a'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_053_LM_a_V1'};
patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'};
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};

first_dir = fullfile(main_dir,'/firstlevel_RS');
input_dir = get_subdir_regex(first_dir,patient_list);

% ROIs = {'Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','Postcentral_L','Postcentral_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10','Caudate_aal_L','Caudate_aal_R','Thalamus_aal_L','Thalamus_aal_R','Putamen_aal_L','Putamen_aal_R','Pallidum_aal_L','Pallidum_aal_R','PPN_YEB_L','PPN_YEB_R'};
%ROIs = {'M1_face_L','M1_face_R','M1_foot_L','M1_foot_R','M1_hand_L','M1_hand_R'};
%dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_pariet_mot_premot_cereb_BG_PPN';
%fileROI = cellstr(char(gfile(dirROI,ROIs))); %fileROI = remove_regex(fileROI,'T1');

dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/FOG_APA_network';
dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/Cognitive_circuit';
dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/Cueing';
dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/Kinect';
dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/Ordi';
[~, roi_group] = fileparts(dirROI);
fileROI = cellstr(char(gfile(dirROI,'.*.nii'))); %fileROI = remove_regex(fileROI,'T1');

dirROI  = get_subdir_regex('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/Grouped_ROIs/','.*');
roi_group = 'Basal_Ganglia-Cueing-Motor-AOT';
fileROI = cellstr(char(vertcat(gfile(dirROI{3},'.*.nii'),gfile(dirROI{1},'.*.nii'),gfile(dirROI{8},'.*.nii'),gfile(dirROI{4},'.*.nii'),gfile(dirROI{7},'.*.nii'),gfile(dirROI{5},'.*.nii'),gfile(dirROI{2},'.*.nii'),gfile(dirROI{6},'.*.nii'))));

ROIs = cell(1,length(fileROI));
for iROI = 1 : length(fileROI)
    [~,roi_name] = fileparts(fileROI{iROI});
    ROIs{iROI} = roi_name;
end

roi_labels = {};

for ir = 1 : length(ROIs)
    [~,roi_name] = fileparts(ROIs{1,ir,1});
    if roi_name(end)=='5'
        roi_name = roi_name(1:end-10);
    end
    
    roi_labels = vertcat(roi_labels,roi_name);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verif VOI
%% get VOIs for each run
nRun = 1;
nRun   = 2;
nSubj = length(input_dir);
%    nROI  = length(ROIs)+1;
nROI = length(fileROI);

voi_file = cell(nSubj,nROI,nRun);

for subj = 1 : nSubj
    for ir = 1 : nROI

       [~,roi_name] = fileparts(fileROI{ir});
       %roi_name = roi_name(1:end-10);

       for iRun = 1 : nRun
           %voi_file{subj,ir,iRun} = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
           voi_file(subj,ir,iRun) = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s_%d.mat',roi_name,iRun));
       end
       
%        model_dir = gpath(V1(subj).getSerie('model'));
%        voi_file(subj,ir,1) = get_subdir_regex_files(model_dir,sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
%        model_dir = gpath(V2(subj).getSerie('model'));
%        voi_file(subj,ir,2) = get_subdir_regex_files(model_dir,sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
       
    end
%     for iRun = 1 : nRun
%        %voi_file{subj,nROI,iRun} = get_subdir_regex_files(input_dir{subj},sprintf('VOI_PCC_%d.mat',iRun));
%        voi_file(subj,nROI,iRun) = get_subdir_regex_files(input_dir{subj},sprintf('VOI_PCC_%d.mat',iRun));
%     end   
end
voi = voi_file;
%fmat = get_subdir_regex_files(input_dir,ROIs);
%fmat{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for subj = 1 : nSubj
    for iRun = 1 : nRun
        f1name = sprintf('/V%d_corr_matrix_%s_%s.csv',iRun,roi_group,patient_list{subj}(12:end));
        f2name = sprintf('/V%d_pval_corr_matrix_%s_%s.csv',iRun,roi_group,patient_list{subj}(12:end));
        fout1 = addprefixtofilenames(input_dir{subj},f1name);
        fout2 = addprefixtofilenames(input_dir{subj},f2name);
        fout1bis = fullfile(first_dir,f1name);
        fout2bis = fullfile(first_dir,f2name);

        sujroi = voi(subj,:,iRun);
        time_course_matrix = {};
        parfor i = 1 : length(sujroi)
            l = load(sujroi{i});
            time_course_matrix{i} = l.Y;
        end
        time_course_matrix = cell2mat(time_course_matrix);

        [cor_mat, pval] = corr(time_course_matrix);

%% text output        

        [~, roiname] = get_parent_path(sujroi);
        roiname = nettoie_dir(change_file_extension(roiname,''))

        Tcor = array2table(cor_mat, 'VariableNames', roiname);
        Tpval = array2table(pval, 'VariableNames', roiname);

        writetable(Tcor,fout1)
        writetable(Tpval,fout2)
        writetable(Tcor,fout1bis)
        writetable(Tpval,fout2bis)

%% visual output        
        t = sprintf('Correlation matrix V%d :',iRun);
        %roilabels = vertcat(ROIs(:),{'PCC'});
        
        figure; imagesc(cor_mat)
        title([t, patient_list_V1{subj}], 'Interpreter', 'none')
        ax = gca;
        x = 0.5:54.5;
        y = 0.5:54.5;
        ax.XTick = x;
        ax.YTick = y;
        ax.XTickLabelRotation = 90;
        ax.YTickLabelRotation = 0;
        ax.XTickLabel = roi_labels;
        ax.YTickLabel = roi_labels;
        ax.TickLabelInterpreter = 'none';
        savefig(sprintf('%sV%d_cormat_%s.fig',input_dir{subj},iRun,roi_group))
        close

    end
    %% Could be a good place to get V1 and V2 table values and make a difference (individual deltas)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORELATION MATRIX PER ROIs DIRECTORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fetch Input & Output dirs
main_dir = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/';
cd (main_dir)

patient_regex_V1 = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c','PARKGAMEII_040.*RE.*_a$','PARKGAMEII_042.*RS.*_a$','PARKGAMEII_046.*HJ.*_c$','PARKGAMEII_053.*LM.*_a$'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_053_LM_a_V1'};

patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'};
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};

input_dir = get_subdir_regex(fullfile(main_dir,'/firstlevel_RS'),patient_list);

dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks';
%dirROI = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks/FOG_APA_network';
ROI_mask = get_subdir_regex(dirROI,roi_group);
ROI_mask = roi_group;

for imask = 1: length(ROI_mask)
    roi_labels = {};
    mask_path = strsplit(ROI_mask{imask},'/');
    mask_name = mask_path{end-1};
    mask_name = ROI_mask;

    fileROI{imask} = cellstr(char(gfile(ROI_mask{imask},'.*')));
    for ir = 1 : length(fileROI{imask})

           [~,roi_name] = fileparts(fileROI{imask}{ir});
           roi_name = roi_name(1:end-10);
           roi_labels = vertcat(roi_labels,roi_name);
    end
            
    
    t = sprintf('V1 Mean correlation matrix %s KINECT:',mask_name);
    roilabels = vertcat(roi_labels(:));

    cd(main_dir)
%    tab =
%    csvread(sprintf('%s/RS_corel_matrix/means/Kinect_%s_delta_mean_relative.csv',main_dir,mask_name(1:end-5)),1,1);
%    %mask
    tab = csvread(sprintf('%s/RS_corel_matrix/means/V1_Kinect_%s_mean.csv',main_dir,mask_name),2,2); %network
    clrLim = [-1 1];
    diamLim = [0.3, 1];

    figure()
    imagesc(tab)
    title([t], 'Interpreter', 'none');
    colormap(gca,'jet');
    c = colorbar;
    caxis([-1 1]);
%            c.Label.String = "Correlation d'activité au repos entre ROI";
    axis equal
    axis tight

    ax = gca;
    ax.XTickLabelRotation = 90;
    ax.YTickLabelRotation = 0;
    ax.XTick = [1:1:length(roilabels)];
    ax.YTick = [1:1:length(roilabels)];
    ax.XTickLabel = roilabels;
    ax.YTickLabel = roilabels;
    ax.TickLabelInterpreter = 'none';


    savefig(sprintf('Mean_V1_Kinect_cormat_%s.fig',mask_name))
    close

  %% BUBBLE VERSION
  clear ax
  clear fh
% Compute center of each circle 
% This assumes the x and y values were not entered in imagesc()
  x = 1 : 1 : size(tab,2); % x edges
  y = 1 : 1 : size(tab,1); % y edges
  [xAll, yAll] = meshgrid(x,y); 
% Set color of each rectangle
% Set color scale
  cmap = parula(256);
  cmap = jet;
  Cscaled = (tab - clrLim(1))/range(clrLim); % always [0:1]
  colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size in the same way we scale the color
  diamSize = Cscaled * range(diamLim) + diamLim(1); 
% diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling

% Create figure
            fh = figure(); 
            ax = axes(fh);
            title([t], 'Interpreter', 'none');
            hold(ax,'on')
            colormap(ax,'jet');
            % Create circles
            theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
            % If you get an error in this line below, "Index in position 1 is invalid",
            % that probably means you didn't set clrLim correctly. 
            h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
                diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
            axis(ax,'equal')
            axis(ax,'tight')
            set(ax,'YDir','Reverse')
            colorbar()
            caxis(clrLim);


            ax = gca;
            ax.XTickLabelRotation = 0;
            ax.YTickLabelRotation = 0;
            ax.XTickLabel = roilabels;
            ax.YTickLabel = roilabels;
            ax.TickLabelInterpreter = 'none';
            
            
            savefig(sprintf('Mean_V1_Kinect_cormat_%s_bubble.fig',mask_name))
            close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = sprintf('V2 Mean correlation matrix %s KINECT:',mask_name);
    roilabels = vertcat(roi_labels(:));

    cd(main_dir)
%    tab =
%    csvread(sprintf('%s/RS_corel_matrix/means/Kinect_%s_delta_mean_relative.csv',main_dir,mask_name(1:end-5)),1,1);
%    %mask
    tab = csvread(sprintf('%s/RS_corel_matrix/means/V2_Kinect_%s_mean.csv',main_dir,mask_name),2,2); %network
    clrLim = [-1 1];
    diamLim = [0.3, 1];

    figure()
    imagesc(tab)
    title([t], 'Interpreter', 'none');
    colormap(gca,'jet');
    c = colorbar;
    caxis([-1 1]);
%            c.Label.String = "Correlation d'activité au repos entre ROI";
    axis equal
    axis tight

    ax = gca;
    ax.XTickLabelRotation = 90;
    ax.YTickLabelRotation = 0;
    ax.XTick = [1:1:length(roilabels)];
    ax.YTick = [1:1:length(roilabels)];
    ax.XTickLabel = roilabels;
    ax.YTickLabel = roilabels;
    ax.TickLabelInterpreter = 'none';


    savefig(sprintf('Mean_V2_Kinect_cormat_%s.fig',mask_name))
    close

  %% BUBBLE VERSION
  clear ax
  clear fh
% Compute center of each circle 
% This assumes the x and y values were not entered in imagesc()
  x = 1 : 1 : size(tab,2); % x edges
  y = 1 : 1 : size(tab,1); % y edges
  [xAll, yAll] = meshgrid(x,y); 
% Set color of each rectangle
% Set color scale
  cmap = parula(256);
  cmap = jet;
  Cscaled = (tab - clrLim(1))/range(clrLim); % always [0:1]
  colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size in the same way we scale the color
  diamSize = Cscaled * range(diamLim) + diamLim(1); 
% diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling

% Create figure
            fh = figure(); 
            ax = axes(fh);
            title([t], 'Interpreter', 'none');
            hold(ax,'on')
            colormap(ax,'jet');
            % Create circles
            theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
            % If you get an error in this line below, "Index in position 1 is invalid",
            % that probably means you didn't set clrLim correctly. 
            h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
                diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
            axis(ax,'equal')
            axis(ax,'tight')
            set(ax,'YDir','Reverse')
            colorbar()
            caxis(clrLim);


            ax = gca;
            ax.XTickLabelRotation = 0;
            ax.YTickLabelRotation = 0;
            ax.XTickLabel = roilabels;
            ax.YTickLabel = roilabels;
            ax.TickLabelInterpreter = 'none';
            
            
            savefig(sprintf('Mean_V2_Kinect_cormat_%s_bubble.fig',mask_name))
            close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = sprintf('V1 Mean correlation matrix %s ORDI:',mask_name);
    roilabels = vertcat(roi_labels(:));

    cd(main_dir)
%    tab =
%    csvread(sprintf('%s/RS_corel_matrix/means/Kinect_%s_delta_mean_relative.csv',main_dir,mask_name(1:end-5)),1,1);
%    %mask
    tab = csvread(sprintf('%s/RS_corel_matrix/means/V1_Ordi_%s_mean.csv',main_dir,mask_name),2,2); %network
    clrLim = [-1 1];
    diamLim = [0.3, 1];

    figure()
    imagesc(tab)
    title([t], 'Interpreter', 'none');
    colormap(gca,'jet');
    c = colorbar;
    caxis([-1 1]);
%            c.Label.String = "Correlation d'activité au repos entre ROI";
    axis equal
    axis tight

    ax = gca;
    ax.XTickLabelRotation = 90;
    ax.YTickLabelRotation = 0;
    ax.XTick = [1:1:length(roilabels)];
    ax.YTick = [1:1:length(roilabels)];
    ax.XTickLabel = roilabels;
    ax.YTickLabel = roilabels;
    ax.TickLabelInterpreter = 'none';


    savefig(sprintf('Mean_V1_Ordi_cormat_%s.fig',mask_name))
    close

  %% BUBBLE VERSION
  clear ax
  clear fh
% Compute center of each circle 
% This assumes the x and y values were not entered in imagesc()
  x = 1 : 1 : size(tab,2); % x edges
  y = 1 : 1 : size(tab,1); % y edges
  [xAll, yAll] = meshgrid(x,y); 
% Set color of each rectangle
% Set color scale
  cmap = parula(256);
  cmap = jet;
  Cscaled = (tab - clrLim(1))/range(clrLim); % always [0:1]
  colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size in the same way we scale the color
  diamSize = Cscaled * range(diamLim) + diamLim(1); 
% diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling

% Create figure
            fh = figure(); 
            ax = axes(fh);
            title([t], 'Interpreter', 'none');
            hold(ax,'on')
            colormap(ax,'jet');
            % Create circles
            theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
            % If you get an error in this line below, "Index in position 1 is invalid",
            % that probably means you didn't set clrLim correctly. 
            h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
                diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
            axis(ax,'equal')
            axis(ax,'tight')
            set(ax,'YDir','Reverse')
            colorbar()
            caxis(clrLim);


            ax = gca;
            ax.XTickLabelRotation = 0;
            ax.YTickLabelRotation = 0;
            ax.XTickLabel = roilabels;
            ax.YTickLabel = roilabels;
            ax.TickLabelInterpreter = 'none';
            
            
            savefig(sprintf('Mean_V1_Ordi_cormat_%s_bubble.fig',mask_name))
            close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = sprintf('V2 Mean correlation matrix %s ORDI:',mask_name);
    roilabels = vertcat(roi_labels(:));

    cd(main_dir)
%    tab =
%    csvread(sprintf('%s/RS_corel_matrix/means/Kinect_%s_delta_mean_relative.csv',main_dir,mask_name(1:end-5)),1,1);
%    %mask
    tab = csvread(sprintf('%s/RS_corel_matrix/means/V2_Ordi_%s_mean.csv',main_dir,mask_name),2,2); %network
    clrLim = [-1 1];
    diamLim = [0.3, 1];

    figure()
    imagesc(tab)
    title([t], 'Interpreter', 'none');
    colormap(gca,'jet');
    c = colorbar;
    caxis([-1 1]);
%            c.Label.String = "Correlation d'activité au repos entre ROI";
    axis equal
    axis tight

    ax = gca;
    ax.XTickLabelRotation = 90;
    ax.YTickLabelRotation = 0;
    ax.XTick = [1:1:length(roilabels)];
    ax.YTick = [1:1:length(roilabels)];
    ax.XTickLabel = roilabels;
    ax.YTickLabel = roilabels;
    ax.TickLabelInterpreter = 'none';


    savefig(sprintf('Mean_V2_Ordi_cormat_%s.fig',mask_name))
    close

  %% BUBBLE VERSION
  clear ax
  clear fh
% Compute center of each circle 
% This assumes the x and y values were not entered in imagesc()
  x = 1 : 1 : size(tab,2); % x edges
  y = 1 : 1 : size(tab,1); % y edges
  [xAll, yAll] = meshgrid(x,y); 
% Set color of each rectangle
% Set color scale
  cmap = parula(256);
  cmap = jet;
  Cscaled = (tab - clrLim(1))/range(clrLim); % always [0:1]
  colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size in the same way we scale the color
  diamSize = Cscaled * range(diamLim) + diamLim(1); 
% diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling

% Create figure
            fh = figure(); 
            ax = axes(fh);
            title([t], 'Interpreter', 'none');
            hold(ax,'on')
            colormap(ax,'jet');
            % Create circles
            theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
            % If you get an error in this line below, "Index in position 1 is invalid",
            % that probably means you didn't set clrLim correctly. 
            h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
                diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
            axis(ax,'equal')
            axis(ax,'tight')
            set(ax,'YDir','Reverse')
            colorbar()
            caxis(clrLim);


            ax = gca;
            ax.XTickLabelRotation = 0;
            ax.YTickLabelRotation = 0;
            ax.XTickLabel = roilabels;
            ax.YTickLabel = roilabels;
            ax.TickLabelInterpreter = 'none';
            
            
            savefig(sprintf('Mean_V2_Ordi_cormat_%s_bubble.fig',mask_name))
            close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = sprintf('Mean correlation matrix %s V1 :',mask_name);
    roilabels = vertcat(roiname(:));
    for ir = 1: length(roilabels)
        roilabels{ir} = roilabels{ir}(5:end-12);
    end
    cd(main_dir)
    tab = csvread(sprintf('%s/RS_corel_matrix/means/V1_%s_mean.csv',main_dir,mask_name(1:end-5)),2,1);
    clrLim = [-1 1];
    diamLim = [0.3, 1];

    f = figure()
    imagesc(tab)
    title([t], 'Interpreter', 'none');
    colormap(gca,'jet');
    c = colorbar;
    caxis([-1 1]);
%            c.Label.String = "Correlation d'activité au repos entre ROI";
    axis equal
    axis tight

    ax = gca;
    ax.XTickLabelRotation = 90;
    ax.YTickLabelRotation = 0;
    ax.XTick = [1:1:length(roilabels)];
    ax.YTick = [1:1:length(roilabels)];
    ax.XTickLabel = roilabels;
    ax.YTickLabel = roilabels;
    ax.TickLabelInterpreter = 'none';


    savefig(sprintf('Mean_V1_cormat_%s.fig',mask_name))
    close
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET ALL THE ROIs in the same place : individual corr matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ROI_all = vertcat(fileROI{1},fileROI{2},fileROI{3},fileROI{4});
ROI_all = fileROI;
% Verif VOI
%% get VOIs for each run

nRun   = 2;
if length(V1) == length(V2)
    nSubj = length(V2);
else
    nSubj = lenght(input_dir);
end

nROI = length(ROI_all);

voi_file = cell(nSubj,nROI,nRun);

for subj = 1 : nSubj
    for ir = 1 : nROI

       [~,roi_name] = fileparts(ROI_all{ir});
       %roi_name = roi_name(1:end-12);
       

          for iRun = 1 : nRun
               voi_file{subj,ir,iRun} = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
%% if proper files exist in input dir
              voi_file(subj,ir,iRun) = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
          end

% %% otherwise
%         V1(subj).gser('model').addVolume(sprintf('VOI_%s_reslice25_1.mat',roi_name),roi_name,1);
%         voi_file(subj,ir,1) = V1(subj).gser('model').gvol(roi_name).toJob;
%         if nRun == 2
%             V2(subj).gser('model').addVolume(sprintf('VOI_%s_reslice25_1.mat',roi_name),roi_name,1);
%             voi_file(subj,ir,2) = V2(subj).gser('model').gvol(roi_name).toJob;
%         end    
    end
end

voi = voi_file;
roi_labels = {};

for ir = 1 : nROI
    [~,roi_name] = fileparts(voi{1,ir,1});
    roi_name = roi_name(5:end-12);
    
    roi_labels = vertcat(roi_labels,roi_name);
end
    
%fmat = get_subdir_regex_files(input_dir,ROIs);
%fmat{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for subj = 1 : nSubj
    for iRun = 1 : nRun
%             f1name = sprintf('/V%d_%s_corr_matrix_%s.csv',iRun,patient_list{subj}(12:end),mask_name);
%             f2name = sprintf('/V%d_%s_pval_corr_matrix_%s.csv',iRun,patient_list{subj}(12:end),mask_name);
        f1name = sprintf('/V%d_%s_corr_matrix_all.csv',iRun,patient_list{subj}(12:end));
        f2name = sprintf('/V%d_%s_pval_corr_matrix_all.csv',iRun,patient_list{subj}(12:end));
        fout1 = addprefixtofilenames(input_dir{subj},f1name);
        fout2 = addprefixtofilenames(input_dir{subj},f2name);

        sujroi = voi(subj,:,iRun);
        time_course_matrix = {};
        parfor i = 1 : length(sujroi)
            l = load(sujroi{i});
            time_course_matrix{i} = l.Y;
        end
        time_course_matrix = cell2mat(time_course_matrix);

        [cor_mat, pval] = corr(time_course_matrix);

%% text output   
        

        Tcor = array2table(cor_mat, 'VariableNames', roi_labels);
        Tpval = array2table(pval, 'VariableNames', roi_labels);

        writetable(Tcor,fout1)
        writetable(Tpval,fout2)

%% visual output

        t = sprintf('Correlation matrix V%d :',iRun);
        roilabels = vertcat(roiname(:));
        for ir = 1: length(roilabels)
            roilabels{ir} = roilabels{ir}(5:end-12);
        end

        clrLim = [-1 1];
        diamLim = [0.3, 1];

        figure()
        imagesc(cor_mat)
        title([t, patient_list{subj}], 'Interpreter', 'none');
        colormap(gca,'jet');
        c = colorbar;
        caxis([-1 1]);
%            c.Label.String = "Correlation d'activité au repos entre ROI";
        axis equal
        axis tight

        ax = gca;
        ax.XTickLabelRotation = 90;
        ax.YTickLabelRotation = 0;
        ax.XTick = [1:1:length(roilabels)];
        ax.YTick = [1:1:length(roilabels)];
        ax.XTickLabel = roilabels;
        ax.YTickLabel = roilabels;
        ax.TickLabelInterpreter = 'none';


        savefig(sprintf('%sV%d_cormat_all.fig',input_dir{subj},iRun))
        close

%             %% BUBBLE VERSION
%                        
%             clear ax
%             clear fh
% % Compute center of each circle 
% % This assumes the x and y values were not entered in imagesc()
%             x = 1 : 1 : size(a,2); % x edges
%             y = 1 : 1 : size(a,1); % y edges
%             [xAll, yAll] = meshgrid(x,y); 
% % Set color of each rectangle
% % Set color scale
%             cmap = parula(256); 
%             Cscaled = (cor_mat - clrLim(1))/range(clrLim); % always [0:1]
%             colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% % Set size of each circle
% % Scale the size in the same way we scale the color
%             diamSize = Cscaled * range(diamLim) + diamLim(1); 
% % diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling
% 
% % Create figure
%             fh = figure(); 
%             ax = axes(fh);
%             title([t, patient_list{subj}], 'Interpreter', 'none');
%             hold(ax,'on')
%             colormap(ax,'parula');
%             % Create circles
%             theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
%             % If you get an error in this line below, "Index in position 1 is invalid",
%             % that probably means you didn't set clrLim correctly. 
%             h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
%                 diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
%             axis(ax,'equal')
%             axis(ax,'tight')
%             set(ax,'YDir','Reverse')
%             colorbar()
%             caxis(clrLim);
% 
% 
%             ax = gca;
%             ax.XTickLabelRotation = 0;
%             ax.YTickLabelRotation = 0;
%             ax.XTickLabel = roilabels;
%             ax.YTickLabel = roilabels;
%             ax.TickLabelInterpreter = 'none';
%             
%             
%             savefig(sprintf('%sV%d_cormat_bubble_%s.fig',input_dir{subj},iRun,mask_name))
%             close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform to Z-score & calculate all connectivity measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        corrMatrix = load(fout1);

        N1 = 1:length(cellstr(char(gfile(ROI_mask{1},'.*'))))
        N2 = 1:length(cellstr(char(gfile(ROI_mask{2},'.*'))))
        N3 = 1:length(cellstr(char(gfile(ROI_mask{3},'.*'))))
        N4 = 1:length(cellstr(char(gfile(ROI_mask{4},'.*'))))
        N5 = 1:length(cellstr(char(gfile(ROI_mask{5},'.*'))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:length(corrMatrix)

            tablemtx = readtable(corrMatrix{i});
            voiname = tablemtx.Properties.VariableNames
            mtx = table2array(tablemtx);
            a = threshold_absolute(mtx,-1);
            z = atanh(a)
            w = abs(z)

            dan = w(N1,N1)
            onw = w(N2,N2)
            rnw = w(N3,N3)
            smn = w(N4,N4)


            zdan(i,1) = mean_connection_strength(dan)
            zonw(i,1) = mean_connection_strength(onw)
            zrnw(i,1) = mean_connection_strength(rnw)
            zsmn(i,1) = mean_connection_strength(smn)

            zornw(i,1) = mean_connection_strength(ornw)

            % a = sum(sum(dan))/((8*(7)/2))
            % W = weight_conversion(w, 'binarize');

            dan_smn = w(N1,N4)
            onw_dan = w(N2,N1)
            rnw_dan = w(N3,N1)

            onw_smn = w(N2,N4)
            rnw_smn = w(N3,N4)
            onw_rnw = w(N2,N3)


            zdan_smn(i,1) =  mean_strength_inter(dan_smn)
            zonw_dan(i,1) =  mean_strength_inter(onw_dan)
            zrnw_dan(i,1) =  mean_strength_inter(rnw_dan)
            zonw_smn(i,1) =  mean_strength_inter(onw_smn)
            zrnw_smn(i,1) =  mean_strength_inter(rnw_smn)


            sw(i) = small_worldness(W,100)    

            Integrity(:,i) = (strengths_und(w))'; % C'est la somme des poids des liens dans un noeud 


        %     modularity(:,i) = modularity_und(W); % c'est le nombre de groupe ou de communauté
        %     Eglobal(i) = efficiency_bin(W);   
        %     Elocal(:,i)=  efficiency_bin(W,1);
        %     clustering(:,i) = clustering_coef_bu(W);


        end
%%% A PART TO BE ADAPTED TO MY DATA
        cout = struct;
        cout = setfield(cout,'suj',sujname)
        cout = setfield(cout,'pool',pool)
        %cout = setfield(cout,'voiname',voiname)
        cout = setfield(cout,'Z_DAN',zdan)
        cout = setfield(cout,'Z_ONW',zonw)
        cout = setfield(cout,'Z_RNW',zrnw)
        cout = setfield(cout,'Z_SMN',zsmn)
        cout = setfield(cout,'Z_ORNW',zornw)


        cout = setfield(cout,'Z_DAN_SMN',zdan_smn)
        cout = setfield(cout,'Z_ONW_DAN',zonw_dan)
        cout = setfield(cout,'Z_RNW_DAN',zrnw_dan)

        cout = setfield(cout,'Z_ONW_SMN', zonw_smn)
        cout = setfield(cout,'Z_RNW_SMN', zrnw_smn)
        cout = setfield(cout,'Z_ORNW_SMN',zornw_smn )
        cout = setfield(cout,'Z_ORNW_DAN',zornw_dan )

        write_result_to_csv(cout,'inter_nw_mesures_natif.csv')


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        fields = fieldnames(cout)
        fields_data = fields(3:end)

        for j =1:length(fields_data)
            myfield = fields_data{j};

            data = getfield(cout,myfield)

            pp = do_statistique_test(data(1:length(patients)), data(length(patients)+1:end));
            pvalue(j)=pp.p


            pval(j)= permutationTest(data(1:length(patients)), data(length(patients)+1:end), 200) % Patients de 1:24 et contrôles 25:48

                %pause
               % statGraphe = setfield(statGraphe,myfield, pval')    
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for nbr =1:(96)


            pp = do_statistique_test(Integrity(nbr,1:length(patients)), Integrity(nbr,length(patients)+1:end))
            zpvalue(nbr)=pp.p


            zpval(nbr)= permutationTest(Integrity(nbr,1:length(patients)), Integrity(nbr,length(patients)+1:end), 200); % Patients de 1:24 et contrôles 25:48


        end

        name = strrep(voiname,'VOI_VOI_w','')
        name = strrep(voiname,'VOI_VOI_','')
        name = strrep(name,'VOI_w','')
        name = strrep(name,'_1','')

        tab = array2table(Integrity','VariableNames',name')
        po = table(pool)
        a = zpval <0.05
        b = zpvalue <0.05
        an = voiname(a')      %% {'VOI_rcLeft_VIIIa_1'}    {'VOI_rcRight_VIIIa_1'}
        bn = voiname(b')

        tab = [po tab]
        writetable(tab,'integrity_nw1_natif.csv')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%


        resultat = readtable('resultat_reseau1_glm6.csv');
        voi = resultat.Roi
        name = strrep(voi,'VOI_ra','')
        name = strrep(name,'VOI_rs','')
        name = strrep(name,'VOI_rt','')
        name = strrep(name,'VOI_rh','')
        name = strrep(name,'VOI_rc','')
        name = strrep(name,'VOI_rw','')
        name = strrep(name,'VOI_ry','')
        name = strrep(name,'_1','')


        metric_network = resultat.Properties.VariableNames
        metric_network = metric_network(2:end-1)


        for nbrV =  1:length(metric_network)
            data_metric = getfield(resultat,metric_network{nbrV})
            mesVoi = voi(data_metric<0.05)
            mesnoms =name((data_metric<0.05))

            metric = struct
            metric = setfield(metric,'suj',sujname)
            metric = setfield(metric,'pool',pool)
            for nbrvoi = 1: length(mesVoi)
                indexVoi = find(strcmp(voiname,mesVoi{nbrvoi}))

                if strcmp(metric_network{nbrV}, 'degree')
                    metric = setfield(metric,strrep(mesnoms{nbrvoi},'-','_'),deg(indexVoi,:)')
                else
                    data_metric2 = eval(metric_network{nbrV})
                     metric = setfield(metric,strrep(mesnoms{nbrvoi},'-','_'),data_metric2(indexVoi,:)')
                end
            end

           write_result_to_csv(metric,['GLM6_' metric_network{nbrV} '.csv']);
        end

    end
end

export_figs('png')
