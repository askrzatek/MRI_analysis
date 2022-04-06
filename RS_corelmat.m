%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS Correlation Matrix script
%% PARKGAME II
%% January 2022
%% A.SKRZATEK
%% inspired from S.OUARAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
patient_regex_V1 = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c','PARKGAMEII_040.*RE.*_a$','PARKGAMEII_042.*RS.*_a$','PARKGAMEII_046.*HJ.*_c$','PARKGAMEII_053.*LM.*_a$'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_053_LM_a_V1'};
patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'};
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};


input_dir = get_subdir_regex(fullfile(main_dir,'/firstlevel_RS'),patient_list);

% ROIs = {'Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','Postcentral_L','Postcentral_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10','Caudate_aal_L','Caudate_aal_R','Thalamus_aal_L','Thalamus_aal_R','Putamen_aal_L','Putamen_aal_R','Pallidum_aal_L','Pallidum_aal_R','PPN_YEB_L','PPN_YEB_R'};
ROIs = {'M1_face_L','M1_face_R','M1_foot_L','M1_foot_R','M1_hand_L','M1_hand_R'};

dirROI  = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_pariet_mot_premot_cereb_BG_PPN';
fileROI = cellstr(char(gfile(dirROI,ROIs))); %fileROI = remove_regex(fileROI,'T1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verif VOI
%% get VOIs for each run

nRun   = 2;
nSubj = length(input_dir);
%    nROI  = length(ROIs)+1;
nROI = length(voi_file)

voi_file = cell(nSubj,nROI,nRun);

for subj = 1 : nSubj
    for ir = 1 : nROI-1

       [~,roi_name] = fileparts(fileROI{ir});
       roi_name = roi_name(1:end-10);

       for iRun = 1 : nRun
           %voi_file{subj,ir,iRun} = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
           voi_file(subj,ir,iRun) = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
       end
    end
    for iRun = 1 : nRun
       %voi_file{subj,nROI,iRun} = get_subdir_regex_files(input_dir{subj},sprintf('VOI_PCC_%d.mat',iRun));
       voi_file(subj,nROI,iRun) = get_subdir_regex_files(input_dir{subj},sprintf('VOI_PCC_%d.mat',iRun));
    end   
end
voi = voi_file;
%fmat = get_subdir_regex_files(input_dir,ROIs);
%fmat{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for subj = 1 : nSubj
    for iRun = 1 : nRun
        f1name = sprintf('/V%d_corr_matrix.csv',iRun);
        f2name = sprintf('/V%d_pval_corr_matrix.csv',iRun);
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

        [~, roiname] = get_parent_path(sujroi);
        roiname = nettoie_dir(change_file_extension(roiname,''))

        Tcor = array2table(cor_mat, 'VariableNames', roiname);
        Tpval = array2table(pval, 'VariableNames', roiname);

        writetable(Tcor,fout1)
        writetable(Tpval,fout2)

%% visual output        
        t = sprintf('Correlation matrix V%d :',iRun);
        roilabels = vertcat(ROIs(:),{'PCC'});
        figure; imagesc(cor_mat)
        title([t, patient_list{subj}], 'Interpreter', 'none')
        ax = gca;
        x = 0.5:54.5;
        y = 0.5:54.5;
        ax.XTick = x;
        ax.YTick = y;
        ax.XTickLabelRotation = 90;
        ax.YTickLabelRotation = 0;
        ax.XTickLabel = roilabels;
        ax.YTickLabel = roilabels;
        savefig(sprintf('%sV%d_cormat.fig',input_dir{subj},iRun))
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
ROI_mask = get_subdir_regex(dirROI,'mask');

for imask = 1: length(ROI_mask)
    fileROI = cellstr(char(gfile(ROI_mask{imask},'.*')));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Verif VOI
    %% get VOIs for each run

    nRun   = 2;
    if length(V1) == length(V2)
        nSubj = length(V2);
    else
        nSubj = lenght(input_dir);
    end
    
    nROI = length(fileROI);

    voi_file = cell(nSubj,nROI,nRun);

    for subj = 1 : nSubj
        for ir = 1 : nROI

           [~,roi_name] = fileparts(fileROI{ir});
           roi_name = roi_name(1:end-10);

%           for iRun = 1 : nRun
               %voi_file{subj,ir,iRun} = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
%% if proper files exist in input dir
%               voi_file(subj,ir,iRun) = get_subdir_regex_files(input_dir{subj},sprintf('VOI_%s.*_%d.mat',roi_name,iRun));
%           end

%% otherwise
            V1(subj).gser('model').addVolume(sprintf('VOI_%s_reslice25_1.mat',roi_name),roi_name,1);
            voi_file(subj,ir,1) = V1(subj).gser('model').gvol(roi_name).toJob;
            if nRun == 2
                V2(subj).gser('model').addVolume(sprintf('VOI_%s_reslice25_1.mat',roi_name),roi_name,1);
                voi_file(subj,ir,2) = V2(subj).gser('model').gvol(roi_name).toJob;
            end    
        end
    end
    
    voi = voi_file;
    mask_path = strsplit(ROI_mask{imask},'/');
    mask_name = mask_path{end-1};
    
    %fmat = get_subdir_regex_files(input_dir,ROIs);
    %fmat{1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    for subj = 1 : nSubj
        for iRun = 1 : nRun
            f1name = sprintf('/V%d_%s_corr_matrix_%s.csv',iRun,patient_list{subj}(12:end),mask_name);
            f2name = sprintf('/V%d_%s_pval_corr_matrix_%s.csv',iRun,patient_list{subj}(12:end),mask_name);
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

            [~, roiname] = get_parent_path(sujroi);
            roiname = nettoie_dir(change_file_extension(roiname,''))
            
            
            Tcor = array2table(cor_mat, 'VariableNames', roiname);
            Tpval = array2table(pval, 'VariableNames', roiname);

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
            
            
            savefig(sprintf('%sV%d_cormat_%s.fig',input_dir{subj},iRun,mask_name))
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

        end
        %% Could be a good place to get V1 and V2 table values and make a difference (individual deltas)    
    end
    sujroi = voi(1,:,1);
    [~, roiname] = get_parent_path(sujroi);
    roiname = nettoie_dir(change_file_extension(roiname,''))
            
    
    t = sprintf('Mean delta correlation matrix %s KINECT:',mask_name);
            roilabels = vertcat(roiname(:));
            for ir = 1: length(roilabels)
                roilabels{ir} = roilabels{ir}(5:end-12);
            end
            cd(main_dir)
            tab = csvread(sprintf('%s/RS_corel_matrix/means/Kinect_%s_delta_mean_relative.csv',main_dir,mask_name(1:end-5)),1,1);
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
            
          
            savefig(sprintf('Mean_rel_Kinect_cormat_%s.fig',mask_name))
            close
            
%             %% BUBBLE VERSION
%             clear ax
%             clear fh
% % Compute center of each circle 
% % This assumes the x and y values were not entered in imagesc()
%             x = 1 : 1 : size(tab,2); % x edges
%             y = 1 : 1 : size(tab,1); % y edges
%             [xAll, yAll] = meshgrid(x,y); 
% % Set color of each rectangle
% % Set color scale
%             cmap = parula(256); 
%             Cscaled = (tab - clrLim(1))/range(clrLim); % always [0:1]
%             colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% % Set size of each circle
% % Scale the size in the same way we scale the color
%             diamSize = Cscaled * range(diamLim) + diamLim(1); 
% % diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling
% 
% % Create figure
%             fh = figure(); 
%             ax = axes(fh);
%             title([t], 'Interpreter', 'none');
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
%             savefig(sprintf('Mean_rel_Kinect_cormat_%s_bubble.fig.fig',mask_name))
%             close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t = sprintf('Mean delta correlation matrix %s ORDI :',mask_name);
            roilabels = vertcat(roiname(:));
            for ir = 1: length(roilabels)
                roilabels{ir} = roilabels{ir}(5:end-12);
            end
            cd(main_dir)
            tab = csvread(sprintf('%s/RS_corel_matrix/means/Ordi_%s_delta_mean_relative.csv',main_dir,mask_name(1:end-5)),2,1);
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
            
          
            savefig(sprintf('Mean_rel_Ordi_cormat_%s.fig',mask_name))
            close
            
%             %% BUBBLE VERSION
%             clear ax
%             clear fh
% % Compute center of each circle 
% % This assumes the x and y values were not entered in imagesc()
%             x = 1 : 1 : size(tab,2); % x edges
%             y = 1 : 1 : size(tab,1); % y edges
%             [xAll, yAll] = meshgrid(x,y); 
% % Set color of each rectangle
% % Set color scale
%             cmap = parula(256); 
%             Cscaled = (tab - clrLim(1))/range(clrLim); % always [0:1]
%             colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% % Set size of each circle
% % Scale the size in the same way we scale the color
%             diamSize = Cscaled * range(diamLim) + diamLim(1); 
% % diamSize = abs(C) * range(diamLim) + diamLim(1); for [-1,1] scaling
% 
% % Create figure
%             fh = figure(); 
%             ax = axes(fh);
%             title([t], 'Interpreter', 'none');
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
%             savefig(sprintf('Mean_rel_Ordi_cormat_%s_bubble.fig.fig',mask_name))
%             close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            
end

export_figs('png')
