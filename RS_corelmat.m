%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS Correlation Matrix script
%% PARKGAME II
%% January 2022
%% A.SKRZATEK
%% inspired from S.OUARAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fetch Input & Output dirs
main_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/';
cd (main_dir)

%% either existing e-file from RS processing
load e
%% or new e-object

V1 = exam(main_dir,'PARKGAME.*.V1_[a,c]$');
V1 = V1(1:11) + V1(13:14) + V1(16:18);
V1.addSerie('S0.*RS$','model$','model',1)
V2 = exam(main_dir,'PARKGAME.*.V2');


%%
patient_regex_V1 = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c','PARKGAMEII_040.*RE.*_a$','PARKGAMEII_042.*RS.*_a$','PARKGAMEII_046.*HJ.*_c$','PARKGAMEII_053.*LM.*_a$'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_053_LM_a_V1'};
patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'};
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};


input_dir = get_subdir_regex(fullfile(main_dir,'/firstlevel_RS'),patient_list);

% ROIs = {'Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','Postcentral_L','Postcentral_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10','Caudate_aal_L','Caudate_aal_R','Thalamus_aal_L','Thalamus_aal_R','Putamen_aal_L','Putamen_aal_R','Pallidum_aal_L','Pallidum_aal_R','PPN_YEB_L','PPN_YEB_R'};
ROIs = {'M1_face_L','M1_face_R','M1_foot_L','M1_foot_R','M1_hand_L','M1_hand_R'};

dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_pariet_mot_premot_cereb_BG_PPN';
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
main_dir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/';
cd (main_dir)
load e

patient_regex_V1 = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c','PARKGAMEII_040.*RE.*_a$','PARKGAMEII_042.*RS.*_a$','PARKGAMEII_046.*HJ.*_c$','PARKGAMEII_053.*LM.*_a$'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_053_LM_a_V1'};
patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'};
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};

input_dir = get_subdir_regex(fullfile(main_dir,'/firstlevel_RS'),patient_list);

dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROIs_masks';
ROI_mask = get_subdir_regex(dirROI,'mask');

for imask = 1: length(ROI_mask)
    fileROI = cellstr(char(gfile(ROI_mask{imask},'.*')));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Verif VOI
    %% get VOIs for each run

    nRun   = 2;
    nSubj = length(input_dir);
%    nROI  = length(ROIs)+1;
    nROI = length(fileROI)

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
            f1name = sprintf('/V%d_corr_matrix_.csv',iRun);
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
end
