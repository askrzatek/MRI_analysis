%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VBM analysis script
%% PARKGAME II
%% January 2022
%% A.SKRZATEK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% Initialisation
CLUSTER = 1;
par.pct = 0;

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti_test');

e_PARKGAME = exam(main_dir,'PARKGAME.*._[a,c]$');

%e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');

%e_PARKGAME = exam(main_dir,'PARKGAME.*._exclu$');
e = e_PARKGAME; % (3:length(e_PARKGAME)); % choose specific

e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^v.*p2.nii','s',1)


%% Run CAT12 segmentation
clear par

%anat segment
fanat = e.getSerie('anat').getVolume('^s').getPath; % anat file

% Retrocompatibility for SPM:Spatial:Segment options
par.GM        = [1 1 1 1]; % warped_space_Unmodulated(wp*) / warped_space_modulated(mwp*) / native_space(p*) / native_space_dartel_import(rp*)
par.WM        = [1 1 1 1];
par.CSF       = [1 0 1 1];
par.TPMC      = [0 0 0 0];

par.bias      = [1 1 0] ;  % native normalize dartel     [0 1]; % bias field / bias corrected image
par.warp      = [1 1]; % warp field native->template / warp field native<-template
par.label     = [1 1 0];
par.las       = [0 0 0];

par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;
par.doSurface = 0;
par.subfolder = 1; % all results in the same subfolder

par.display = 0;
par.pct     = 0;
par.redo    = 0;

if CLUSTER
    par.run     = 0;
    par.sge     = 1;
else
    par.run     = 1;
    par.sge     = 0;
end

j_segment = job_do_segmentCAT12(fanat,par);

e.addSerie('t1mpr_S256_0_8iso_p2$','mri','cat',1)
e.getSerie('cat').addVolume('^mwp1.*p2.nii','mwp1',1)
e.gser('anat_T1').addVolume('^cat_v.*.xml','tiv',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smooth with 8

clear par
par.redo         = 0;
if CLUSTER
    par.run      = 0;
    par.sge      = 1;
else
    par.run      = 1;
    par.sge      = 0;
end
par.auto_add_obj = 0;
par.smooth       = [8 8 8];
par.prefix       = 's';
par.jobname      = 'spm_smooth_8';

img = e.gser('cat').gvol('^mwp1').removeEmpty;

job_smooth(img, par);

% save smoothed outputs

e.gser('cat').addVolume('^smwp1','smwp1',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make symbolic links in VBM firstlevel directory (MANUAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET INSPIRATION FROM RS_secondlevel_ROI_par_group.m
%% IT IS THE BEST WAY TO AUTOMATISE THIS PART : we need to execute classical stat tests (paired t-test, 2 sample t-test, etc.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd (main_dir)
stat_dir = fullfile(main_dir,'VBM');
Stat = exam(stat_dir,'PARK');

%% symbolic link creation ?? all analysis directly within the original dirs or we have to create a firstlevel vbm dir with the right subject names
%% Create OUTPUT-INPUT STRUCTURE for firstlevel data

Input   = fullfile(stat_dir,'firstlevel_vbm');

patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_053_LM_a_V1'};

for ipatient = 1: length(patient_list_V1)
    mkdir(Input, patient_list_V1{ipatient});
end

patients_dir_V1 = get_subdir_regex(Input, patient_list_V1); %%outdirs
patients_dir = get_subdir_regex(Input, patient_list); %%outdirs


%% ONCE THE INPUT STRUCTURE EXISTS
% separate files in analysis groups : V1_all, PARK_c, PARK_a
VBM_a = exam(Input,'PARK.*a$');
VBM_c = exam(Input,'PARK.*[c,DD]$');
VBM_all = exam(Input,'PARK.*a$') + exam(Input,'PARK.*[c,DD]$');
V1_all = exam(Input,'PARK.*a$') + exam(Input,'PARK.*a_V1') + exam(Input,'PARK.*[c,DD]$') + exam(Input,'PARK.*c_V1');

Sessions    = {'V1','V2'}%,'V2-V1'};
VBM_a.addSerie('report','TIV',1);
VBM_a.addSerie('mri','mri',1);

VBM_c.addSerie('report','TIV',1);
VBM_c.addSerie('mri','mri',1);

VBM_all.addSerie('report','TIV',1);
VBM_all.addSerie('mri','mri',1);

for iS = 1 : length(Sessions)
%PARK_a
    
    VBM_a.getSerie('mri').addVolume(sprintf('^smwp1_V%d',iS),sprintf('smwp1_V%d',iS));
    VBM_a.getSerie('TIV').addVolume(sprintf('^cat_v_V%d',iS),sprintf('TIV_V%d',iS));    
    
%PARK_c
    
    VBM_c.getSerie('mri').addVolume(sprintf('^smwp1_V%d',iS),sprintf('smwp1_V%d',iS));
    VBM_c.getSerie('TIV').addVolume(sprintf('^cat_v_V%d',iS),sprintf('TIV_V%d',iS));    
    
%all V1 & V2

    VBM_all.getSerie('mri').addVolume(sprintf('^smwp1_V%d',iS),sprintf('smwp1_V%d',iS));
    VBM_all.getSerie('TIV').addVolume(sprintf('^cat_v_V%d',iS),sprintf('TIV_V%d',iS));
    
end

% FIRSTLEVEL for ALL PARKs add volumes for each ROI's contrast

V1_all.addSerie('report','TIV',1);
V1_all.gser('TIV').addVolume('cat_v_V1','TIV_V1',1)
V1_all.addSerie('mri','mri',1);
V1_all.gser('mri').addVolume('smwp1_V1','smwp1_V1',1)

%% Regressors definition V1
%% AGE
covars_V1{1} = [70 % 001_NB_a
74              % 002_BM
72              % 003_SM_c
64              % 007_SD
76              % 008_JR
68              % 023_LJ
79              % 025_CA
68              % 028_PC_c
72              % 033_DD
61              % 039_KM
71              % 040_RE ?
72              % 042_RS ?
75              % 043_PD
56              % 044_CK
62              % 046_HJ ?
57              % 047_BF
66              % 048_SB
66              % 052_HJ
59];            % 053_LM
         
%% GENDER
covars_V1{2} = [1  % 001_NB_a
1               % 002_BM
1               % 003_SM_c
1               % 007_SD
2               % 008_JR
2               % 023_LJ
1               % 025_CA
2               % 028_PC_c
2               % 033_DD
1               % 039_KM
1               % 040_RE ?
2               % 042_RS ?
2               % 043_PD
2               % 044_CK
1               % 046_HJ
1               % 047_BF
2               % 048_SB
2               % 052_HJ
1];             % 053_LM

%% CAT12 TIV estimation
%% V1
%freport = e.getSerie('anat_T1').getVolume('tiv').toJob;
%catreport_xml   = {freport{1:2},freport{5},freport{7:8},freport{11:13},freport{17},freport{19:22},freport{24},freport{26},freport{28:29},freport{32:33}};
catreport_xml   = V1_all.gser('anat_T1').gvol('tiv').toJob;

clear par
par.fname       = 'TIV_V1';
par.foutdir     = main_dir;
par.run         = 1;

job_cat_TIV_estimate(catreport_xml,par);
load TIV_V1.txt

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
load e

patient_regex_V1 = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c','PARKGAMEII_040.*RE.*_a$','PARKGAMEII_042.*RS.*_a$','PARKGAMEII_046.*HJ.*_c$','PARKGAMEII_053.*LM.*_a$'};
patient_list_V1 = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c','PARKGAMEII_040_RE_a_V1','PARKGAMEII_042_RS_a_V1','PARKGAMEII_046_HJ_c_V1','PARKGAMEII_053_LM_a_V1'};
patient_regex = {'PARKGAMEII.*NB.*_a','PARKGAMEII.*BM.*_a','PARKGAMEII.*SM.*_c','PARKGAMEII.*SD.*_a','PARKGAMEII.*JR.*_a','PARKGAMEII.*LJ.*_c','PARKGAMEII.*CA.*_a','PARKGAMEII.*PC.*_c','PARKGAMEII.*DD.*','PARKGAMEII.*KM.*_a','PARKGAMEII.*PD.*_a','PARKGAMEII.*CK.*_c','PARKGAMEII.*BF.*_c','PARKGAMEII.*SB.*_a','PARKGAMEII_052.*HJ.*_c'};
patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a','PARKGAMEII_052_HJ_c'};


input_dir = get_subdir_regex(fullfile(main_dir,'/firstlevel_RS'),patient_list);

ROIs = {'Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','Postcentral_L','Postcentral_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10','Caudate_aal_L','Caudate_aal_R','Thalamus_aal_L','Thalamus_aal_R','Putamen_aal_L','Putamen_aal_R','Pallidum_aal_L','Pallidum_aal_R','PPN_YEB_L','PPN_YEB_R'};
% ROIs = {'SMA_face_L','SMA_face_R','SMA_foot_L','SMA_foot_R','SMA_hand_L','SMA_hand_R'};
dirROI  = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/ROI_pariet_mot_premot_cereb_BG_PPN';
fileROI = cellstr(char(gfile(dirROI,ROIs))); %fileROI = remove_regex(fileROI,'T1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verif VOI
%% get VOIs for each run

nRun   = 2;
nSubj = length(input_dir);
nROI  = length(ROIs)+1;

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
