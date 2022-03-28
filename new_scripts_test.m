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
Sessions = exam(Input,'^V[1,2]');
Sessions.addSerie('report','TIV',1);
Sessions.addSerie('mri','mri',1);

Sessions.gser('mri').addVolume('^smwp1.*a.nii$','smwp1_a',8);
Sessions.gser('mri').addVolume('^smwp1.*c.nii$','smwp1_c',7);

Sessions.gser('TIV').addVolume('^cat_rmv.*a.xml','TIV_a',8);
Sessions.gser('TIV').addVolume('^cat_rmv.*c.xml','TIV_c',7);


%% CAT12 TIV estimation
%% V1
%freport = e.getSerie('anat_T1').getVolume('tiv').toJob;
%catreport_xml   = {freport{1:2},freport{5},freport{7:8},freport{11:13},freport{17},freport{19:22},freport{24},freport{26},freport{28:29},freport{32:33}};
catreport_xml   = cellstr(char(Sessions.gser('TIV').gvol('TIV').toJob));

clear par
par.fname       = 'TIV_V1';
par.foutdir     = main_dir;
par.run         = 1;

job_cat_TIV_estimate(catreport_xml,par);
load TIV_V1.txt

%% V1 per group
catreport_xml_V1_a     = cellstr(char(Sessions(1).gser('TIV').gvol('TIV_a') .toJob));
catreport_xml_V1_c     = cellstr(char(Sessions(1).gser('TIV').gvol('TIV_c').toJob));

clear par
par.fname       = 'TIV_V1_a';
par.foutdir     = main_dir;
par.run         = 1;
job_cat_TIV_estimate(catreport_xml_V1_a,par);

clear par
par.fname       = 'TIV_V1_c';
par.foutdir     = main_dir;
par.run         = 1;
job_cat_TIV_estimate(catreport_xml_V1_c,par);

%% V2 per group
catreport_xml_V2_a     = cellstr(char(Sessions(2).gser('TIV').gvol('TIV_a') .toJob));
catreport_xml_V2_c     = cellstr(char(Sessions(2).gser('TIV').gvol('TIV_c').toJob));

clear par
par.fname       = 'TIV_V2_a';
par.foutdir     = main_dir;
par.run         = 1;
job_cat_TIV_estimate(catreport_xml_V2_a,par);

clear par
par.fname       = 'TIV_V2_c';
par.foutdir     = main_dir;
par.run         = 1;
job_cat_TIV_estimate(catreport_xml_V2_c,par);

TIV_reg{1,1} = load('TIV_V1_a.txt');
TIV_reg{1,2} = load('TIV_V1_c.txt');
TIV_reg{2,1} = load('TIV_V2_a.txt');
TIV_reg{2,2} = load('TIV_V2_c.txt');

%% Regressors definition
%% AGE
covars{1,1} = [70 % 001_NB_a
74              % 002_BM
64              % 007_SD
76              % 008_JR
79              % 025_CA
61              % 039_KM
75              % 043_PD
66              % 048_SB
];
covars{1,2} = [72              % 003_SM_c
68              % 023_LJ
68              % 028_PC_c
72              % 033_DD
56              % 044_CK
57              % 047_BF
66              % 052_HJ
];            
         
%% GENDER
covars{2,1} = [1  % 001_NB_a
1               % 002_BM
1               % 007_SD
2               % 008_JR
1               % 025_CA
1               % 039_KM
2               % 043_PD
2               % 048_SB
];
covars{2,2} = [1               % 003_SM_c
2               % 023_LJ
2               % 028_PC_c
2               % 033_DD
2               % 044_CK
1               % 047_BF
2               % 052_HJ
];             

%% Standard analysis paired t-test, 2 sample t-test and regression to come

VBM_a_V1 = cellstr(char(Sessions(1).gser('mri').gvol('smwp1_a').toJob));  
VBM_a_V2 = cellstr(char(Sessions(2).gser('mri').gvol('smwp1_a').toJob));  
VBM_c_V1 = cellstr(char(Sessions(1).gser('mri').gvol('smwp1_c').toJob));  
VBM_c_V2 = cellstr(char(Sessions(2).gser('mri').gvol('smwp1_c').toJob));  

for i = 1 : length(VBM_a_V1)
   VBM_a{i} = {VBM_a_V1{i}, VBM_a_V2{i}};
end
for i = 1 : length(VBM_c_V1)
   VBM_c{i} = {VBM_c_V1{i}, VBM_c_V2{i}};
end

anal_dir = fullfile(main_dir,'deltas');
outdirs = {cellstr(sprintf('%s/PARK_a_V1vsV2_paired_vbm',anal_dir)), cellstr(sprintf('%s/PARK_c_V1vsV2_paired_vbm',anal_dir))};
groups = {VBM_a, VBM_c};

for igroup = 1 : length(groups)
    i = 0;
    for icov = 1: length(covars{1,igroup})
        i = i+1;
        covars_double{1,igroup}(i) = covars{1,igroup}(icov);
        covars_double{2,igroup}(i) = covars{2,igroup}(icov);
        covars_double{3,igroup}(i) = TIV_reg{1,igroup}(icov);
        i = i+1;
        covars_double{1,igroup}(i) = covars{1,igroup}(icov);
        covars_double{2,igroup}(i) = covars{2,igroup}(icov);
        covars_double{3,igroup}(i) = TIV_reg{2,igroup}(icov);
    end
end

% define models per group
addpath /home/anna.skrzatek/MRI_analysis/

par.run = 1;
par.sge = 0;
par.jobname = 'job_VBM_paired_secondlevel_auto';

secondlevel_paired_vbm_matlabbatch(groups,outdirs,covars_double,par)

for iout = 1 : length(outdirs)
    
    fspm = addsuffixtofilenames(outdirs{iout}, '/SPM.mat');

    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_paired_secondlevel_vbm_est';
    job_first_level_estimate(fspm,par)
    
end

%% Contrast creation for each SPM.mat

% t-statistics
    V1_V2 = [ 1 -1];
    V2_V1 = [-1  1];

for iout = 1 : length(outdirs)
    
    fspm = addsuffixtofilenames(outdirs{iout}, '/SPM.mat');

%% Contrast names
    contrast_t.names = {
        'VBM V1>V2'
        'VBM V1<V2'
        }';

    %% Contrast values
    contrast_t.values = {
        V1_V2
        V2_V1
        }';

    %% Contrast type
    contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

    contrast.names  = [contrast_t.names];
    contrast.values = [contrast_t.values];
    contrast.types  = [contrast_t.types];

    %% Contrast : write
    clear par

    par.sge = 0;
    par.run = 1;
    par.display = 0;
    par.jobname = 'spm_write_vbm_con';

    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;
     
    job_first_level_contrast(fspm,contrast,par);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 SAMPLE T-TEST // shouldn't apply because of too small sample sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear groups
clear fspm
clear par

outdirs = {cellstr(sprintf('%s/PARK_a_vs_c_2sample_vbm',anal_dir))};
groups = {VBM_a_V2, VBM_c_V2};
covars_2s = {covars{1,:}; covars{2,:}; TIV_reg{2,1}, TIV_reg{2,2}};

par.run = 1;
par.sge = 0;
par.jobname = 'job_RS_2s_secondlevel_auto_covars';

secondlevel_vbm_2s_matlabbatch(groups,outdirs,covars_2s,par)

%% 2 sample model estimate

fspm = fullfile(outdirs{:},'SPM.mat');
par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.jobname  = 'spm_2s_secondlevel_vbm_est';
job_first_level_estimate(fspm,par)

%% 2 sample model write constrasts

% t-statistics
    Kinect_Ordi = [ 1 -1];
    Ordi_Kinect = [-1  1];
    
%% Contrast names
    contrast_t.names = {
        'VBM V2 K>O'
        'VBM V2 K<O'
        }';

    %% Contrast values
    contrast_t.values = {
        Kinect_Ordi
        Ordi_Kinect
        }';

    %% Contrast type
    contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

    contrast.names  = [contrast_t.names];
    contrast.values = [contrast_t.values];
    contrast.types  = [contrast_t.types];

    %% Contrast : write
    clear par

    par.sge = 0;
    par.run = 1;
    par.display = 0;
    par.jobname = 'spm_write_vbm_con';

    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;
     
    job_first_level_contrast(fspm,contrast,par);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REGRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covarsV1_reg = {vertcat(covars{1,:}),vertcat(covars{2,:}),vertcat(TIV_reg{1,1},TIV_reg{1,2})};
covarsV2_reg = {vertcat(covars{1,:}),vertcat(covars{2,:}),vertcat(TIV_reg{2,1},TIV_reg{2,2})};

mkdir(stat_dir,'full_resliced_multiple_regression')
MultiRegDir = fullfile(stat_dir,'multiple_regression')

MultiRegDirC = fullfile(MultiRegDir, '/VBM_clinic');
%   Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII_AXIAL','UPDRSIII_SUP'};
Clinic = {'AXIAL','GABS','UPDRSIII'};
    
MultiRegDirG = fullfile(MultiRegDir, '/VBM_gait/Spontaneous');
%   Gait  = {'APA_AP','DA','Step_Size'};
Gait  = {'APA_AP','DA'};

MultiRegDir = horzcat({MultiRegDirG}, {MultiRegDirC});
model = horzcat({Gait},{Clinic});
for imodel = 1:2
    targetmodel = model{imodel};
    for ivar = 1 : length(targetmodel)
       mkdir(MultiRegDir{imodel},sprintf('%s',targetmodel{ivar}));
       MultiStatObjG{imodel}{ivar} = exam(MultiRegDir{imodel}, sprintf('^%s',targetmodel{ivar})); % the Spontaneous gait condition is one subject shorter - needs a separate processing
       MultiStatObjC{imodel}{ivar} = exam(MultiRegDir{imodel},sprintf('^%s$',targetmodel{ivar})); % we have two of them containing UPDRSIII & AXIAL
    end
end

MultiStat = MultiStatObjG{1}{1} + MultiStatObjG{1}{2} + MultiStatObjC{2}{1} + MultiStatObjC{2}{2} + MultiStatObjC{2}{3};
MultiStat.explore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V1 regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spontan APA_AP
targetsAPA1 = [34.7171835889
20.8184436256
14.0398756669
26.6625447276
30.8778153637
34.7289858365
36.4238611957
30.6084803157
37.3914711347
33.3035609785
41.7655772122
26.2586959347
32.1145722308
35.176597489
31.2547254406906];

%% Spontan DA
targetsDA1 = [0.229615923
0.2486013264
0.2995684834
0.2920857032
0.2572407702
0.2558116756
0.1843411258
0.2517826825
0.1624504829
0.2447760219
0.2448816538
0.2531833616
0.2184634975
0.1773839842
0.245088527771041];
        
%% AXIAL
targetAxial1 = [2
4
7
7
10
4
1
3
3
4
8
7
3
2
11];
           
%% GABS
targetGabs1 =  [20
28
31
33
43
28
19
29
19
18
32
26
14
18
39];   

%% UPDRS III
targetUPDRS1 = [17
22
28
25
36
13
20
15
23
33
36
50
20
31
34];        
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V2 - V1 regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spontan APA_AP
targetsAPA= [23.1555713831
            1.40458578
            2.0220748182
            -2.0984142453
            13.9773798135
            3.4770321678
            21.0407162019
            9.1372754791
            12.9359263434
            7.3001966719
            9.4773158702
            0.2727108575
            -0.3172168418
            9.8691469678
            -2.09655588609517];

%% Spontan DA
targetsDA = [-0.0552436383
            -0.0173607532
            -0.0288778442
            0.0046389714
            0.0163977306
            -0.023941678
            -0.0038788317
            -0.0083970006
            0.0099967642
            0.0302428424
            -0.0200398573
            0.0146434635
            0.0547832049
            0.0004295012
            0.0375939849624061];

%% AXIAL
targetAxial = [0
               -1
               0
               -1
               -6
               0
               -1
               -2
               0
               1
               -2
               2
               -2
               0
               1];
           
%% GABS
targetGabs =  [-4
               -7
               0
               -2
               1
               -11
               -12
               -8
               -4
               11
               1
               1
               -1
               -5
               -8];   

%% UPDRS III
targetUPDRS = [1
               -3
               3
               -3
               -7
               -4
               -13
               -5
               -5
               -1
               -6
               -3
               -6
               4
               9];        
           
%% Creating groups structure
MultiStat.mkdir('V1');
MultiStat.addSerie('V1','V1');
MultiStat.mkdir('V2_V1');
MultiStat.addSerie('V2_V1','V2_V1');
multioutdirs = MultiStat.getSerie('.*') .toJob;

%% Multiregression V1 & V2-V1 all in    
%    VBM_a = {VBM_a_V1, VBM_a_V2_V1};    
%    VBM_c = {VBM_c_V1, VBM_c_V2_V1};
    VBM_a = {VBM_a_V1, VBM_a_V2};
    VBM_c = {VBM_c_V1, VBM_c_V2};
    covarsV2_V1_reg = vertcat(covarsV2_reg{1}, covarsV2_reg{2},covarsV2_reg{3}-covarsV1_reg{3});
    covars_multi = {covarsV1_reg ;covarsV2_V1_reg};
    target_regressors.name  = {MultiStat.name};
    APA         = {targetsAPA1; targetsAPA};
    DA          = {targetsDA1; targetsDA};
    Axial       = {targetAxial1; targetAxial};
    GABS        = {targetGabs1; targetGabs};
    UPDRSIII    = {targetUPDRS1, targetUPDRS};
    
    target_regressors.value = {APA,DA, Axial,GABS,UPDRSIII};
    
    par.run = 1;
    par.sge = 0;
    par.jobname = 'VBM_multireg_model_spec';
    par.nb_cond = length(multioutdirs);
    par.nb_ses = 1; %length(multioutdirs{1});
    
    multiregression_vbm_model_spec(multioutdirs, VBM_a, VBM_c, covars_multi, target_regressors, par)

%% multiregression model estimate

for iout = 1 : length(MultiStat)
    fspm = fullfile(multioutdirs{iout}(1),'SPM.mat');
    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_multireg_vbm_est';
    job_first_level_estimate(fspm,par)
end

%% regression model write constrasts
% t-statistics
    PosCorrelation = [0 0 0 1] ;
    NegCorrelation = [0 0 0 -1] ;

for iout = 1 : length(MultiStat)    
    %% Contrast names
        contrast_T.names = {
            sprintf('Pos_regression_%s_on_V1',target_regressors.name{iout})
            sprintf('Neg regression_%s_on_V1',target_regressors.name{iout})
            }';

        %% Contrast values
        contrast_T.values = {
            PosCorrelation
            NegCorrelation
            }';

        %% Contrast type
        contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

        contrast.names  = [contrast_T.names];
        contrast.values = [contrast_T.values];
        contrast.types  = [contrast_T.types];

        %% Contrast : write
        clear par

        par.sge = 0;
        par.run = 1;
        par.display = 0;
%        par.jobname = sprintf('spm_write_%s_%s_con',target_regressors.name{iout},roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(fspm,contrast,par);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This part I did through the interface : flexible factorial design 
%% - there is a batch saved too but I can't get the contrasts right
%% the output is in VBM/factorial_flex 2 groups 2 timepoints

cd (main_dir)
flexoutdir = fullfile(stat_dir,'factorial_flex');
fspm = addsuffixtofilenames(flexoutdir,'/SPM.mat');

%% variables according to Glascher & Gitelman (2008)

n1 = 8; %(number of subjects in group 1) 
n2 = 7; %(number of subjects in group 2)
nc = 2; %(number of levels in condition factor) 
ng = 2; %(number of groups) 
MEc = [1:nc]-mean(1:nc); %(main effect of condition, here: [-1 0 1]) 
MEg = [1 -1]; %(main effect of group: Group 1 > Group 2)  

%%
% F-test

Group_main      = [ones(1,nc)/nc -ones(1,nc)/nc ones(1,n1)/n1 -ones(1,n2)/n2];
Session_main    = [ MEc MEc zeros(1,n1+n2)];
SxG_inter       = [ MEc -MEc zeros(1,n1+n2)];

V1_K       = [ 1 0 0 0 ones(1,n1)/n1 zeros(1,n2)];
V2_K       = [ 0 1 0 0 ones(1,n1)/n1 zeros(1,n2)];
V1_O       = [ 0 0 1 0 ones(1,n1)/n1 zeros(1,n2)];
V2_O       = [ 0 0 0 1 ones(1,n1)/n1 zeros(1,n2)];

%% Contrast names
contrast_F.names = {
            'Main group effect'
            'Main session effect'
            'Session x Group interaction'
            'V1 effect Kinect'
            'V2 effect Kinect'
            'V1 effect Ordi'
            'V2 effect Ordi'
            }';

%% Contrast values
contrast_F.values = {
    Group_main
    Session_main
    SxG_inter
    V1_K
    V2_K
    V1_O
    V2_O
    }';

%% flex factorial classic instruction
%% t-test
%
% V1_V2_K     = [0 0 1 -1];  %Session V1 effect in group K
% V2_V1_K     = [0 0 -1 1];  %Session V2 effect in group K
% V1_V2_O     = [0 0 1 -1];  %Session V1 effect in group O
% V2_V1_O     = [0 0 -1 1];  %Session V2 effect in group O
% V1_V2_KO    = [1 -1 1 -1]; %Session V1 effect in groups mixed
% V2_V1_KO    = [-1 1 -1 1]; %Session V2 effect in groups mixed
% V1_V2_K_O   = [1 -1 -1 1]; %Interaction between 2 factors (group x session)
% 
% %% Contrast names
% contrast_t.names = {
%             'Session V1 > V2 Kinect'
%             'Session V1 < V2 Kinect'
%             'Session V1 > V2 Ordi'
%             'Session V1 < V2 Ordi'
%             'Session V1 > V2 K&O'
%             'Session V1 < V2 K&O'
%             'Session x Group interaction'
%             }';
%         
% %% Contrast values
% contrast_t.values = {
%     V1_V2_K
%     V2_V1_K
%     V1_V2_O
%     V2_V1_O
%     V1_V2_KO
%     V2_V1_KO
%     V1_V2_K_O
%     }';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% F-test
%
% Session_K       = [eye(2)-1/2 zeros(2)];
% Session_O       = [zeros(2) eye(2)-1/2];
% Session_main    = [eye(2)-1/2 eye(2)-1/2];
% SxG_inter       = [eye(2)-1/2 1/2-eye(2)];

% %% Contrast names
% contrast_F.names = {
%             'Session effect Kinect'
%             'Session effect Ordi'
%             'Main session effect'
%             'Session x Group interaction'
%             }';
% 
% %% Contrast values
% contrast_F.values = {
%     Session_K
%     Session_O
%     Session_main
%     SxG_inter
%     }';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contrast type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
contrast_t.types = cat(1,repmat({'T'},[1 length(contrast_t.names)]));

contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

%% Contrast mix
contrast.names  = [contrast_F.names contrast_t.names];
contrast.values = [contrast_F.values contrast_t.values];
contrast.types  = [contrast_F.types contrast_t.types];

%% Contrast : write
clear par

par.sge = 0;
par.run = 1;
par.display = 0;

% par.sessrep = 'both';
par.sessrep = 'none';

par.delete_previous = 1;
par.report          = 0;

job_first_level_contrast({fspm},contrast,par);

