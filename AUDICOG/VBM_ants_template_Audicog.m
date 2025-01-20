%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ants template creation + VBM analysis script : combined scripts from A.SKRZATEK & S.OUARAB
%% AUDICOG
%% June 2024
%% A.SKRZATEK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% Initialisation
CLUSTER = 1;
par.pct = 0;

main_dir = fullfile('/network/iss/cenir/analyse/irm/studies/AUDICOG','DATA/Non_chirurgicaux');
cd (main_dir)

suj = exam(main_dir,'20.*._AUDICOG_Sujet.*');
e = suj; % (3:length(e_PARKGAME)); % choose specific

e.addSerie('S.*T1w$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^v.*T1w.nii','s',1)


%% Run CAT12 segmentation

%anat segment
fanat = e.getSerie('anat').getVolume('^s').getPath; % anat file

clear par 

par.subfolder = 1;         % write in subfolder, just for better organisation (mri and report folders)

par.GM        = [0 0 1 2]; % index 3 = native_space (p1), value 1 = yes. // index 4 = native_space_dartel_import (rp1) value 2 = Affine
par.WM        = [0 0 1 2]; %                        (p2)                 //           rp2
par.CSF       = [0 0 1 0]; % just                  (p3)            
par.TPMC      = [0 0 0 0]; % just            (p4,p5,p6) don't need them for VBM

par.label     = [0 0 0] ;  % don't need label map 
par.bias      = [0 0 0] ;  % don't save the bias field corrected  + SANLM (global) T1
par.las       = [0 0 0] ;  % don't save the bias field corrected  + SANLM (local)  T1
par.warp      = [1 1]   ;  % warp field native to MNI (y_) / MNI to native(iy_) for normalize or de-normalize images 
par.redo      = 1;

if CLUSTER
    par.run          = 0;
    par.sge          = 1;   % using cluster 
    par.mem          = '16G';
    par.jobname      = 'VBM_SEG_CAT12';
else
    par.run     = 1;
    par.sge     = 0;
end

% addpath('/network/iss/cenir/analyse/irm/users/salim.ouarab/dev_matvol/iNprogress')
job_do_segmentCAT12(fanat,par)



e.addSerie('T1w$','mri','cat',1)
e.getSerie('cat').addVolume('^rp1','gm_dartel',1)
e.getSerie('cat').addVolume('^rp2','wm_dartel',1);
% frmg = gpath(e.getSerie('cat').getVolume('gm_dartel'));
% frmb = gpath(e.getSerie('cat').getVolume('wm_dartel'));

vbm = gpath(e.getSerie('cat'));
frmg = gfile(vbm,'^rp1');
frmb = gfile(vbm,'^rp2');

clear par
par.run          = 0;
par.sge          = 1;   % using cluster 
par.mem          = '16G';
par.jobname      = 'Dartel_Template';

job_do_dartel_template(frmg,frmb,par);

% e.getSerie('cat').addVolume('^mwp1.*p2.nii','mwp1',1)
% e.gser('anat_T1').addVolume('^cat_v.*.xml','tiv',1)

% generate u_ (flow field) for each suj and 6 templates in the first suj
% folder - to be moved to the mqin_dir manually or with command

%% Normalize MG 
% change folder if using CAT12 and take p1.*nii
e.getSerie('cat').addVolume('^u.*.nii','ffield',1)
e.getSerie('cat').addVolume('^p1.*.nii','gm',1)

template  = gfile(main_dir,'Template_6');
ffield    = gfile(vbm,'^u.*nii');
mg        = gfile(vbm,'^p1.*nii');

clear par
par.preserve = 1;         % modulation [mw ]
par.fwhm     = [8 8 8];   % smooth     [smw ] default kernel size
par.redo     = 1;
par.run          = 0;
par.sge          = 1;   % using cluster 
par.mem          = '16G';
par.jobname      = 'Dartel_Normalize_CAT12';

job_dartel_normalize(template,ffield,mg,par)


%%%% ===> Check normalized data with "chekreg"

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
Sessions = exam(Input,'^PARKGAMEII.*');
Sessions.addSerie('report','TIV',1);
Sessions.addSerie('mri','mri',1);

Sessions.gser('mri').addVolume('^swp1v_V1.*.nii$','smwp1',1);
Sessions.gser('TIV').addVolume('^cat_v.*.xml','TIV',1);

addpath('/home/anna.skrzatek/MRI_analysis/')

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VBM V2-V1 individual statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cocovars{1,:} = [covars{1,1};covars{1,2}];
cocovars{2,:} = [covars{2,1};covars{2,2}];
coTIV_reg = [TIV_reg{:,1};TIV_reg{:,2}];

for iS = 1 : length(patients_dir)

    cocovars_double{1,iS}(1) = cocovars{1}(iS);
    cocovars_double{2,iS}(1) = cocovars{2}(iS);
    cocovars_double{3,iS}(1) = coTIV_reg(iS,1);
    
    cocovars_double{1,iS}(2) = cocovars{1}(iS);
    cocovars_double{2,iS}(2) = cocovars{2}(iS);
    cocovars_double{3,iS}(2) = coTIV_reg(iS,2);

end

for iS = 1:length(patients_dir)
    mkdir(patients_dir{iS},'V2_V1')
end

out_a = exam(Input,'^PARKGAME.*_a$');
out_c = exam(Input,'^PARKGAME.*[D,c]$');
out = out_a + out_c;
out.addSerie('V2_V1','V2_V1',1);
cooutdir = gpath(out.gser('V2_V1'));

subj = [VBM_a, VBM_c];

par.run = 1;

firstlevel_2s_vbm_matlabbatch(subj,cooutdir,cocovars_double,par)

%% Single paired t-test model estimate

for iout = 10 : length(cooutdir)
    
    fspm = addsuffixtofilenames(cooutdir(iout), '/SPM.mat');

    %fspm = addsuffixtofilenames(cooutdir(15), '/SPM.mat');

    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = 'spm_single_paired_secondlevel_vbm_est';
    job_first_level_estimate(fspm,par)
    
end

%% Contrast creation for each SPM.mat : single paired t-test
% t-statistics
V1 = [1 0];    
V2 = [0 1];
%%%%%% ???

for iout = 1 : length(cooutdir)
    
    fspm = addsuffixtofilenames(cooutdir(iout), '/SPM.mat');

%% Contrast names
    contrast_t.names = {
        'VBM V1'
        'VBM V2'
        'VBM V1<V2'
        'VBM V1>V2'
        }';

    %% Contrast values
    contrast_t.values = {
        V1
        V2
        V2-V1
        V1-V2
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
%% REGRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covarsV1_reg = {vertcat(covars{1,:}),vertcat(covars{2,:}),vertcat(TIV_reg{1,1},TIV_reg{1,2})};
covarsV2_reg = {vertcat(covars{1,:}),vertcat(covars{2,:}),vertcat(TIV_reg{2,1},TIV_reg{2,2})};

%% AGE
covars{1} = [70 % 001_NB_a
74              % 002_BM
64              % 007_SD
76              % 008_JR
79              % 025_CA
61              % 039_KM
71              % 040_RE ?
72              % 042_RS ?
75              % 043_PD
66              % 048_SB
59              % 053_LM ?
72              % 003_SM_c
68              % 023_LJ
68              % 028_PC_c
72              % 033_DD
56              % 044_CK
62              % 046_HJ ?
57              % 047_BF
66];            % 052_HJ
         
%% GENDER
covars{2} = [1  % 001_NB_a
1               % 002_BM
1               % 007_SD
2               % 008_JR
1               % 025_CA
1               % 039_KM
1               % 040_RE ?
2               % 042_RS ?
2               % 043_PD
2               % 048_SB
1               % 053_LM ?
1               % 003_SM_c
2               % 023_LJ
2               % 028_PC_c
2               % 033_DD
2               % 044_CK
1               % 046_HJ ?
1               % 047_BF
2];             % 052_HJ

%% TIV
covars{3} = vertcat(TIV_reg{1,1},TIV_reg{1,2});

%% Target regressors definition & their outdirs
% define or create output directory
mkdir(main_dir,'full_resliced_multiple_regression_V1')
MultiRegDir = fullfile(main_dir,'full_resliced_multiple_regression_V1')

if ACTION
    MultiRegDirC = fullfile(MultiRegDir, '/ACT_clinic_V1');
    Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII_AXIAL','UPDRSIII_SUP'};
    
    MultiRegDirG = fullfile(MultiRegDir, '/ACT_gait_V1/Spontaneous');
    Gait  = {'APA_AP','DA','Step_Size'};
    %Gait_up  = {'Rapid','Spontaneous'}; % the Rapid gait condition is one subject shorter - needs a separate processing
end

if RS
    MultiRegDirC = fullfile(MultiRegDir, '/RS_clinic_V1');
    Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII_AXIAL','UPDRSIII_SUP'};
    
    MultiRegDirG = fullfile(MultiRegDir, '/RS_gait/Spontaneous');
    Gait  = {'APA_AP','DA','Step_Size'};
    %Gait_up  = {'Rapid','Spontaneous'}; % the Rapid gait condition is one subject shorter - needs a separate processing
end

%RegDir = horzcat({RegDirG}, {RegDirC});
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

MultiStat = MultiStatObjG{1}{1} + MultiStatObjG{1}{2} + MultiStatObjG{1}{3} + MultiStatObjC{2}{1} + MultiStatObjC{2}{2} + MultiStatObjC{2}{3} + MultiStatObjC{2}{4} + MultiStatObjC{2}{5} ;
MultiStat.explore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V1 regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spontan APA_AP
targetsAPA1 = [34.7171835889    % 001_NB
20.8184436256                   % 002_BM
14.0398756669                   % 007_SD
26.6625447276                   % 008_JR
30.8778153637                   % 025_CA
34.7289858365                   % 039_KM
12.809081479101                 % 040_RE ?
38.1662614156685                % 042_RS ?
36.4238611957                   % 043_PD
30.6084803157                   % 048_SB
55.1500901593812                % 053_LM ?
37.3914711347                   % 003_SM_c
33.3035609785                   % 023_LJ
41.7655772122                   % 028_PC
26.2586959347                   % 033_DD
32.1145722308                   % 044_CK
34.2056328066679                % 046_HJ ?
35.176597489                    % 047_BF
31.2547254406906];              % 052_HJ

%% Spontan DA
targetsDA1 = [0.229615923       % 001_NB
0.2486013264                    % 002_BM
0.2995684834                    % 007_SD
0.2920857032                    % 008_JR
0.2572407702                    % 025_CA
0.2558116756                    % 039_KM
0.26562558007923                % 040_RE ?
0.197120826259196               % 042_RS ?
0.1843411258                    % 043_PD
0.2517826825                    % 048_SB
0.210950764006791               % 053_LM ?
0.1624504829                    % 003_SM_c
0.2447760219                    % 023_LJ
0.2448816538                    % 028_PC
0.2531833616                    % 033_DD
0.2184634975                    % 044_CK
0.283773951006549               % 046_HJ ?
0.1773839842                    % 047_BF
0.245088527771041];             % 052_HJ

%% Spontan Step_Size
targetsSS1 = [257.3421060198    % 001_NB
461.4812761193                  % 002_BM
249.806860919                   % 007_SD
388.3478308826                  % 008_JR
163.3475358607                  % 025_CA
275.7546680306                  % 039_KM
162.962895971799                % 040_RE ?
415.737032074417                % 042_RS ?
315.6876356933                  % 043_PD
371.2412257681                  % 048_SB
406.869666413503                % 053_LM ?
230.2511033831                  % 003_SM_c
302.5132704199                  % 023_LJ
262.8219167263                  % 028_PC
340.6669320703                  % 033_DD
421.525704218                   % 044_CK
336.690629708893                % 046_HJ ?
447.2221013178                  % 047_BF
240.696518058572];              % 052_HJ

%% Rapid APA_AP
targetrAPA1 = [56.1860003179    % 001_NB
42.202244009                    % 002_BM
19.9279281562                   % 007_SD
38.3572741643                   % 008_JR
43.1886233663                   % 039_KM
% 040_RE ?
% 042_RS ?
63.3029794047                   % 043_PD
57.0992020159                   % 048_SB
% 053_LM
61.9964561812                   % 003_SM_c
42.2617161059                   % 023_LJ
53.1842418093                   % 028_PC
51.3589704894                   % 033_DD
61.1835358307                   % 044_CK
% 046_HJ
51.1603131288                    % 047_BF
];  % 052_HJ

%% Rapid DA
targetrDA1 = [0.1781937653      % 001_NB
0.1892911011                    % 002_BM
0.2379586858                    % 007_SD
0.2200210944                    % 008_JR
0.1808432371                    % 039_KM
% RE ?
% RS ?
0.1457373272                    % 043_PD
0.2024911505                    % 048_SB
0.1217863129                    % 003_SM_c
0.2044142615                    % 023_LJ


0.1992642898                    % 028_PC
0.1682406621                    % 033_DD
0.1644736842                    % 044_CK
% HJ
0.1373089983                    % 047_BF
];  % 052_HJ
        
%% Rapid Step_Size
targetrSS1 = [385.33658393      % 001_NB
400.9512381886                  % 002_BM
326.7006392201                  % 007_SD
428.223371688                   % 008_JR
398.4146093581                  % 039_KM
% RE ?
% RS ?
502.8718930147                  % 043_PD
484.1753627556                  % 048_SB
200.2581921664                  % 003_SM_c
424.4445720736                  % 023_LJ


407.5565252707                  % 028_PC
513.2893254134                  % 033_DD
507.5458941188                  % 044_CK
% HJ
517.2496431319                    % 047_BF
];  % 052_HJ

%% AXIAL
targetAxial1 = [2               % 001_NB
4                               % 002_BM
7                               % 007_SD
7                               % 008_JR
10                              % 025_CA
4                               % 039_KM
5                               % 040_RE ?
14                              % 042_RS ?
1                               % 043_PD
3                               % 048_SB
1                               % 053_LM ?
3                               % 003_SM_c
4                               % 023_LJ
8                               % 028_PC
7                               % 033_DD
3                               % 044_CK
7                               % 046_HJ ?
2                               % 047_BF
11];                            % 052_HJ

%% GABS
targetGabs1 =  [20              % 001_NB
28                              % 002_BM
31                              % 007_SD
33                              % 008_JR
43                              % 025_CA
28                              % 039_KM
47                              % 040_RE ?
22                              % 042_RS ?
19                              % 043_PD
29                              % 048_SB
31                              % 053_LM ?
19                              % 003_SM_c
18                              % 023_LJ
32                              % 028_PC
26                              % 033_DD
14                              % 044_CK
21                              % 046_HJ ?
18                              % 047_BF
39];                            % 052_HJ

%% UPDRS III
targetUPDRS1 = [17              % 001_NB
22                              % 002_BM
28                              % 007_SD
25                              % 008_JR
36                              % 025_CA
13                              % 039_KM
18                              % 040_RE ?
23                              % 042_RS ?
20                              % 043_PD
15                              % 048_SB
8                               % 053_LM ?
23                              % 003_SM_c
33                              % 023_LJ
36                              % 028_PC
50                              % 033_DD
20                              % 044_CK
37                              % 046_HJ ?
31                              % 047_BF
34];                            % 052_HJ

%% UPDRSIII-AXIAL
targetUPDRSIII_Axial1 = [15     % 001_NB
18                              % 002_BM
21                              % 007_SD
18                              % 008_JR
26                              % 025_CA
9                               % 039_KM
13                              % 040_RE ?
9                               % 042_RS ?
19                              % 043_PD
12                              % 048_SB
7                               % 053_LM ?
20                              % 003_SM_c
29                              % 023_LJ
28                              % 028_PC
43                              % 033_DD
17                              % 044_CK
30                              % 046_HJ ?
29                              % 047_BF
23];                            % 052_HJ

%% UPDRSIII-MEMBRES SUP
targetUPDRSIII_Sup1 = [4       % 001_NB
7                              % 002_BM
12                             % 007_SD
6                              % 008_JR
5                              % 025_CA
13                             % 039_KM
7                              % 040_RE ?
4                              % 042_RS ?
12                             % 043_PD
6                              % 048_SB
5                              % 053_LM ?
14                             % 003_SM_c
2                              % 023_LJ
5                              % 028_PC
5                              % 033_DD
8                              % 044_CK
16                             % 046_HJ ?
4                              % 047_BF
11];                           % 052_HJ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
           
%% Creating groups structure
MultiStat.mkdir('V1');
MultiStat.addSerie('V1','V1');
MultiStat.mkdir('V2_V1');
MultiStat.addSerie('V2_V1','V2_V1');
multioutdirs = MultiStat.getSerie('.*') .toJob;

%% Multiregression V1 & V2-V1 all in    
%     VBM_a = {VBM_a_V1, VBM_a_V2};
%     VBM_c = {VBM_c_V1, VBM_c_V2};
%     covarsV2_V1_reg = vertcat(covarsV2_reg{1}, covarsV2_reg{2},covarsV2_reg{3}-covarsV1_reg{3});
%     covars_multi = {covarsV1_reg ;covarsV2_V1_reg};
%     
%     target_regressors.name  = {MultiStat.name};
%     
%     APA         = {targetsAPA1; targetsAPA};
%     DA          = {targetsDA1; targetsDA};
%     Axial       = {targetAxial1; targetAxial};
%     GABS        = {targetGabs1; targetGabs};
%     UPDRSIII    = {targetUPDRS1, targetUPDRS};
%
    VBM_a = {VBM_a_V1};
    VBM_c = {VBM_c_V1};
    covars_multi = covars;
    
    APA = {targetsAPA1};
    target_regressors.name = {MultiStat.name};
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

