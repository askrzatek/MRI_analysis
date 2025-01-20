%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script for creating the regression models between VBM & other statistical measures
%% ex. Clinical, Gait or Game Data regression for each fMRI condition
%% Anna SKRZATEK May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% Initialise
addpath /home/anna.skrzatek/MRI_analysis/

main_dir = fullfile('/network/iss/cenir/analyse/irm/users/anna.skrzatek','/nifti_test');
cd (main_dir)

% define input directory
InputVBM = fullfile(main_dir,'/VBM');


%% V1 REGRESSORS
%% AGE
covars_V1{1} = [70 % 001_NB_a
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
covars_V1{2} = [1  % 001_NB_a
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
%covars_V1{3} = vertcat(TIV_reg{1,1});
covars_V1{3} = [1398.22
1337.26
1344.44
1222.67
1270.20
1587.18
1284.63
1626.88
1594.15
1530.38
1329.31
1487.37
1426.57
1586.94
1587.62
1347.26
1657.54
1440.08
1268.36
];

%% Target regressors definition & their outdirs
% define or create output directory
mkdir(main_dir, '/full_VBM_clinic');
mkdir(main_dir, '/full_VBM_gait');
RegDirC = fullfile(main_dir, '/full_VBM_clinic');
RegDirG = fullfile(main_dir, '/full_VBM_gait');

MultiRegDir = fullfile(main_dir,'full_resliced_multiple_regression_delta');
mkdir(MultiRegDir, '/full_VBM_clinic')
mkdir(MultiRegDir, '/full_VBM_gait/Spontaneous')

MultiRegDirC = fullfile(MultiRegDir, '/full_VBM_clinic');
MultiRegDirG = fullfile(MultiRegDir, '/full_VBM_gait/Spontaneous');

RegDir = horzcat({RegDirG}, {RegDirC});
MultiRegDir = horzcat({MultiRegDirG}, {MultiRegDirC});

Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII_AXIAL','UPDRSIII_SUP'};
Gait  = {'APA_AP','DA','Step_Size'};

model = horzcat({Gait},{Clinic});

for imodel = 1:2
    targetmodel = model{imodel};
    for ivar = 1 : length(targetmodel)
       mkdir(RegDir{imodel},sprintf('%s',targetmodel{ivar})); % universal creation of directories depending on the chosen model
       mkdir(MultiRegDir{imodel},sprintf('%s',targetmodel{ivar}));
        % StatDir{ivar} = fullfile(RegDir,Gait{ivar});
        StatObjG{imodel}{ivar} = exam(RegDir{imodel}, sprintf('%s',targetmodel{ivar})); % the Spontaneous gait condition is one subject shorter - needs a separate processing
        StatObjC{imodel}{ivar} = exam(RegDir{imodel},sprintf('^%s$',targetmodel{ivar})); % we have two of them containing UPDRSIII & AXIAL
        MultiStatObjG{imodel}{ivar} = exam(MultiRegDir{imodel}, sprintf('^%s',targetmodel{ivar})); % the Spontaneous gait condition is one subject shorter - needs a separate processing
        MultiStatObjC{imodel}{ivar} = exam(MultiRegDir{imodel},sprintf('^%s$',targetmodel{ivar})); % we have two of them containing UPDRSIII & AXIAL
    end
end

Stat = StatObjG{1}{1} + StatObjG{1}{2} + StatObjG{1}{3} + StatObjC{2}{1} + StatObjC{2}{2} + StatObjC{2}{3} + StatObjC{2}{4} + StatObjC{2}{5};
Stat.explore

MultiStat = MultiStatObjG{1}{1} + MultiStatObjG{1}{2} + MultiStatObjG{1}{3} + MultiStatObjC{2}{1} + MultiStatObjC{2}{2} + MultiStatObjC{2}{3} + MultiStatObjC{2}{4} + MultiStatObjC{2}{5};
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
12.8090818479101 % 040_RE
38.1668614156625 % 042_RS
36.4238611957
30.6084803157
55.1500901593812 % 053_LM
37.3914711347
33.3035609785
41.7655772122
26.2586959347
32.1145722308
34.2056328066679 % 046_HJ
35.176597489
31.2547254406906];

%% Spontan DA
targetsDA1 = [0.229615923
0.2486013264
0.2995684834
0.2920857032
0.2572407702
0.2558116756
0.266562558007923 % 040_RE
0.197120826259196 % 042_RS
0.1843411258
0.2517826825
0.21095076400679 % 053_LM
0.1624504829
0.2447760219
0.2448816538
0.2531833616
0.2184634975
0.283773951006549 % 046_HJ
0.1773839842
0.245088527771041];
        
%% Spontan Step_Length
targetsSS1 = [257.3421060198
461.4812761193
249.806860919
388.3478308826
163.3475358607
275.7546680306
162.962895971799 % 040_RE
415.737032074417 % 042_RS
315.6876356933
371.2412257681
406.869666413503 % 053_LM
230.2511033831
302.5132704199
262.8219167263
340.6669320703
421.525704218
336.690629708893 % 046_HJ
447.2221013178
240.696518058572];        

%% AXIAL
targetAxial1 = [2
4
7
7
10
4
5 % 040_RE
14 % 042_RS
1
3
1 % 053_LM
3
4
8
7
3
7 % 046_HJ
2
11];
           
%% GABS
targetGabs1 =  [20
28
31
33
43
28
47 % 040_RE
22 % 042_RS
19
29
31 % 053_LM
19
18
32
26
14
21 % 046_HJ
18
39];   

%% UPDRS III
targetUPDRS1 = [17
22
28
25
36
13
20 % 040_RE
15 % 042_RS
8
18
23 % 053_LM
23
33
36
50
20
37 % 046_HJ
31
34];        
           
%% UPDRSIII-AXIAL
targetUPDRSIII_Axial1 = [15
18
21
18
26
9
13 % 040_RE
9 % 042_RS
19
12
7 % 053_LM
20
29
28
43
17
30 % 046_HJ
29
23];

%% UPDRSIII-MEMBRES SUP
targetUPDRSIII_Sup1 = [4
7
6
5
12
2
7 % 040_RE
4 % 042_RS
5
4
4 % 053_LM
12
13
6
14
5
16 % 046_HJ
8
11];
          
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

targetsAPAK= [23.1555713831
            1.40458578
            2.0220748182
            -2.0984142453
            13.9773798135
            3.4770321678
            21.0407162019
            9.1372754791];
        
targetsAPAO= [12.9359263434
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

targetsDAK = [-0.0552436383
            -0.0173607532
            -0.0288778442
            0.0046389714
            0.0163977306
            -0.023941678
            -0.0038788317
            -0.0083970006];
        
targetsDAO = [0.0099967642
            0.0302428424
            -0.0200398573
            0.0146434635
            0.0547832049
            0.0004295012
            0.0375939849624061];        
        
%% Spontan Step_Size
targetsSS = [150.103370391
            -123.5867048881
            -9.7077513532
            4.697752465
            45.0325100751
            8.0160415648
            131.2529377267
            -70.3945812668
            60.0510168179
            32.1219007707
            88.6692758751
            4.4760244141
            -117.5319247696
            -25.6892601563
            -25.0872759090923];
        
targetsSSK = [150.103370391
            -123.5867048881
            -9.7077513532
            4.697752465
            45.0325100751
            8.0160415648
            131.2529377267
            -70.3945812668];
        
targetsSSO = [60.0510168179
            32.1219007707
            88.6692758751
            4.4760244141
            -117.5319247696
            -25.6892601563
            -25.0872759090923];        

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
           
targetAxialK = [0
               -1
               0
               -1
               -6
               0
               -1
               -2];
           
targetAxialO = [0
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
           
targetGabsK =  [-4
               -7
               0
               -2
               1
               -11
               -12
               -8];
           
targetGabsO =  [-4
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
           
targetUPDRSK = [1
               -3
               3
               -3
               -7
               -4
               -13
               -5];
           
targetUPDRSO = [-5
               -1
               -6
               -3
               -6
               4
               9];
           
%% UPDRSIII-AXIAL
targetUPDRSIII_Axial = [1
                        -2
                        3
                        -2
                        -1
                        -4
                        -12
                        -3
                        -5
                        -2
                        -4
                        -5
                        -4
                        4
                        8];
                    
targetUPDRSIII_AxialK = [1
                        -2
                        3
                        -2
                        -1
                        -4
                        -12
                        -3];
                    
targetUPDRSIII_AxialO = [-5
                        -2
                        -4
                        -5
                        -4
                        4
                        8];
          
%% UPDRSIII-MEMBRES SUP
targetUPDRSIII_Sup = [2
                      0
                      2
                      -1
                      0
                      -1
                      -3
                      -1
                      -5
                      4
                      -3
                      0
                      -3
                       1
                      0];
                  
targetUPDRSIII_SupK = [2
                      0
                      2
                      -1
                      0
                      -1
                      -3
                      -1];
                  
targetUPDRSIII_SupO = [-5
                      4
                      -3
                      0
                      -3
                       1
                      0];
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VBM
%% Choice of VBM sessions //scans, outdir structure & names 
Sessions    = {'V1','delta_a','delta_c'}; %

% input different groups distinction and paths for the cons we need
% add volumes for each contrast - each individual

VBMObj1 = exam(InputVBM,'firstlevel'); % check directory names
VBMObj1.addSerie('^PARKGAME.*_a','V1_a',11);
VBMObj1.addSerie('^PARKGAME.*_[c,DD]','V1_c',8);
VBMObj1.gser('V1').addVolume('^smwp1.*nii$','vbm',1);


VBMObj = exam(InputVBM,'firstlevel'); % check directory names
% VBMObj.addSerie('^V1$','mri','V1',1);
% VBMObj.addSerie('^V2$','mri','V2',1);
VBMObj.addSerie('PARK.*a$','^V2_V1$','delta_a',8);
VBMObj.addSerie('PARK.*[c,DD]$','^V2_V1$','delta_c',7);

% VBMObj.gser('V1').addVolume('^smwp1.*PARK.*a.nii$','a',8);
% VBMObj.gser('V1').addVolume('^smwp1.*PARK.*[c,DD].nii$','c',7);
% 
% VBMObj.gser('V2').addVolume('^smwp1.*PARK.*a.nii$','a',8);
% VBMObj.gser('V2').addVolume('^smwp1.*PARK.*[c,DD].nii$','c',7);

VBMObj.gser('delta_a').addVolume('swp1.*','a');
VBMObj.getSerie('delta_c').addVolume('swp1.*','c');
VBMObj.addSerie('PARK.*a$','^report$','TIV_a',8);
VBMObj.addSerie('PARK.*[c,DD]$','^report$','TIV_c',7);
VBMObj.gser('TIV').addVolume('cat_v.*V1.xml','TIV_V1',1);
VBMObj.gser('TIV').addVolume('cat_v.*V2.xml','TIV_V2',1);

multicons_a = VBMObj1.getSerie('^V1_a').getVolume('vbm') .toJob;
multigcons_a = cellstr(multicons_a{:});
multicons_c =  VBMObj1.getSerie('^V1_c').getVolume('vbm') .toJob;
multigcons_c = cellstr(multicons_c{:});

cons_a = VBMObj.getSerie('delta_a').getVolume('a') .toJob;
gcons_a = cellstr(cons_a{:});
cons_c = VBMObj.getSerie('delta_c').getVolume('c') .toJob;
gcons_c = cellstr(cons_c{:});


Stat.mkdir('VBM_V2_V1_multireg');
Stat.addSerie('VBM_V2_V1_multi','VBM_V2_V1');
Stat.mkdir('VBM_V2_V1_multireg_a');
Stat.addSerie('VBM_V2_V1_multireg_a','VBM_V2_V1_a');
Stat.mkdir('VBM_V2_V1_multireg_c');
Stat.addSerie('VBM_V2_V1_multireg_c','VBM_V2_V1_c');

MultiStat.mkdir('VBM_V1');
MultiStat.addSerie('VBM_V1','VBM_V1');

outdirs = Stat.getSerie('.*') .toJob; % 7 x 1 x 12 cells
outdirs_a = Stat.getSerie('.*_a') .toJob;
outdirs_c = Stat.getSerie('.*_c') .toJob;
multioutdirs = MultiStat.getSerie('.*') .toJob;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regressors definition
%% AGE
covars{1} = [70 % 001_NB
74 % 002_BM
64 % 007_SD
76 % 008_JR
79 % 025_CA
61 % 039_KM
75 % 043_PD
66 % 048_SB
72 % 003_SM_c
68 % 023_LJ
68 % 028_PC
72 % 033_DD
56 % 044_CK
57 % 047_BF
66]; %052_HJ
         

covars_a{1} = [70 % 001_NB
74 % 002_BM
64 % 007_SD
76 % 008_JR
79 % 025_CA
61 % 039_KM
75 % 043_PD
66]; % 048_SB

covars_c{1} = [72 % 003_SM_c
68 % 023_LJ
68 % 028_PC
72 % 033_DD
56 % 044_CK
57 % 047_BF
66]; %052_HJ

%% GENDER
covars{2} = [1
1
1
2
1
1
2
2
1
2
2
2
2
1
2];

covars_a{2} = [1
1
1
2
1
1
2
2];

covars_c{2} = [1
2
2
2
2
1
2];

%% TIV
%% V1 per group

catreport_xml_V1     = cellstr(char(VBMObj.gser('TIV').gvol('TIV_V1') .toJob));
catreport_xml_V2     = cellstr(char(VBMObj.gser('TIV').gvol('TIV_V2') .toJob));

clear par
par.fname       = 'TIV_V1_n15';
par.foutdir     = main_dir;
par.run         = 1;
job_cat_TIV_estimate(catreport_xml_V1,par);

clear par
par.fname       = 'TIV_V2_n15';
par.foutdir     = main_dir;
par.run         = 1;
job_cat_TIV_estimate(catreport_xml_V2,par);


TIV_reg{1} = load('TIV_V1_n15.txt');
TIV_reg{2} = load('TIV_V2_n15.txt');


covars{3} = TIV_reg{1};
% covars{3} = TIV_reg{2} - TIV_reg{1};

covars_a{3} = TIV_reg{1}(1:8);
covars_c{3} = TIV_reg{1}(9:15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Multiregression V2-V1 (n=15)    
    
par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'VBM_delta_reg_model_spec';
par.mem = 10000;
par.nb_cond = length(outdirs);
par.nb_cons = length(outdirs{1});
target_regressors.name  = {Stat.name};
target_regressors.value = {targetsAPA,targetsDA,targetsSS, targetAxial,targetGabs,targetUPDRS,targetUPDRSIII_Axial,targetUPDRSIII_Sup}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

deltaregression_vbm_model_spec(outdirs, gcons_a, gcons_c, covars, target_regressors, par)
multiregression_vbm_model_spec(outdirs, {gcons_a}, {gcons_c}, covars, target_regressors, par);


%% Models estimation
clear par
cd (main_dir)

for iout = 1 : length(outdirs)
    fspm = addsuffixtofilenames(outdirs{iout}, 'SPM.mat');
%     multifspm = addsuffixtofilenames(multioutdirs{iout}, 'SPM.mat');

    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname = 'VBM_delta_reg_est';
    par.mem = 10000;

    par.jobname  = sprintf('spm_reg_model_est_%s',target_regressors.name{iout});
    job_first_level_estimate(fspm,par)

%     par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{iout});
%     job_first_level_estimate(multifspm,par)
end

%% Contrast creation for each SPM.mat
% F-statistics

    Diff_effect = [0 0 0 0 1 -1];
    Main_effect = [0 0 0 0 0 1
                   0 0 0 0 0 -1];

for iout = 1 : length(outdirs)
%fspm = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    modest = addsuffixtofilenames(outdirs{iout},'SPM.mat');

    %% Contrast names
    contrast_F.names = {
        sprintf('Main effect_%s_on_delta_VBM',target_regressors.name{iout})
        sprintf('Diff effect_%s_on_delta_VBM',target_regressors.name{iout})}';

    %% Contrast values
    contrast_F.values = {
        Main_effect
        Diff_effect}';

    %% Contrast type
    contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

    contrast.names  = [contrast_F.names];
    contrast.values = [contrast_F.values];
    contrast.types  = [contrast_F.types];

    %% Contrast : write
    clear par

    par.sge = 0;
    par.run = 1;
    par.display = 0;
    par.jobname = sprintf('spm_write_%s_VBM_delta_con',target_regressors.name{iout});

    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;

    job_first_level_contrast(modest,contrast,par);

end

%% Contrast creation for each SPM.mat
% F-statistics
    PosCorrelation = [0 0 0 0 1] ;
    NegCorrelation = [0 0 0 0 -1] ;

for iout = 1 : length(outdirs)
    modest = addsuffixtofilenames(outdirs{iout},'SPM.mat');
        
    %% Contrast names
    contrast_T.names = {
        sprintf('Pos correlation_%s_on_VBM_delta',target_regressors.name{iout})
        sprintf('Neg correlation_%s_on_VBM_delta',target_regressors.name{iout})}';

    %% Contrast values
    contrast_T.values = {
        PosCorrelation
        NegCorrelation}';

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
    par.jobname = sprintf('spm_write_%s_con',target_regressors.name{iout});

    % par.sessrep = 'both';
    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;

    job_first_level_contrast(modest,contrast,par);
        
end

%% Multiregression V2-V1 (n=15) per group    
%% KINECT    
par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'VBM_delta_reg_model_spec_a';
par.mem = 10000;
par.nb_cond = length(outdirs_a);
par.nb_cons = length(outdirs_a{1});
target_regressors_a.name  = {Stat.name};
target_regressors_a.value = {targetsAPAK,targetsDAK,targetsSSK, targetAxialK,targetGabsK,targetUPDRSK,targetUPDRSIII_AxialK,targetUPDRSIII_SupK}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

multiregression_vbm_model_spec(outdirs_a, gcons_a, covars_a, target_regressors_a, par)

%% Models estimation
clear par
cd (main_dir)

for iout = 1 : length(outdirs_a)
    fspm_a = addsuffixtofilenames(outdirs_a{iout}, 'SPM.mat');
%     multifspm = addsuffixtofilenames(multioutdirs{iout}, 'SPM.mat');

    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname = 'VBM_delta_reg_est';
    par.mem = 10000;

    par.jobname  = sprintf('spm_reg_model_est_%s',target_regressors_a.name{iout});
    job_first_level_estimate(fspm_a,par)

%     par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{iout});
%     job_first_level_estimate(multifspm,par)
end


%% Contrast creation for each SPM.mat
% F-statistics
    PosCorrelation = [0 0 0 0 1] ;
    NegCorrelation = [0 0 0 0 -1] ;

for iout = 1 : length(outdirs_a)
    modest = addsuffixtofilenames(outdirs_a{iout},'SPM.mat');
        
    %% Contrast names
    contrast_T.names = {
        sprintf('Pos correlation_%s_on_VBM_delta',target_regressors_a.name{iout})
        sprintf('Neg correlation_%s_on_VBM_delta',target_regressors_a.name{iout})}';

    %% Contrast values
    contrast_T.values = {
        PosCorrelation
        NegCorrelation}';

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
    par.jobname = sprintf('spm_write_%s_con',target_regressors_a.name{iout});

    % par.sessrep = 'both';
    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;

    job_first_level_contrast(modest,contrast,par);
        
end

%% ORDINATEUR    
par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.jobname = 'VBM_delta_reg_model_spec_c';
par.mem = 10000;
par.nb_cond = length(outdirs_c);
par.nb_cons = length(outdirs_c{1});
target_regressors_c.name  = {Stat.name};
target_regressors_c.value = {targetsAPAO,targetsDAO,targetsSSO, targetAxialO,targetGabsO,targetUPDRSO,targetUPDRSIII_AxialO,targetUPDRSIII_SupO}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

multiregression_vbm_model_spec(outdirs_c, gcons_c, covars_c, target_regressors_c, par)

%% Models estimation
clear par
cd (main_dir)

for iout = 1 : length(outdirs_c)
    fspm_c = addsuffixtofilenames(outdirs_c{iout}, 'SPM.mat');
%     multifspm = addsuffixtofilenames(multioutdirs{iout}, 'SPM.mat');

    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname = 'VBM_delta_reg_est';
    par.mem = 10000;

    par.jobname  = sprintf('spm_reg_model_est_%s',target_regressors_c.name{iout});
    job_first_level_estimate(fspm_c,par)

%     par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{iout});
%     job_first_level_estimate(multifspm,par)
end


%% Contrast creation for each SPM.mat
% F-statistics
    PosCorrelation = [0 0 0 0 1] ;
    NegCorrelation = [0 0 0 0 -1] ;

for iout = 1 : length(outdirs_c)
    modest = addsuffixtofilenames(outdirs_c{iout},'SPM.mat');
        
    %% Contrast names
    contrast_T.names = {
        sprintf('Pos correlation_%s_on_VBM_delta',target_regressors_c.name{iout})
        sprintf('Neg correlation_%s_on_VBM_delta',target_regressors_c.name{iout})}';

    %% Contrast values
    contrast_T.values = {
        PosCorrelation
        NegCorrelation}';

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
    par.jobname = sprintf('spm_write_%s_con',target_regressors_c.name{iout});

    % par.sessrep = 'both';
    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;

    job_first_level_contrast(modest,contrast,par);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiregression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiregression V1 (n=19)    
    
par.run = 1;
par.sge = 0;
par.sge_queu = 'normal,bigmem';
par.jobname = 'VBM_delta_reg_model_spec';
par.mem = 10000;
par.nb_cond = length(multioutdirs);
par.nb_cons = length(multioutdirs{1});
target_regressors.name  = {MultiStat.name};
%target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1,targetUPDRSIII_Sup1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)
target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1,targetUPDRSIII_Sup1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)
%target_regressors.value = {targetsAPA, targetsDA, targetAxial, targetGabs};
%cellstr(multicons_a{1}.path)
%cons_c

multiregression_vbm_model_spec(multioutdirs, [cellstr(char(multigcons_a{:}));cellstr(char(multigcons_c{:}))], covars_V1, target_regressors, par);
%multiregression_vbm_model_spec(multioutdirs, VBM_a, VBM_c, covars_multi, target_regressors, par)
%regression_model_spec(multioutdirs, {multigcons_a}, {multigcons_c}, covars_V1, target_regressors, par)


%% Models estimation
clear par
cd (main_dir)

for iout = 1 : length(multioutdirs)
    
    multifspm = addsuffixtofilenames(multioutdirs{iout}, 'SPM.mat');

    par.run = 1;
    par.sge = 0;
    par.sge_queu = 'normal,bigmem';
    par.jobname = 'VBM_delta_reg_est';
    par.mem = 10000;

    par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{iout});
    job_first_level_estimate(multifspm,par)
end

%% Contrast creation for each SPM.mat
% F-statistics
    PosCorrelation = [0 0 0 0 1] ;
    NegCorrelation = [0 0 0 0 -1] ;

for iout = 1 : length(multioutdirs)
    modest = addsuffixtofilenames(multioutdirs{iout},'SPM.mat');
        
    %% Contrast names
    contrast_T.names = {
        sprintf('Pos correlation_%s_on_VBM_V1',target_regressors.name{iout})
        sprintf('Neg correlation_%s_on_VBM_V1',target_regressors.name{iout})}';

    %% Contrast values
    contrast_T.values = {
        PosCorrelation
        NegCorrelation}';

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
    par.jobname = sprintf('spm_write_%s_con',target_regressors.name{iout});

    % par.sessrep = 'both';
    par.sessrep = 'none';

    par.delete_previous = 1;
    par.report          = 0;

    job_first_level_contrast(modest,contrast,par);
        
end

