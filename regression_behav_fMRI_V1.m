%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script for creating the regression models between fMRI & other statistical measures
%% ex. Clinical, Gait or Game Data regression for each fMRI condition
%% Anna SKRZATEK May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% Initialise

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','/nifti_test');
cd (main_dir)

ACTION = 0;
RS = 1;

% define input directory
patients = exam(main_dir,'PARKGAME.*V1_[a,c]$');
patients.addSerie('rsmodel_ts_tapas','firstlevel',1);
InputfMRI = patients.gser('firstlevel') .toJob; %% we will need contrasts from 8 to 11
InputRS   = fullfile(main_dir,'firstlevel_RS');


%% Regressors definition
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
59              % 053_LM
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
1               % 053_LM
1               % 003_SM_c
2               % 023_LJ
2               % 028_PC_c
2               % 033_DD
2               % 044_CK
1               % 046_HJ
1               % 047_BF
2];             % 052_HJ

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
55.1500901593812                % 053_LM
37.3914711347                   % 003_SM_c
33.3035609785                   % 023_LJ
41.7655772122                   % 028_PC
26.2586959347                   % 033_DD
32.1145722308                   % 044_CK
34.2056328066679                % 046_HJ
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
0.210950764006791               % 053_LM
0.1624504829                    % 003_SM_c
0.2447760219                    % 023_LJ
0.2448816538                    % 028_PC
0.2531833616                    % 033_DD
0.2184634975                    % 044_CK
0.283773951006549               % 046_HJ
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
406.869666413503                % 053_LM
230.2511033831                  % 003_SM_c
302.5132704199                  % 023_LJ
262.8219167263                  % 028_PC
340.6669320703                  % 033_DD
421.525704218                   % 044_CK
336.690629708893                % 046_HJ
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
1                               % 053_LM
3                               % 003_SM_c
4                               % 023_LJ
8                               % 028_PC
7                               % 033_DD
3                               % 044_CK
7                               % 046_HJ
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
31                              % 053_LM
19                              % 003_SM_c
18                              % 023_LJ
32                              % 028_PC
26                              % 033_DD
14                              % 044_CK
21                              % 046_HJ
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
8                               % 053_LM
23                              % 003_SM_c
33                              % 023_LJ
36                              % 028_PC
50                              % 033_DD
20                              % 044_CK
37                              % 046_HJ
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
7                               % 053_LM
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
5                              % 053_LM
14                             % 003_SM_c
2                              % 023_LJ
5                              % 028_PC
5                              % 033_DD
8                              % 044_CK
16                             % 046_HJ ?
4                              % 047_BF
11];                           % 052_HJ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACTION
%% Choice of fMRI conditions //scans, condition names 
if ACTION
    Conditions  = {'IL','IR','RL','RR'};
    Sessions    = {'V1'};

% add volumes for each contrast - each individual
% groups combined
    patients.getSerie('firstlevel').addVolume('con.*08','RL_V1')
    patients.getSerie('firstlevel').addVolume('con.*09','RR_V1')
    
    patients.getSerie('firstlevel').addVolume('con.*10','IL_V1')
    patients.getSerie('firstlevel').addVolume('con.*11','IR_V1')
    
end

if RS
    %% get the RS inputfiles in the ref directory
    patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_025_CA_a','PARKGAMEII_039_KM_a','PARKGAMEII_040_RE_a','PARKGAMEII_042_RS_a','PARKGAMEII_043_PD_a','PARKGAMEII_048_SB_a','PARKGAMEII_053_LM_a','PARKGAMEII_003_SM_c','PARKGAMEII_023_LJ_c','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_044_CK_c','PARKGAMEII_046_HJ_c','PARKGAMEII_047_BF_c','PARKGAMEII_052_HJ_c'};

    for ipatient = 1: length(patient_list)
        mkdir(InputRS, patient_list{ipatient});
        patients_dir{ipatient} = get_subdir_regex(InputRS, patient_list{ipatient}); %%outdirs
    end
    
    % input different groups distinction and paths for the cons we need
    %RSObj = exam(InputRS,'firstlevel'); % check directory names
%     RSObj.addSerie('PARK.*a$','a',8)
%     RSObj.addSerie('PARK.*[c,DD]$','c',6)
    
    RSObj = exam(InputRS,'PARK.*[a,DD,c]$');
    
    ROIs = {'Caudate_L','Caudate_R','Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Pallidum_L','Pallidum_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','PCC','Postcentral_L','Postcentral_R','PPN_L','PPN_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','Putamen_L','Putamen_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Thalamus_L','Thalamus_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10'};
%     ROIs = {'SMA_face_L','SMA_face_R','SMA_foot_L','SMA_foot_R','SMA_hand_L','SMA_hand_R'};
    Conditions  = ROIs;
    Sessions    = {'V1'};
    for iR = 1 : length(ROIs)
        for iS = 1 : length(Sessions)
            %mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
            RSObj.mkdir(sprintf('%s_%s',ROIs{iR},Sessions{iS}));
            RSObj.addSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS}),sprintf('%s_%s',ROIs{iR},Sessions{iS}));
            
            % GROUPS COMBINED add volumes for each contrast - each individual
            RSObj.getSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS})).addVolume('con.*01',sprintf('%s_%s',ROIs{iR},Sessions{iS}))
            
        end
    end

    % do a rsync -rlv copy of all individual model dirs adding V1 or V2 at the end
    % accordingly - then the organisation would be almost the same as for the
    % double run - otherwise - prepare the doublerun analysis for RS data - all
    % you need is separate sessions V1 & V2 and create a structure of 2 runs
    % per subject cf. firstlevel_tedana_doublerun & RS_firstlevel_Cecile

end

%% this part will probably become in common with the RS processing
n = 1;
c = 1;
for iC = 1 : length(Conditions)
    MultiStat.mkdir(sprintf('%s',Conditions{iC}));
    MultiStat.addSerie(sprintf('^%s$',Conditions{iC}),sprintf('%s',Conditions{iC}));
    if ACTION
        multicons{c} = patients.getSerie('firstlevel').getVolume(sprintf('^%s_V1$',Conditions{iC})) %.toJob;
        multigcons{c} = cellstr(multicons{c}(:) .toJob)
    end
    if RS
        multiRS{c}   = RSObj.getSerie(sprintf('%s_V1',Conditions{iC})).getVolume(sprintf('^%s_V1$',Conditions{iC})) %.toJob;
        multigRS{c} = cellstr(multiRS{c}(:) .toJob)            
    end
    c = c + 1;
    iS = 1;
        % cons a & cons_c to be found with the
        % sprintf('%s_%s',Conditions{iC},Sessions{iS}) match and a
        % structure needed for the output of it (not to lose any cons by replacing them)
        % we can use the getSerie.getVol .toJob with the regex of the above 
    if ACTION    
        cons{n} = patients.getSerie('firstlevel').getVolume(sprintf('^%s_%s$',Conditions{iC},Sessions{iS})) %.toJob;
        gcons{n} = cellstr(cons{n}(:) .toJob)
    end
    if RS
        RS_all{n}   = RSObj.getSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS})).getVolume(sprintf('%s_%s',Conditions{iC},Sessions{iS})) %.toJob;
        gRS{n} = cellstr(RS_all{n}(:) .toJob)
    end
    n = n + 1;

end
if ACTION
    multioutdirs = MultiStat.getSerie('[R,I]') .toJob;
end
if RS
    multioutdirs = MultiStat.getSerie('.*') .toJob;
end

addpath /home/anna.skrzatek/MRI_analysis/

if ACTION
%% Models specification 
% Multiregression V1 all in    
    
    par.run = 1;
    par.sge = 0;
    par.jobname = 'ACT_multireg_model_spec';
    par.nb_cond = length(multioutdirs);
    par.nb_cons = length(multioutdirs{1});
    target_regressors.name  = {MultiStat.name};
    target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1, targetUPDRSIII_Sup1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)
    
    multiregressionV1_model_spec(multioutdirs, multigcons, covars, target_regressors, par)
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS
%% Models specification
if RS
% Multiregression V1 all in    
    
    par.run = 0;
    par.sge = 1;
    par.jobname = 'RS_multireg_model_spec';
    par.nb_cond = length(multioutdirs);
    par.nb_cons = length(multioutdirs{1});
    target_regressors.name  = {MultiStat.name};
    target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1, targetUPDRSIII_Sup1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

    %cellstr(multicons_a{1}.path)
    %cons_c
    
    multiregressionV1_model_spec(multioutdirs, multigRS, covars, target_regressors, par)

%%
 
end

%% Models estimation
clear par
cd (main_dir)

for iout = 1 : length(multioutdirs)
    multifspm = addsuffixtofilenames(multioutdirs{iout}, 'SPM.mat');

    par.run = 1;
    %par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    
    par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{iout});
    job_first_level_estimate(multifspm,par)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiregression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contrast creation for each SPM.mat
% F-statistics
    PosCorrelation = [0 0 0 1] ;
    NegCorrelation = [0 0 0 -1] ;

for iout = 1 : length(multioutdirs)
    modest = addsuffixtofilenames(multioutdirs{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_T.names = {
            sprintf('Pos correlation_%s_on_%s',target_regressors.name{iout},roilabel)
            sprintf('Neg correlation_%s_on_%s',target_regressors.name{iout},roilabel)}';

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
        par.jobname = sprintf('spm_write_%s_%s_con',target_regressors.name{iout},roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(modest(iroi),contrast,par);
        
        MultiStat(iout).getSerie(roilabel).addVolume('spmT_0001','positive',1)
        MultiStat(iout).getSerie(roilabel).addVolume('spmT_0002','negative',1)
        poscorel = MultiStat(iout).getSerie(roilabel).getVolume('positive') .toJob
        negcorel = MultiStat(iout).getSerie(roilabel).getVolume('negative') .toJob
        %mask{iroi} = cellstr(fullfile(multioutdirs{iout}{iroi},'mask.nii'));
        
        %% pTFCE toolbox for all con_001 & con_002 in our outdirs
%% % sadly we still don't know how to transform variable to img - computation works, but no file is created

        addpath /network/lustre/iss01/cenir/software/irm/spm12/toolbox/pTFCE/
        
        [PpTFCE_Z, PpTFCE_p] = pTFCE_adapt(modest{iroi}, char(poscorel));
        [NpTFCE_Z, NpTFCE_p] = pTFCE_adapt(modest{iroi}, char(negcorel));
        
    end
end