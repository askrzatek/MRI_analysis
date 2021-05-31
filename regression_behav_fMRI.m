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

ACTION = 1;
RS = 0;

% define input directory
InputfMRI = fullfile(main_dir,'/doublerun_sts_tapas_resliced');
InputRS   = fullfile(main_dir,'firstlevel_RS');


%% Regressors definition
%% AGE
covars{1} = [70
74
64
76
79
61
75
66
72
68
68
72
56
57];
         
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
1];

%% Target regressors definition & their outdirs
% define or create output directory
mkdir(main_dir,'resliced_multiple_regression')
MultiRegDir = fullfile(main_dir,'resliced_multiple_regression')

if ACTION
    RegDirC = fullfile(main_dir, '/resliced_ACT_clinic');
    MultiRegDirC = fullfile(MultiRegDir, '/ACT_clinic');
    Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII_AXIAL'};
    
    RegDirG = fullfile(main_dir, '/resliced_ACT_gait/Spontaneous');
    MultiRegDirG = fullfile(MultiRegDir, '/ACT_gait/Spontaneous');
    Gait  = {'APA_AP','DA','Step_Size'};
    %Gait_up  = {'Rapid','Spontaneous'}; % the Rapid gait condition is one subject shorter - needs a separate processing
end

if RS
    RegDirC = fullfile(main_dir, '/resliced_RS_clinic');
    MultiRegDirC = fullfile(MultiRegDir, '/RS_clinic');
    Clinic = {'AXIAL','GABS','UPDRSIII','UPDRSIII_AXIAL'};
    
    RegDirG = fullfile(main_dir, '/resliced_RS_gait/Spontaneous');
    MultiRegDirG = fullfile(MultiRegDir, '/RS_gait/Spontaneous');
    Gait  = {'APA_AP','DA','Step_Size'};
    %Gait_up  = {'Rapid','Spontaneous'}; % the Rapid gait condition is one subject shorter - needs a separate processing
end

RegDir = horzcat({RegDirG}, {RegDirC});
MultiRegDir = horzcat({MultiRegDirG}, {MultiRegDirC});
model = horzcat({Gait},{Clinic});
for imodel = 1:2
    targetmodel = model{imodel};
    for ivar = 1 : length(targetmodel)
       mkdir(RegDir{imodel},sprintf('%s',targetmodel{ivar})); % universal creation of directories depending on the chosen model
       mkdir(MultiRegDir{imodel},sprintf('%s',targetmodel{ivar}));
        % StatDir{ivar} = fullfile(RegDir,Gait{ivar});
        StatObjG{imodel}{ivar} = exam(RegDir{imodel}, sprintf('^[r,s]%s',targetmodel{ivar})); % the Spontaneous gait condition is one subject shorter - needs a separate processing
        StatObjC{imodel}{ivar} = exam(RegDir{imodel},sprintf('^%s$',targetmodel{ivar})); % we have two of them containing UPDRSIII & AXIAL
        MultiStatObjG{imodel}{ivar} = exam(MultiRegDir{imodel}, sprintf('^%s',targetmodel{ivar})); % the Spontaneous gait condition is one subject shorter - needs a separate processing
        MultiStatObjC{imodel}{ivar} = exam(MultiRegDir{imodel},sprintf('^%s$',targetmodel{ivar})); % we have two of them containing UPDRSIII & AXIAL
    end
end

Stat = StatObjG{1}{1} + StatObjG{1}{2} + StatObjG{1}{3} + StatObjC{2}{1} + StatObjC{2}{2} + StatObjC{2}{3} + StatObjC{2}{4};
Stat.explore

MultiStat = MultiStatObjG{1}{1} + MultiStatObjG{1}{2} + MultiStatObjG{1}{3} + MultiStatObjC{2}{1} + MultiStatObjC{2}{2} + MultiStatObjC{2}{3} + MultiStatObjC{2}{4};
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
35.176597489];

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
0.1773839842];
        
%% Spontan Step_Size
targetsSS1 = [257.3421060198
461.4812761193
249.806860919
388.3478308826
163.3475358607
275.7546680306
315.6876356933
371.2412257681
230.2511033831
302.5132704199
262.8219167263
340.6669320703
421.525704218
447.2221013178];        

%% Rapid APA_AP
targetrAPA1 = [56.1860003179
42.202244009
19.9279281562
38.3572741643
43.1886233663
63.3029794047
57.0992020159
61.9964561812
42.2617161059
53.1842418093
51.3589704894
61.1835358307
51.1603131288];
%% Rapid DA
targetrDA1 = [0.1781937653
0.1892911011
0.2379586858
0.2200210944
0.1808432371
0.1457373272
0.2024911505
0.1217863129
0.2044142615
0.1992642898
0.1682406621
0.1644736842
0.1373089983];
        
%% Rapid Step_Size
targetrSS1 = [385.33658393
400.9512381886
326.7006392201
428.223371688
398.4146093581
502.8718930147
484.1753627556
200.2581921664
424.4445720736
407.5565252707
513.2893254134
507.5458941188
517.2496431319];        

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
2];
           
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
18];   

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
31];        
           
%% UPDRSIII-AXIAL
targetUPDRSIII_Axial1 = [15
18
21
18
26
9
19
12
20
29
28
43
17
29];
          
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
            9.8691469678];

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
            0.0004295012];
        
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
            -25.6892601563];        

%% Rapid APA_AP
targetrAPA= [2.3987452788
            -7.7653542778
            5.7803799484
            18.0505768971
            5.3649924749
            14.7903939438
            3.5444038378
            6.3540641136
            12.0672220288
            12.3199655882
            -2.705904908
            -15.8568741964
            11.3793264819];
%% Rapid DA
targetrDA = [-0.017665519
            -0.0123177476
            -0.021931571
            -0.0201645819
            0.0013393699
            0.0327550936
            0.00503371
            -0.0117976353
            0.0112543503
            -0.0114459536
            0.004957791
            0.0274072277
            0.0137945671];
        
%% Rapid Step_Size
targetrSS = [42.2100452713
            56.8306003683
            -16.3692344288
            38.2421823408
            -29.4672299092
            12.7926055023
            -12.2378975762
            77.203440484
            22.4424770158
            80.2671310068
            -98.195711282
            -66.8538622411
            -30.2765261511];        

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
               0];
           
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
               -5];   

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
               4];        
           
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
                        4];
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACTION
%% Choice of fMRI conditions //scans, condition names 
if ACTION
    Conditions  = {'IL','IR','RL','RR'};
    Sessions    = {'V1','V2','V2_V1'};
% input different groups distinction and paths for the cons we need
    fMRIObj = exam(InputfMRI,'firstlevel'); % check directory names
    fMRIObj.addSerie('PARK.*a$','a',8)
    fMRIObj.addSerie('PARK.*[c,DD]$','c',6)

% add volumes for each contrast - each individual
% group A
    fMRIObj.getSerie('a').addVolume('con.*13','RL_V1')
    fMRIObj.getSerie('a').addVolume('con.*14','RL_V2')
    fMRIObj.getSerie('a').addVolume('con.*15','RR_V1')
    fMRIObj.getSerie('a').addVolume('con.*16','RR_V2')
    
    fMRIObj.getSerie('a').addVolume('con.*17','IL_V1')
    fMRIObj.getSerie('a').addVolume('con.*18','IL_V2')
    fMRIObj.getSerie('a').addVolume('con.*19','IR_V1')
    fMRIObj.getSerie('a').addVolume('con.*20','IR_V2')
    
    fMRIObj.getSerie('a').addVolume('con.*21','RR_V2_V1')
    fMRIObj.getSerie('a').addVolume('con.*22','IR_V2_V1')
    fMRIObj.getSerie('a').addVolume('con.*23','RL_V2_V1')
    fMRIObj.getSerie('a').addVolume('con.*24','IL_V2_V1')
    
%group C
    fMRIObj.getSerie('c').addVolume('con.*13','RL_V1')
    fMRIObj.getSerie('c').addVolume('con.*14','RL_V2')
    fMRIObj.getSerie('c').addVolume('con.*15','RR_V1')
    fMRIObj.getSerie('c').addVolume('con.*16','RR_V2')
    
    fMRIObj.getSerie('c').addVolume('con.*17','IL_V1')
    fMRIObj.getSerie('c').addVolume('con.*18','IL_V2')
    fMRIObj.getSerie('c').addVolume('con.*19','IR_V1')
    fMRIObj.getSerie('c').addVolume('con.*20','IR_V2')
    
    fMRIObj.getSerie('c').addVolume('con.*21','RR_V2_V1')
    fMRIObj.getSerie('c').addVolume('con.*22','IR_V2_V1')
    fMRIObj.getSerie('c').addVolume('con.*23','RL_V2_V1')
    fMRIObj.getSerie('c').addVolume('con.*24','IL_V2_V1')
    
end

if RS
    %% get the RS inputfiles in the ref directory
    patient_list = {'PARKGAMEII_001_NB_a','PARKGAMEII_002_BM_a','PARKGAMEII_003_SM_c','PARKGAMEII_007_SD_a','PARKGAMEII_008_JR_a','PARKGAMEII_023_LJ_c','PARKGAMEII_025_CA_a','PARKGAMEII_028_PC_c','PARKGAMEII_033_DD','PARKGAMEII_039_KM_a','PARKGAMEII_043_PD_a','PARKGAMEII_044_CK_c','PARKGAMEII_047_BF_c','PARKGAMEII_048_SB_a'};

    for ipatient = 1: length(patient_list)
        mkdir(InputRS, patient_list{ipatient});
        patients_dir{ipatient} = get_subdir_regex(InputRS, patient_list{ipatient}); %%outdirs
    end
    
    % input different groups distinction and paths for the cons we need
    %RSObj = exam(InputRS,'firstlevel'); % check directory names
%     RSObj.addSerie('PARK.*a$','a',8)
%     RSObj.addSerie('PARK.*[c,DD]$','c',6)
    
    RSObj_a = exam(InputRS,'PARK.*a$');
    RSObj_c = exam(InputRS,'PARK.*[c,DD]$');
    
    ROIs = {'Caudate_L','Caudate_R','Cereb3_L','Cereb3_R','Cereb4_5_L','Cereb4_5_R','Cereb6_L','Cereb6_R','Cereb7b_L','Cereb7b_R','Cereb8_L','Cereb8_R','Cereb9_L','Cereb9_R','Cereb10_L','Cereb10_R','Cereb_Crus1_L','Cereb_Crus1_R','Cereb_Crus2_L','Cereb_Crus2_R','Cuneus_L','Cuneus_R','Insula_L','Insula_R','Pallidum_L','Pallidum_R','Parietal_Inf_L','Parietal_Inf_R','Parietal_Sup_L','Parietal_Sup_R','PCC','Postcentral_L','Postcentral_R','PPN_L','PPN_R','Precentral_L','Precentral_R','Precuneus_L','Precuneus_R','Putamen_L','Putamen_R','SMA_Ant_L','SMA_Ant_R','SMA_Post_L','SMA_Post_R','Thalamus_L','Thalamus_R','Vermis1_2','Vermis3','Vermis4_5','Vermis6','Vermis7','Vermis8','Vermis9','Vermis10'};
    Conditions  = ROIs;
    Sessions    = {'V1','V2'}%,'V2_V1'};
    for iR = 1 : length(ROIs)
        for iS = 1 : length(Sessions)
            %mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
            RSObj_a.mkdir(sprintf('%s_%s',ROIs{iR},Sessions{iS}));
            RSObj_a.addSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS}),sprintf('%s_%s',ROIs{iR},Sessions{iS}));
            
            % GROUP A add volumes for each contrast - each individual
            RSObj_a.getSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS})).addVolume('con.*01',sprintf('%s_%s',ROIs{iR},Sessions{iS}))
            
            RSObj_c.mkdir(sprintf('%s_%s',ROIs{iR},Sessions{iS}));
            RSObj_c.addSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS}),sprintf('%s_%s',ROIs{iR},Sessions{iS}));
            
            % GROUP C add volumes for each contrast - each individual
            RSObj_c.getSerie(sprintf('%s_%s',ROIs{iR},Sessions{iS})).addVolume('con.*01',sprintf('%s_%s',ROIs{iR},Sessions{iS}))
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
        multicons_a{c} = fMRIObj.getSerie('a').getVolume(sprintf('^%s_V1$',Conditions{iC})) %.toJob;
        multigcons_a{c} = cellstr(multicons_a{c}(:) .toJob)
        multicons_c{c} = fMRIObj.getSerie('c').getVolume(sprintf('^%s_V1$',Conditions{iC})) %.toJob;
        multigcons_c{c} = cellstr(multicons_c{c}(:) .toJob)
    end
    if RS
        multiRS_a{c}   = RSObj_a.getSerie(sprintf('%s_V1$',Conditions{iC})).getVolume(sprintf('^%s_V1$',Conditions{iC})) %.toJob;
        multigRS_a{c} = cellstr(multiRS_a{c}(:) .toJob)            
        multiRS_c{c}  = RSObj_c.getSerie(sprintf('%s_V1$',Conditions{iC})).getVolume(sprintf('^%s_V1$',Conditions{iC})) %.toJob;
        multigRS_c{c} = cellstr(multiRS_c{c}(:) .toJob)
    end
    c = c + 1;
    for iS = 1 : length(Sessions)
        %mkdir(StatDir{:},char(sprintf('%s_%s',Conditions{iC},Sessions{iS})))
        Stat.mkdir(sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        Stat.addSerie(sprintf('^%s_%s$',Conditions{iC},Sessions{iS}),sprintf('%s_%s',Conditions{iC},Sessions{iS}));
        
        % cons a & cons_c to be found with the
        % sprintf('%s_%s',Conditions{iC},Sessions{iS}) match and a
        % structure needed for the output of it (not to lose any cons by replacing them)
        % we can use the getSerie.getVol .toJob with the regex of the above 
        if ACTION    
            cons_a{n} = fMRIObj.getSerie('a').getVolume(sprintf('^%s_%s$',Conditions{iC},Sessions{iS})) %.toJob;
            gcons_a{n} = cellstr(cons_a{n}(:) .toJob)
            cons_c{n} = fMRIObj.getSerie('c').getVolume(sprintf('^%s_%s$',Conditions{iC},Sessions{iS})) %.toJob;
            gcons_c{n} = cellstr(cons_c{n}(:) .toJob)
        end
        if RS
            RS_a{n}   = RSObj_a.getSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS})).getVolume(sprintf('%s_%s',Conditions{iC},Sessions{iS})) %.toJob;
            gRS_a{n} = cellstr(RS_a{n}(:) .toJob)            
            RS_c{n}  = RSObj_c.getSerie(sprintf('%s_%s',Conditions{iC},Sessions{iS})).getVolume(sprintf('%s_%s',Conditions{iC},Sessions{iS})) %.toJob;
            gRS_c{n} = cellstr(RS_c{n}(:) .toJob)
        end
        n = n + 1;
    end
end
if ACTION
    outdirs = Stat.getSerie('[R,I]') .toJob; % 7 x 1 x 12 cells
    multioutdirs = MultiStat.getSerie('[R,I]') .toJob;
end
if RS
    outdirs = Stat.getSerie('.*') .toJob; % 7 x 1 x 108 cells
    multioutdirs = MultiStat.getSerie('.*') .toJob;
end

addpath /home/anna.skrzatek/MRI_analysis/

if ACTION
%% Models specification
    
    par.run      = 1;
    par.sge      = 0;
    
    par.jobname = 'ACT_reg_model_spec';
    par.nb_cond = length(outdirs);
    par.nb_cons = length(outdirs{1});
    target_regressors.name  = {Stat.name};
    target_regressors.value = {targetsAPA,targetsDA,targetsSS, targetAxial,targetGabs,targetUPDRS,targetUPDRSIII_Axial}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)
%     target_regressors.value = {targetrAPA,targetrDA,targetrSS};
    %cellstr(cons_a{1}.path)
    %cons_c
    
    regression_model_spec(outdirs, gcons_a, gcons_c, covars, target_regressors, par)
    
%% Multiregression V1 all in    
    
    par.run = 0;
    par.sge = 1;
    par.jobname = 'ACT_multireg_model_spec';
    par.nb_cond = length(multioutdirs);
    par.nb_cons = length(multioutdirs{1});
    target_regressors.name  = {MultiStat.name};
    target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

    %cellstr(multicons_a{1}.path)
    %cons_c
    
    multiregression_model_spec(multioutdirs, multigcons_a, multigcons_c, covars, target_regressors, par)
    
    
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS
%% Models specification
if RS
    
    par.run      = 1;
    par.sge      = 0;
    
    par.jobname = 'RS_reg_model_spec';
    par.nb_cond = length(outdirs);
    par.nb_cons = length(outdirs{1});
    target_regressors.name  = {Stat.name};
    target_regressors.value = {targetsAPA,targetsDA,targetsSS, targetAxial,targetGabs,targetUPDRS,targetUPDRSIII_Axial}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

    regression_model_spec(outdirs, gRS_a, gRS_c, covars, target_regressors, par)
    
    %% Multiregression V1 all in    
    
    par.run = 0;
    par.sge = 1;
    par.jobname = 'RS_multireg_model_spec';
    par.nb_cond = length(multioutdirs);
    par.nb_cons = length(multioutdirs{1});
    target_regressors.name  = {MultiStat.name};
    target_regressors.value = {targetsAPA1,targetsDA1,targetsSS1, targetAxial1,targetGabs1,targetUPDRS1,targetUPDRSIII_Axial1}; % get variables from the variable name regex or get all variables in the same structure before and then just search by their name index (being the same as their index in the structure)

    %cellstr(multicons_a{1}.path)
    %cons_c
    
    multiregression_model_spec(multioutdirs, multigRS_a, multigRS_c, covars, target_regressors, par)


%%
 
end

%% Models estimation
clear par
cd (main_dir)

for iout = 1 : length(outdirs)
    fspm = addsuffixtofilenames(outdirs{iout}, 'SPM.mat');
    multifspm = addsuffixtofilenames(multioutdirs{iout}, 'SPM.mat');

    par.run = 1;
    %par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.jobname  = sprintf('spm_reg_model_est_%s',target_regressors.name{iout});
    job_first_level_estimate(fspm,par)

    par.jobname  = sprintf('spm_multireg_model_est_%s',target_regressors.name{iout});
    job_first_level_estimate(multifspm,par)
end

%% Contrast creation for each SPM.mat
% F-statistics
    Diff_effect = [0 0 0 0 1 -1];
    Main_effect = [0 0 0 0 1 0
                   0 0 0 0 0 1];

for iout = 1 : length(outdirs)
    %fspm = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    modest = addsuffixtofilenames(outdirs{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_F.names = {
            sprintf('Main effect_%s_on_%s',target_regressors.name{iout},roilabel)
            sprintf('Diff effect_%s_on_%s',target_regressors.name{iout},roilabel)}';

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
        par.jobname = sprintf('spm_write_%s_%s_con',target_regressors.name{iout},roilabel);

        % par.sessrep = 'both';
        par.sessrep = 'none';

        par.delete_previous = 1;
        par.report          = 0;

        job_first_level_contrast(modest(iroi),contrast,par);
        
        Stat(iout).getSerie(roilabel).addVolume('spmF_0001','main',1)
        Stat(iout).getSerie(roilabel).addVolume('spmF_0002','diff',1)
        mainef = Stat(iout).getSerie(roilabel).getVolume('main') .toJob
        diffef = Stat(iout).getSerie(roilabel).getVolume('diff') .toJob
        mask{iroi} = cellstr(fullfile(outdirs{iout}{iroi},'mask.nii'));
        
        %% pTFCE toolbox for all con_001 & con_002 in our outdirs
        addpath /network/lustre/iss01/cenir/software/irm/spm12/toolbox/pTFCE/
        
%         [mpTFCE_Z, mpTFCE_p] = pTFCE_adapt(modest{iroi}, char(mainef))
%         [dpTFCE_Z, dpTFCE_p] = pTFCE_adapt(modest{iroi}, char(diffef))
%         still testing : error Assignment has fewer non-singleton rhs dimensions than non-singleton subscripts Error in pTFCE (line 179) PVC(:,:,:,hi)=pvc;

% %       [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask, Rd, V, Nh, Zest, C, verbose)        
% %       load (char(modest{iroi}))
% %       rD = SPM.xVol.R(4);
% %       V = SPM.xVol.S;
% %       [pTFCE_Z, pTFCE_p] = pTFCE(char(mainef),mask, rD, V)
% %       [pTFCE_Z, pTFCE_p] = pTFCE(char(diffef),mask, rD, V)
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiregression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contrast creation for each SPM.mat
% F-statistics
    Correlation = [1 1 0 0] ;

for iout = 1 : length(multioutdirs)
    modest = addsuffixtofilenames(multioutdirs{iout},'SPM.mat');
    for iroi = 1 : length(modest)
        parts = strsplit(char(modest(iroi)), '/');
        roilabel = parts{end-1};   
    
        %% Contrast names
        contrast_T.names = {
            sprintf('Correlation_%s_on_%s',target_regressors.name{iout},roilabel)}';

        %% Contrast values
        contrast_T.values = {
            Correlation}';

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
        
        Stat(iout).getSerie(roilabel).addVolume('spmT_0001','corel',1)
        corelef = Stat(iout).getSerie(roilabel).getVolume('corel') .toJob
        mask{iroi} = cellstr(fullfile(outdirs{iout}{iroi},'mask.nii'));
        
        %% pTFCE toolbox for all con_001 & con_002 in our outdirs
        addpath /network/lustre/iss01/cenir/software/irm/spm12/toolbox/pTFCE/
        
%         [mpTFCE_Z, mpTFCE_p] = pTFCE_adapt(modest{iroi}, char(corelef))
%         still testing : error Assignment has fewer non-singleton rhs dimensions than non-singleton subscripts Error in pTFCE (line 179) PVC(:,:,:,hi)=pvc;

% %       [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask, Rd, V, Nh, Zest, C, verbose)        
% %       load (char(modest{iroi}))
% %       rD = SPM.xVol.R(4);
% %       V = SPM.xVol.S;
% %       [pTFCE_Z, pTFCE_p] = pTFCE(char(corelef),mask, rD, V)
        
    end
end