%% Script to Permutation testing for RS seed-to-seed analysis

clc

addpath('/network/iss/cenir/analyse/irm/users/salim.ouarab/toolbox/functions/')

%% Uploading the table
main_dir = '/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed'

%% SALIENCE NETWORK
% Xmat = readtable(fullfile(main_dir, 'RS_seed2seed/Xmat_Tin_Con_ROI2ROI.csv'));
% tt = get_subdir_regex_files(main_dir,'AUDICOG_17_01_25_distribution_ROI2ROI_Salience.csv');
%Xmat = readtable(fullfile(main_dir, 'AUDICOG_17_01_25_distribution_ROI2ROI_Salience.csv'));
Xmat = readtable(fullfile(main_dir, 'AUDICOG_23_01_25_distribution_ROI2ROI_SN_Tinnitus_Alert_AudioACT.csv'));
Xmat = readtable(fullfile(main_dir, 'AUDICOG_20_05_25_distribution_ROI2ROI_SN_Tinnitus_Alert_AudioACT.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Tin_FC = Xmat.(varname)(1:25);
    Con_FC = Xmat.(varname)(26:50);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Tin_FC, Con_FC, 10000 );%, 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

AUDICOG_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed/')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_salience_net_permutation_test.csv');
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_permutation_test_23_01_25.csv');
%writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_permutation_test_13_05_25.csv');
writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_permutation_test_20_05_25.csv');

%% CENETWORK
Xmat = readtable(fullfile(main_dir, 'AUDICOG_23_01_25_distribution_ROI2ROI_CEN_Tinnitus_Alert_AudioACT.csv'));
Xmat = readtable(fullfile(main_dir, 'AUDICOG_20_05_25_distribution_ROI2ROI_CEN_Tinnitus_Alert_AudioACT.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Tin_FC = Xmat.(varname)(1:25);
    Con_FC = Xmat.(varname)(26:50);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Tin_FC, Con_FC, 10000);% , 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

AUDICOG_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed/')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_CEN_permutation_test_23_01_25.csv');
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_CEN_permutation_test_13_05_25.csv');
writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_CEN_permutation_test_20_05_25.csv');

%% DMNETWORK
Xmat = readtable(fullfile(main_dir, 'AUDICOG_23_01_25_distribution_ROI2ROI_DMN_Tinnitus_Alert_AudioACT.csv'));
Xmat = readtable(fullfile(main_dir, 'AUDICOG_20_05_25_distribution_ROI2ROI_DMN_Tinnitus_Alert_AudioACT.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Tin_FC = Xmat.(varname)(1:25);
    Con_FC = Xmat.(varname)(26:50);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Tin_FC, Con_FC, 10000 );%, 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

AUDICOG_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed/')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_DMN_permutation_test_23_01_25.csv');
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_DMN_permutation_test_13_05_25.csv');
writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_DMN_permutation_test_20_05_25.csv');

%% 3 NETWORKS
Xmat = readtable(fullfile(main_dir, 'AUDICOG_23_01_25_distribution_ROI2ROI_SN_CEN_DMN.csv'));
Xmat = readtable(fullfile(main_dir, 'AUDICOG_20_05_25_distribution_ROI2ROI_SN_CEN_DMN.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Tin_FC = Xmat.(varname)(1:25);
    Con_FC = Xmat.(varname)(26:50);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Tin_FC, Con_FC, 10000);% , 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

AUDICOG_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed/')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_CEN_DMN_permutation_test_23_01_25.csv');
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_CEN_DMN_permutation_test_13_05_25.csv');
writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_CEN_DMN_permutation_test_20_05_25.csv');

%% CONTROL NETWORKS
Xmat = readtable(fullfile(main_dir, 'AUDICOG_03_02_25_distribution_ROI2ROI_SMN_Visual.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Tin_FC = Xmat.(varname)(1:25);
    Con_FC = Xmat.(varname)(26:50);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Tin_FC, Con_FC, 1000 , 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

AUDICOG_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed/')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SMN_Visual_permutation_test_03_02_25.csv');
writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SMN_Visual_permutation_test_13_05_25.csv');

%% ALL NETWORKS
Xmat = readtable(fullfile(main_dir, 'AUDICOG_03_02_25_distribution_ROI2ROI_SN_CEN_DMN_SMN_Visual.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Tin_FC = Xmat.(varname)(1:25);
    Con_FC = Xmat.(varname)(26:50);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Tin_FC, Con_FC, 100000 );%, 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

AUDICOG_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/studies/AUDICOG/Results/RS_seed2seed/')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_CEN_DMN_SMN_Visual_permutation_test_03_02_25.csv');
writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_SN_CEN_DMN_SMN_Visual_permutation_test_13_05_25.csv');

