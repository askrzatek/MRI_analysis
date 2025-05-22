%% Script to Permutation testing for RS seed-to-seed analysis

clc

addpath('/network/iss/cenir/analyse/irm/users/salim.ouarab/toolbox/functions/')

%% Uploading the table
main_dir = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/Results';

%% CEREBELLUM all categories
% Xmat = readtable(fullfile(main_dir, 'RS_seed2seed/Xmat_Tin_Con_ROI2ROI.csv'));
% tt = get_subdir_regex_files(main_dir,'AUDICOG_17_01_25_distribution_ROI2ROI_Salience.csv');
%Xmat = readtable(fullfile(main_dir, 'AUDICOG_17_01_25_distribution_ROI2ROI_Salience.csv'));
Xmat = readtable(fullfile(main_dir, 'DysPPN_05_05_25_distribution_ROI2ROI_PPN_Cerebellum.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Dys_FC = Xmat.(varname)(1:19);
    Con_FC = Xmat.(varname)(20:41);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Dys_FC, Con_FC, 10000);% , 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

DysPPN_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/Results')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_salience_net_permutation_test.csv');
% writetable(DysPPN_permut_test_tab,'DysPPN_seed2seed_PPN_Cerebellum_permutation_test_05_05_25.csv');
writetable(DysPPN_permut_test_tab,'DysPPN_seed2seed_PPN_Cerebellum_permutation_test_20_05_25.csv');

%% CORTEX reduced categories /without Limbic and Primary
% Xmat = readtable(fullfile(main_dir, 'RS_seed2seed/Xmat_Tin_Con_ROI2ROI.csv'));
% tt = get_subdir_regex_files(main_dir,'AUDICOG_17_01_25_distribution_ROI2ROI_Salience.csv');
%Xmat = readtable(fullfile(main_dir, 'AUDICOG_17_01_25_distribution_ROI2ROI_Salience.csv'));
Xmat = readtable(fullfile(main_dir, 'DysPPN_05_05_25_distribution_ROI2ROI_PPN_Cortical.csv'));

pval = [];
obs_d = [];
es = [];

for i = 2:length(Xmat.Properties.VariableNames)
    varname = Xmat.Properties.VariableNames{i};
    Dys_FC = Xmat.(varname)(1:19);
    Con_FC = Xmat.(varname)(20:41);
    [p(i-1), observeddifference(i-1), effectsize(i-1)] = permutationTest(Dys_FC, Con_FC, 10000 );%, 'showprogress', 250);

    pval = [pval ; p(i-1)];
    obs_d = [obs_d ; observeddifference(i-1)];
    es = [es ; effectsize(i-1)];
end

ROInames = Xmat.Properties.VariableNames
ROInames = ROInames(2:length(ROInames))

DysPPN_permut_test_tab = table(ROInames.',pval,obs_d,es);

cd ('/network/iss/cenir/analyse/irm/users/anna.skrzatek/DYS_PPN/RSFC/Results')
% writetable(AUDICOG_permut_test_tab,'AUDICOG_seed2seed_salience_net_permutation_test.csv');
% writetable(DysPPN_permut_test_tab,'DysPPN_seed2seed_PPN_Cortex_permutation_test_05_05_25.csv');
writetable(DysPPN_permut_test_tab,'DysPPN_seed2seed_PPN_Cortex_permutation_test_20_05_25.csv');

