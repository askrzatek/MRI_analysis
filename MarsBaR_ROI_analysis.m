% Extract Mean Percent Signal Change script.  Problem indicated with comment.
% MarsBar percent activation script
% Before running this do the following
% 1) make sure you have imported your ROIs from WFU and saved them
% 2) make sure you have spm and MarsBaR open and running

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% BASIC STEPS IN ROI ANALYSIS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spm_name = '/my/path/SPM.mat';
% roi_file = '/my/path/my_roi.mat';
% % Make marsbar design 
% objectD = mardo(spm_name);
% % Make marsbar ROI 
% objectR = maroi(roi_file);
% % Fetch data into marsbar data 
% objectY = get_marsy(R, D, 'mean');
% % Get contrasts from original 
% designxCon = get_contrasts(D);
% % Estimate design on ROI 
% dataE = estimate(D, Y);
% % Put contrasts from original design back into design 
% objectE = set_contrasts(E, xCon);
% % get design betasb = betas(E);
% % get stats and stuff for all contrasts into statistics 
% structuremarsS = compute_contrasts(E, 1:length(xCon));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 

%%% #1 Start marsbar to make sure spm_get works

addpath /network/lustre/iss02/cenir/software/irm/spm12/toolbox/marsbar/
marsbar('on')
% Set up the SPM defaults, just in case
spm('defaults', 'fmri');

%% #2 Initialisation

main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek','nifti_test');
%roi_model_dir = fullfile(char(main_dir), 'secondlevel_ACTIVATION_PARK_S1');
%roi_model_dir = fullfile(char(main_dir), 'secondlevel_tedana_ac_PARK');
roi_model_dir = fullfile(char(main_dir), 'secondlevel_sts_tapas_PARK');

% roi_group.regex = {'LEFT_REAL_S1.*_roi.mat', 'LEFT_IMAGINARY_S1.*_roi.mat', 'RIGHT_REAL_S1.*_roi.mat', 'RIGHT_IMAGINARY_S1.*_roi.mat', 'Masked_LEFT_IMAGINARY_REAL_S1.*_roi.mat', 'Masked_RIGHT_IMAGINARY_REAL_S1.*_roi.mat'};
% roi_group.name = {'Main_spe_LEFT_REAL_S1', 'Main_spe_LEFT_IMAGINARY_S1', 'Main_spe_RIGHT_REAL_S1', 'Main_spe_RIGHT_IMAGINARY_S1', 'Masked_LEFT_IMAGINARY-REAL_S1', 'Masked_RIGHT_IMAGINARY-REAL_S1'};

%roi_group.regex = {'Main_spe_LEFT_REAL_S1.*_roi.mat', 'Main_spe_LEFT_IMAGINARY_S1.*_roi.mat', 'Main_spe_RIGHT_REAL_S1.*_roi.mat', 'Main_spe_RIGHT_IMAGINARY_S1.*_roi.mat', 'Main_conj_LEFT_IMAGINARY_REAL_S1.*_roi.mat', 'Main_conj_RIGHT_IMAGINARY_REAL_S1.*_roi.mat'};
%roi_group.name = {'Main_spe_LEFT_REAL_S1', 'Main_spe_LEFT_IMAGINARY_S1', 'Main_spe_RIGHT_REAL_S1', 'Main_spe_RIGHT_IMAGINARY_S1', 'Main_conj_LEFT_IMAGINARY_REAL_S1', 'Main_conj_RIGHT_IMAGINARY_REAL_S1'};

% roi_group.regex = {'LEFT_REAL.*_roi.mat', 'LEFT_IMAGINARY.*_roi.mat', 'RIGHT_REAL.*_roi.mat', 'RIGHT_IMAGINARY.*_roi.mat', 'Masked_LEFT_IMAGINARY_REAL.*_roi.mat', 'Masked_RIGHT_IMAGINARY_REAL.*_roi.mat', 'Masked_LEFT_REAL_IMAGINARY.*_roi.mat', 'Masked_RIGHT_REAL_IMAGINARY.*_roi.mat', 'Conj_LEFT.*_roi.mat', 'Conj_RIGHT.*_roi.mat'};
% roi_group.name = {'Main_spe_LEFT_REAL', 'Main_spe_LEFT_IMAGINARY', 'Main_spe_RIGHT_REAL', 'Main_spe_RIGHT_IMAGINARY', 'Masked_LEFT_IMAGINARY-REAL', 'Masked_RIGHT_IMAGINARY-REAL', 'Masked_LEFT_REAL-IMAGINARY', 'Masked_RIGHT_REAL_IMAGINARY', 'Conj_LEFT_REAL_IMA', 'Conj_RIGHT_REAL_IMA'};

roi_group.regex = {'atlas_1_roi.mat','atlas_2_roi.mat','atlas_3_roi.mat','atlas_4_roi.mat','atlas_5_roi.mat','atlas_6_roi.mat','atlas_7_roi.mat','atlas_8_roi.mat','atlas_9_roi.mat','atlas_10_roi.mat','atlas_11_roi.mat','atlas_12_roi.mat','atlas_13_roi.mat','atlas_14_roi.mat','atlas_15_roi.mat','atlas_16_roi.mat','atlas_17_roi.mat','atlas_18_roi.mat','atlas_19_roi.mat','atlas_20_roi.mat','atlas_21_roi.mat','atlas_22_roi.mat','atlas_23_roi.mat','atlas_24_roi.mat','atlas_25_roi.mat','atlas_26_roi.mat','atlas_27_roi.mat','atlas_28_roi.mat','atlas_29_roi.mat','atlas_30_roi.mat','atlas_31_roi.mat','atlas_32_roi.mat','atlas_33_roi.mat','atlas_34_roi.mat','atlas_35_roi.mat','atlas_36_roi.mat','atlas_37_roi.mat','atlas_38_roi.mat','atlas_39_roi.mat','atlas_40_roi.mat','atlas_41_roi.mat','atlas_42_roi.mat','atlas_43_roi.mat','atlas_44_roi.mat','atlas_45_roi.mat','atlas_46_roi.mat','atlas_47_roi.mat','atlas_48_roi.mat','atlas_49_roi.mat','atlas_50_roi.mat','atlas_51_roi.mat','atlas_52_roi.mat','atlas_53_roi.mat','atlas_54_roi.mat','atlas_55_roi.mat','atlas_56_roi.mat','atlas_57_roi.mat','atlas_58_roi.mat','atlas_59_roi.mat','atlas_60_roi.mat','atlas_61_roi.mat','atlas_62_roi.mat','atlas_63_roi.mat','atlas_64_roi.mat','atlas_65_roi.mat','atlas_66_roi.mat','atlas_67_roi.mat','atlas_68_roi.mat','atlas_69_roi.mat','atlas_70_roi.mat','atlas_71_roi.mat','atlas_72_roi.mat','atlas_73_roi.mat','atlas_74_roi.mat','atlas_75_roi.mat','atlas_76_roi.mat','atlas_77_roi.mat','atlas_78_roi.mat','atlas_79_roi.mat','atlas_80_roi.mat','atlas_81_roi.mat','atlas_82_roi.mat','atlas_83_roi.mat','atlas_84_roi.mat','atlas_85_roi.mat','atlas_86_roi.mat','atlas_87_roi.mat','atlas_88_roi.mat','atlas_89_roi.mat','atlas_90_roi.mat','atlas_91_roi.mat','atlas_92_roi.mat','atlas_93_roi.mat','atlas_94_roi.mat','atlas_95_roi.mat','atlas_96_roi.mat','atlas_97_roi.mat','atlas_98_roi.mat','atlas_99_roi.mat','atlas_100_roi.mat','atlas_101_roi.mat','atlas_102_roi.mat','atlas_103_roi.mat','atlas_104_roi.mat','atlas_105_roi.mat','atlas_106_roi.mat','atlas_107_roi.mat','atlas_108_roi.mat','atlas_109_roi.mat','atlas_110_roi.mat','atlas_111_roi.mat','atlas_112_roi.mat','atlas_113_roi.mat','atlas_114_roi.mat','atlas_115_roi.mat','atlas_116_roi.mat','atlas_117_roi.mat','atlas_118_roi.mat','atlas_119_roi.mat','atlas_120_roi.mat','atlas_121_roi.mat','atlas_122_roi.mat','atlas_123_roi.mat','atlas_124_roi.mat','atlas_125_roi.mat','atlas_126_roi.mat','atlas_127_roi.mat','atlas_128_roi.mat','atlas_129_roi.mat','atlas_130_roi.mat','atlas_131_roi.mat','atlas_132_roi.mat'}
roi_group.name = {'FP r (Frontal Pole Right)','FP l (Frontal Pole Left)','IC r (Insular Cortex Right)','IC l (Insular Cortex Left)','SFG r (Superior Frontal Gyrus Right)','SFG l (Superior Frontal Gyrus Left)','MidFG r (Middle Frontal Gyrus Right)','MidFG l (Middle Frontal Gyrus Left)','IFG tri r (Inferior Frontal Gyrus, pars triangularis Right)','IFG tri l (Inferior Frontal Gyrus, pars triangularis Left)','IFG oper r (Inferior Frontal Gyrus, pars opercularis Right)','IFG oper l (Inferior Frontal Gyrus, pars opercularis Left)','PreCG r (Precentral Gyrus Right)','PreCG l (Precentral Gyrus Left)','TP r (Temporal Pole Right)','TP l (Temporal Pole Left)','aSTG r (Superior Temporal Gyrus, anterior division Right)','aSTG l (Superior Temporal Gyrus, anterior division Left)','pSTG r (Superior Temporal Gyrus, posterior division Right)','pSTG l (Superior Temporal Gyrus, posterior division Left)','aMTG r (Middle Temporal Gyrus, anterior division Right)','aMTG l (Middle Temporal Gyrus, anterior division Left)','pMTG r (Middle Temporal Gyrus, posterior division Right)','pMTG l (Middle Temporal Gyrus, posterior division Left)','toMTG r (Middle Temporal Gyrus, temporooccipital part Right)','toMTG l (Middle Temporal Gyrus, temporooccipital part Left)','aITG r (Inferior Temporal Gyrus, anterior division Right)','aITG l (Inferior Temporal Gyrus, anterior division Left)','pITG r (Inferior Temporal Gyrus, posterior division Right)','pITG l (Inferior Temporal Gyrus, posterior division Left)','toITG r (Inferior Temporal Gyrus, temporooccipital part Right)','toITG l (Inferior Temporal Gyrus, temporooccipital part Left)','PostCG r (Postcentral Gyrus Right)','PostCG l (Postcentral Gyrus Left)','SPL r (Superior Parietal Lobule Right)','SPL l (Superior Parietal Lobule Left)','aSMG r (Supramarginal Gyrus, anterior division Right)','aSMG l (Supramarginal Gyrus, anterior division Left)','pSMG r (Supramarginal Gyrus, posterior division Right)','pSMG l (Supramarginal Gyrus, posterior division Left)','AG r (Angular Gyrus Right)','AG l (Angular Gyrus Left)','sLOC r (Lateral Occipital Cortex, superior division Right)','sLOC l (Lateral Occipital Cortex, superior division Left)','iLOC r (Lateral Occipital Cortex, inferior division Right)','iLOC l (Lateral Occipital Cortex, inferior division Left)','ICC r (Intracalcarine Cortex Right)','ICC l (Intracalcarine Cortex Left)','MedFC (Frontal Medial Cortex)','SMA r (Juxtapositional Lobule Cortex -formerly Supplementary Motor Cortex- Right)','SMA L(Juxtapositional Lobule Cortex -formerly Supplementary Motor Cortex- Left)','SubCalC (Subcallosal Cortex)','PaCiG r (Paracingulate Gyrus Right)','PaCiG l (Paracingulate Gyrus Left)','AC (Cingulate Gyrus, anterior division)','PC (Cingulate Gyrus, posterior division)','Precuneous (Precuneous Cortex)','Cuneal r (Cuneal Cortex Right)','Cuneal l (Cuneal Cortex Left)','FOrb r (Frontal Orbital Cortex Right)','FOrb l (Frontal Orbital Cortex Left)','aPaHC r (Parahippocampal Gyrus, anterior division Right)','aPaHC l (Parahippocampal Gyrus, anterior division Left)','pPaHC r (Parahippocampal Gyrus, posterior division Right)','pPaHC l (Parahippocampal Gyrus, posterior division Left)','LG r (Lingual Gyrus Right)','LG l (Lingual Gyrus Left)','aTFusC r (Temporal Fusiform Cortex, anterior division Right)','aTFusC l (Temporal Fusiform Cortex, anterior division Left)','pTFusC r (Temporal Fusiform Cortex, posterior division Right)','pTFusC l (Temporal Fusiform Cortex, posterior division Left)','TOFusC r (Temporal Occipital Fusiform Cortex Right)','TOFusC l (Temporal Occipital Fusiform Cortex Left)','OFusG r (Occipital Fusiform Gyrus Right)','OFusG l (Occipital Fusiform Gyrus Left)','FO r (Frontal Operculum Cortex Right)','FO l (Frontal Operculum Cortex Left)','CO r (Central Opercular Cortex Right)','CO l (Central Opercular Cortex Left)','PO r (Parietal Operculum Cortex Right)','PO l (Parietal Operculum Cortex Left)','PP r (Planum Polare Right)','PP l (Planum Polare Left)','HG r (Heschl Gyrus Right)','HG l (Heschl Gyrus Left)','PT r (Planum Temporale Right)','PT l (Planum Temporale Left)','SCC r (Supracalcarine Cortex Right)','SCC l (Supracalcarine Cortex Left)','OP r (Occipital Pole Right)','OP l (Occipital Pole Left)','Thalamus r','Thalamus l','Caudate r','Caudate l','Putamen r','Putamen l','Pallidum r','Pallidum l','Hippocampus r','Hippocampus l','Amygdala r','Amygdala l','Accumbens r','Accumbens l','Brain-Stem','Cereb1 l (Cerebelum Crus1 Left)','Cereb1 r (Cerebelum Crus1 Right)','Cereb2 l (Cerebelum Crus2 Left)','Cereb2 r (Cerebelum Crus2 Right)','Cereb3 l (Cerebelum 3 Left)','Cereb3 r (Cerebelum 3 Right)','Cereb45 l (Cerebelum 4 5 Left)','Cereb45 r (Cerebelum 4 5 Right)','Cereb6 l (Cerebelum 6 Left)','Cereb6 r (Cerebelum 6 Right)','Cereb7 l (Cerebelum 7b Left)','Cereb7 r (Cerebelum 7b Right)','Cereb8 l (Cerebelum 8 Left)','Cereb8 r (Cerebelum 8 Right)','Cereb9 l (Cerebelum 9 Left)','Cereb9 r (Cerebelum 9 Right)','Cereb10 l (Cerebelum 10 Left)','Cereb10 r (Cerebelum 10 Right)','Ver12 (Vermis 1 2)','Ver3 (Vermis 3)','Ver45 (Vermis 4 5)','Ver6 (Vermis 6)','Ver7 (Vermis 7)','Ver8 (Vermis 8)','Ver9 (Vermis 9)','Ver10 (Vermis 10)'}

out_tab = {'id', 'group', 'session', 'roi', 'contrast', 'value', 'T', 'pval', 'pvalC'};
nrow = 1; % row index/number in out_tab - keep calm and count
ncol = length(out_tab);
%pct_tab = {'id', 'group', 'session', 'roi', 'event1_name', 'event2_name', 'event3_name', 'event4_name', 'event5_name'};
pct_tab = {'id', 'group', 'session', 'roi', 'event1_name', 'event2_name', 'event3_name', 'event4_name', 'event5_name', 'event6_name'};
pctrow = 1;

% dirgroup = {'PARKGAME','REMINARY'};
% dirgroup = {'PARK_a', 'PARK_c'};
dirgroup = {'PARKGAME_all'};

par.group = 1;

%% #3a THE ROIS LOOP - getting the existing contrast-related ROIs

for roic =1 : length(roi_group.name)
    
    par.conname = roi_group.regex{roic};
    par.subdir = 'ANOVA2x2_LxT_S1';

    %test_batch(par); % batch spm for loading a contrast and saving it in .mat in model dir %turns out not useful at all

    %close all

    %%
    if par.group == 1 % 
       roi_model_dir = fullfile(char(roi_model_dir), dirgroup{1});
       par.group = 0;
    end
    model_dir = fullfile(char(roi_model_dir), par.subdir);
    
    %rois_dir = get_subdir_regex(char(model_dir),'k10_p01');
    %rois_dir = get_subdir_regex(char(model_dir),'k10_p01');
    %rois_dir = get_subdir_regex(char(model_dir),'rois_S1_p05');
    %rois_dir = get_subdir_regex(char(roi_model_dir),'rois_all_S1_p001_k10');
    rois_dir = fullfile(char(main_dir),'rois_atlas')
    
    spm_name = fullfile(char(model_dir),'SPM.mat');
    if isempty(spm_select('list',rois_dir, par.conname))
        % a skip variable to create and display at the end of the process
        % the same if condition to be made in the pct section
        
        %sprintf('%s unexisting for this parameters',par.conname)
    else
        roi_files = fullfile(rois_dir, spm_select('list',rois_dir, par.conname));
        R = maroi(roi_files);

    %% #3b THE ROIS LOOP - creating ROIs from contrasts

    %%
    %% version to build the ROIs from image by binarizing the spmT/F image.nii
    %V = spm_vol(roi_files);
    %roi_name = char(addsuffixtofilenames(par.conname, '_ROI'));
    %R = maroi_image(struct('vol',V, 'binarize', 0, 'func', 'img'));
    %R = maroi_matrix(R);
    %saveroi(R, roi_name);

    %R = maroi(roi_files)

    %% #4 Individual stat design part 

        % set design from file /use for the individual loop
        %spm_names = spm_select([1 Inf], 'SPM.mat', ''); useful only if SPM.mat undefined or not found

        %stat_dir = get_subdir_regex(main_dir, '.*PARKGAME.*1_[a,c]$');  % regular expression for group-session folder to choose
        stat_dir = get_subdir_regex(main_dir, '.*PARKGAME.*_[a,c]$');  % regular expression for group-session folder to choose
        model_dir = get_subdir_regex(stat_dir, 'smodel_ts_tapas'); % regular expression for model folder to choose

        %stat_dir = get_subdir_regex(main_dir, '.*PARKGAME.*1_a$');
        %model_dir = get_subdir_regex(stat_dir, 'model_tedana');

        spm_names = fullfile(model_dir, spm_select('list', model_dir, 'SPM.mat'));

        %% #5 SUBJECTS LOOP - giving us individual statistics

        for subj = 1 : length(spm_names)
        %subj = 1;
            D = mardo(spm_names{subj});
            xCon = get_contrasts(D);

            Y = get_marsy(R{:}, D, 'mean');
            E = estimate(D, Y);
            E = set_contrasts(E, xCon);
            b = betas(E);
            %marsS = compute_contrasts (E, 1:length(xCon));
            [rep_strs, marsS, marsD, changef] = stat_table(E, 1:length(xCon));

            %% #6 TAB CREATION - Creating a file and writing stats in it : initialisation

            subj_names = get_parent_path(model_dir,1);
            subj_name = subj_names{subj}(length(subj_names{subj})-20 :(length(subj_names{subj})));
            subj_group = subj_name(length(subj_name));
            subj_session = subj_name(length(subj_name)-2);
            %subj_session = subj_name(length(subj_name));

            %% prepare the outf structure for output csv file
            %clear outf
            outf.suj = {subj_name};
            outf.group = subj_group;
            outf.session = subj_session;

            % ROI loop in the TAB

            for ic = 8:11
                outf.contrast = marsS.rows{ic}.name;
                for iroi = 1:length(marsS.columns)
                    nrow = nrow+1;
                    out_tab{nrow,1} = subj_name(1:5);
                    out_tab{nrow,2} = subj_group;
                    out_tab{nrow,3} = subj_session;
                    out_tab{nrow,4} = marsS.columns{iroi};
                    out_tab{nrow,5} = marsS.rows{ic}.name; % one in xCon(8:11).name range
                    out_tab{nrow,6} = marsS.con(ic,iroi);
                    out_tab{nrow,7} = marsS.stat(ic,iroi);
                    out_tab{nrow,8} = marsS.P(ic,iroi);
                    out_tab{nrow,9} = marsS.Pc(ic,iroi);
    %                 res(nrow,1) = subj_name;
    %                 res(nrow,2) = subj_group;
    %                 res(nrow,3) = subj_session;
    %                 res(nrow,4) = marsS.columns{iroi};
    %                 res(nrow,5) = marsS.rows{ic}.name;
    %                 res(nrow,6) = marsS.con(ic,iroi);
    %                 res(nrow,7) = marsS.stat(ic,iroi);
    %                 res(nrow,8) = marsS.P(ic,iroi);
    %                 res(nrow,9) = marsS.Pc(ic,iroi);
    %                 outf.roi = marsS.columns{iroi};
    %                 outf.value = marsS.con(ic,iroi);
    %                 outf.T = marsS.stat(ic,iroi);
    %                 outf.pval = marsS.P(ic,iroi);
    %                 outf.pvalC = marsS.Pc(ic,iroi);
    %                 
                    %write_result_to_csv(outf,'test3.csv', {'group'; 'roix'; 'contrastx'; 'value'; 'T'; 'pval'; 'pvalC'}) % exceeding dimensions
                    %write_result_to_csv(outf,'test3.csv')
                end
            end
            %write_result_to_csv(outf,'test3.csv', {'group'; 'roix'; 'contrastx'; 'value'; 'T'; 'pval'; 'pvalC'}) % exceeding dimensions

            %% save to txt file
            %cd (roi_model_dir)

            %fic_name = strcat(subj_name, '_', con_list.name{c}, '_ROI_stats.txt');
            %fid = fopen(fic_name,'a+');
            %for sno = 1:numel(rep_strs)
            %    fprintf(fid,'%s\n', rep_strs{sno});
            %end
            %fclose(fid);
    
    %% #8 PERCENT SIGNAL ESTIMATE TABLE :

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%here we could call a function or write a code to make a symbolic link of this file in secondlevel_ACTIVATION dir with subj_name as prefix %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Correction needed for the 'no-surviving-voxels' case
            [e_specs, e_names] = event_specs(E);
            n_events = size(e_specs, 2);
            dur = 0;
            pct_ev = cell(n_events,1);
            pctrow = pctrow+1;
            % Return percent signal estimate for all events in design
            for e_s = 1:n_events % with respectively Rest, LR, LI, RR, RI, Instruction
                pct_ev{e_s} = event_signal(E, e_specs(:,e_s), dur);
                pct_tab{1,e_s+4} = e_names{e_s};
                for t = 1:length(pct_ev{e_s}) %% CAREFUL HERE !!! CRUSH IF SOME CONS DON'T HAVE SURVIVING VOXELS
                    pct_tab{pctrow,1} = subj_name(1:5);
                    pct_tab{pctrow,2} = subj_group;
                    pct_tab{pctrow,3} = subj_session;
                    pct_tab{pctrow,4} = marsS.columns{t}; % why did it work before with the e_s as iterator ?! it definitely should not
                    pct_tab{pctrow,e_s+4} = pct_ev{e_s}(t);
                end
            end
        end
    end
end
cd (roi_model_dir)

%% #9 SAVE TABLE TO TXT : creating the out_tab to txt

[r,l] = size(out_tab);
nfid = fopen('stat_tab_PARK_S1_rois_aal.txt','w'); % nope either = format inaccessible
formatSpec = '%s;%s;%s;%s;%s;%f;%f;%f;%f\n';

fprintf(nfid,'%s;%s;%s;%s;%s;%s;%s;%s;%s\n',out_tab{1,:});
for k= 2:(r-1)
    fprintf(nfid,formatSpec,out_tab{k,1:end});
end
fprintf(nfid,'%s;%s;%s;%s;%s;%f;%f;%f;%f',out_tab{r,1:end});

fclose(nfid);

% %TEST %% Convert cell to a table and use first row as variable names
% res = [];
% for ncol = 1:length(var_tab)
%     res(nrow,ncol) = out_tab{nrow,ncol};
% end
% 
% T = cell2table(res{2:end,:},'VariableNames',res{1,:}) % nope, must
% %be a 2-D cell array
% %% Write the table to a CSV file
% %writetable(T,'test_tab.csv')

%write_result_to_csv(outf, 'stats_table.csv')

%% #10 SAVE SIGNAL PERCENT ESTIMATE TABLE TO TXT

%%SAME CASE US ABOVE @ #8 : if corrected than uncomment

%% save pct_tab to txt
[r,l] = size(pct_tab);
pctfid = fopen('signal_per_event_tab_PARK_S1_rois_aal.txt','w'); % nope either = format inaccessible
formatSpec = '%s;%s;%s;%s;%f;%f;%f;%f;%f;%f\n';

fprintf(nfid,'%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n',pct_tab{1,:});
for k= 2:(r-1)
    fprintf(pctfid,formatSpec,pct_tab{k,1:end});
end
fprintf(pctfid,'%s;%s;%s;%s;%f;%f;%f;%f;%f;%f',pct_tab{r,1:end});

fclose(pctfid);
