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

%%% Start marsbar to make sure spm_get works
addpath /network/lustre/iss01/cenir/software/irm/spm12/toolbox/marsbar/
marsbar('on')
% Set up the SPM defaults, just in case
spm('defaults', 'fmri');

%% Initialisation

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti');
roi_model_dir = fullfile(char(main_dir), 'secondlevel_ACTIVATION_PARK_S1');
roi_group.regex = {'Main_spe_LEFT_REAL_S1.*_roi.mat', 'Main_spe_LEFT_IMAGINARY_S1.*_roi.mat', 'Main_spe_RIGHT_REAL_S1.*_roi.mat', 'Main_spe_RIGHT_IMAGINARY_S1.*_roi.mat', 'Main_conj_LEFT_IMAGINARY_REAL_S1.*_roi.mat', 'Main_conj_RIGHT_IMAGINARY_REAL_S1.*_roi.mat'};
roi_group.name = {'Main_spe_LEFT_REAL_S1', 'Main_spe_LEFT_IMAGINARY_S1', 'Main_spe_RIGHT_REAL_S1', 'Main_spe_RIGHT_IMAGINARY_S1', 'Main_conj_LEFT_IMAGINARY_REAL_S1', 'Main_conj_RIGHT_IMAGINARY_REAL_S1'};
out_tab = {'id', 'group', 'session', 'roi', 'contrast', 'value', 'T', 'pval', 'pvalC'};
nrow = 1; % row index/number in out_tab - keep calm and count

% dirgroup = fullfile(char(dirstat), {'PARKGAME_a', 'PARKGAME_c'});
for roic =1 : length(roi_group.name)
    par.conname = roi_group.regex{roic};
    par.subdir = 'ANOVA2x2_LxT';

    %test_batch(par); % batch spm for loading a contrast and saving it in .mat in model dir %turns out not useful at all

    %close all

    %%
    model_dir = fullfile(char(roi_model_dir), par.subdir);
    rois_dir = get_subdir_regex(char(model_dir),'rois');

    %rois_dir = r_mkdir(char(model_dir),'rois')
    spm_name = fullfile(char(model_dir),'SPM.mat');
    roi_files = fullfile(rois_dir, spm_select('list',rois_dir, par.conname));

    R = maroi(roi_files);

    % create ROIs from contrasts

    %%
    %% version to build the ROIs from image by binarizing the spmT/F image.nii
    %V = spm_vol(roi_files);
    %roi_name = char(addsuffixtofilenames(par.conname, '_ROI'));
    %R = maroi_image(struct('vol',V, 'binarize', 0, 'func', 'img'));
    %R = maroi_matrix(R);
    %saveroi(R, roi_name);

    %R = maroi(roi_files)

    %% Individual stat design part 

    % set design from file /use for the individual loop
    %spm_names = spm_select([1 Inf], 'SPM.mat', ''); useful only if SPM.mat undefined or not found
    stat_dir = get_subdir_regex(main_dir, '.*PARKGAME.*1_\w{1}$');
    model_dir = get_subdir_regex(stat_dir, 'model_tedana');
    spm_names = fullfile(model_dir, spm_select('list', model_dir, 'SPM.mat'));

    %% subject loop giving us individual statistics
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

        %% Creating a file and writing stats in it : initialisation
        subj_names = get_parent_path(model_dir,1);
        subj_name = subj_names{subj}(length(subj_names{subj})-21 :(length(subj_names{subj})));
        subj_group = subj_name(length(subj_name));
        subj_session = subj_name(length(subj_name)-2);
        
        
        %% prepare the outf structure for output csv file
        %clear outf
        outf.suj = {subj_name};
        outf.group = subj_group;
        outf.session = subj_session;
        
        for ic = 8:11
            outf.contrast = marsS.rows{ic}.name;
            for iroi = 1:length(roi_files)
                nrow = nrow+1;
                out_tab{nrow,1} = subj_name;
                out_tab{nrow,2} = subj_group;
                out_tab{nrow,3} = subj_session;
                out_tab{nrow,4} = marsS.columns{iroi};
                out_tab{nrow,5} = marsS.rows{ic}.name; % one in xCon(8:11).name range
                out_tab{nrow,6} = marsS.con(ic,iroi);
                out_tab{nrow,7} = marsS.stat(ic,iroi);
                out_tab{nrow,8} = marsS.P(ic,iroi);
                out_tab{nrow,9} = marsS.Pc(ic,iroi);
                
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
%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%here we should call a function or write a code to make a symbolic link of this file in secondlevel_ACTIVATION dir with subj_name as prefix %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [e_specs, e_names] = event_specs(E);
        n_events = size(e_specs, 2);
        dur = 0;
        pct_ev = cell(n_events,1);

        % Return percent signal estimate for all events in design
        for e_s = 1:n_events % with respectively Rest, LR, LI, RR, RI
            pct_ev{e_s} = event_signal(E, e_specs(:,e_s), dur);
        end    
    end
end
% %% save out_tab to csv
% [r,l] = size(out_tab);
% nfid = fopen('test_tab.csv','wt'); % nope either = format inaccessible
% 
% for k=1:r
%     fprintf(out_tab,'%s, ',out_tab{k,1:end-1})
%     fprintf(out_tab,'%s\n',out_tab{k,end})
% end
% fclose(nfid);

%% Convert cell to a table and use first row as variable names
%T = cell2table(out_tab{2:end,:},'VariableNames',out_tab{1,:}) % nope, must
%be a 2-D cell array
%% Write the table to a CSV file
%writetable(T,'test_tab.csv')

%write_result_to_csv(outf, 'stats_table.csv')
        