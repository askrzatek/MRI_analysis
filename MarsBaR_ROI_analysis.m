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

%%% Start marsbar to make sure spm_get works
addpath /network/lustre/iss01/cenir/software/irm/spm12/toolbox/marsbar/
marsbar('on')
% Set up the SPM defaults, just in case
spm('defaults', 'fmri');

%% Initialisation

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti');
roi_model_dir = fullfile(char(main_dir), 'secondlevel_ACTIVATION_PARK_S1');
con_list.regex = {'Main_spe_LEFT_REAL_S1.*_roi.mat', 'Main_spe_LEFT_IMAGINARY_S1.*_roi.mat', 'Main_spe_RIGHT_REAL_S1.*_roi.mat', 'Main_spe_RIGHT_IMAGINARY_S1.*_roi.mat', 'Main_conj_LEFT_IMAGINARY_REAL_S1.*_roi.mat', 'Main_conj_RIGHT_IMAGINARY_REAL_S1.*_roi.mat'};
con_list.name = {'Main_spe_LEFT_REAL_S1', 'Main_spe_LEFT_IMAGINARY_S1', 'Main_spe_RIGHT_REAL_S1', 'Main_spe_RIGHT_IMAGINARY_S1', 'Main_conj_LEFT_IMAGINARY_REAL_S1', 'Main_conj_RIGHT_IMAGINARY_REAL_S1'};

% dirgroup = fullfile(char(dirstat), {'PARKGAME_a', 'PARKGAME_c'});
for c =1 : length(con_list.name)
    par.conname = con_list.regex{c};
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

        %% Creating a file and writing stats in it (in subject's model_dir)
        subj_name = get_parent_path(model_dir,1);
        subj_name = subj_name{subj}(length(subj_name{subj})-21 :(length(subj_name{subj})));
        fic_name = strcat(subj_name, '_', con_list.name{c}, '_ROI_stats.txt');
        cd (roi_model_dir)
        fid = fopen(fic_name,'a+');
        for sno = 1:numel(rep_strs)
            fprintf(fid,'%s\n', rep_strs{sno});
        end
        fclose(fid);
%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%here we should call a function or write a code to make a symbolic link of this file in secondlevel_ACTIVATION dir with subj_name as prefix %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [e_specs, e_names] = event_specs(E);
        n_events = size(e_specs, 2);
        dur = 0;
        pct_ev = cell(n_events,1);

        % Return percent signal estimate for all events in design
        for e_s = 1:n_events
            pct_ev{e_s} = event_signal(E, e_specs(:,e_s), dur);
        end    
    end
end