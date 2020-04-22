% Extract Mean Percent Signal Change script.  Problem indicated with comment.
*
% MarsBar percent activation script
% Before running this do the following
% 1) make sure you have imported your ROIs from WFU and saved them
% 2) make sure you have spm and MarsBaR open and running



%%% Start marsbar to make sure spm_get works
addpath /network/lustre/iss01/cenir/software/irm/spm12/toolbox/marsbar/
marsbar('on')
% Set up the SPM defaults, just in case
spm('defaults', 'fmri');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASIC STEPS IN ROI ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spm_name = '/my/path/SPM.mat';
roi_file = '/my/path/my_roi.mat';
% Make marsbar design 
objectD = mardo(spm_name);
% Make marsbar ROI 
objectR = maroi(roi_file);
% Fetch data into marsbar data 
objectY = get_marsy(R, D, 'mean');
% Get contrasts from original 
designxCon = get_contrasts(D);
% Estimate design on ROI 
dataE = estimate(D, Y);
% Put contrasts from original design back into design 
objectE = set_contrasts(E, xCon);
% get design betasb = betas(E);
% get stats and stuff for all contrasts into statistics 
structuremarsS = compute_contrasts(E, 1:length(xCon));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti');
dirstat = fullfile(char(main_dir), 'secondlevel_ACTIVATION_PARK_S1');
% dirgroup = fullfile(char(dirstat), {'PARKGAME_a', 'PARKGAME_c'});

model_dir = fullfile(char(dirstat), 'ANOVA2x2_LxT');
%rois_dir = fullfile(char(model_dir),'rois');
rois_dir = r_mkdir(char(model_dir),'rois')
spm_name = fullfile(char(model_dir),'SPM.mat');
roi_files = fullfile(model_dir,'roi.mat'); %paths to roi files - in cellstructure

% set design from file
%spm_names = spm_select(1, 'SPM.mat', 'Select SPM.mat'); useful only if SPM.mat undefined or not found

D = mardo(spm_name);

% extract ROI data if exist
%roi_files = spm_select([1 Inf], 'mat', 'Select ROI ');

% create ROIs from contrasts

%%

R = cellstr(roi_files);
R = maroi(R);


% Fetch data into marsbar data object
mY  = get_marsy(R{1:length(roi_files)}, D, 'mean');
y = summary_data(mY);

% Get contrasts from original design
xCon = get_contrasts(D);


	Y = get_marsy(R{i}, D, 'mean')
	% Estimate design on ROI data --> each ROI separately or can we pool covariance estimate across ROIs (SPM tells its unlikely to be valid)
	E = estimate(D, Y{i});

	% Put contrasts from original design back into design object
	E = set_contrasts(E, xCon);

	% get design betas
	b = betas(E);

	% get stats and stuff for all contrasts into statistics structure
	marsS = compute_contrasts(E, 1:length(xCon));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get definitions of all events in model --> impossible because it is a second level analysis data, ergo no event time info left --> we should use 1st level anal here%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[e_specs, e_names] = event_specs(E);
	n_events = size(e_specs, 2);
	dur = 0;

% --> we still have the y - summary data though: maybe this could be useful
	
	% Return percent signal estimate for all events in design
	for e_s = 1:n_events
  		pct_ev(e_s) = event_signal(E, e_specs(:,e_s), dur);
	end
