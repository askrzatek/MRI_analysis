function conn_batch_singlesubject_test
% This example prompts the user to select one anatomical volume and one first-level SPM.mat file 
% and performs all connectivity analyses using all the default atlas ROIs
%
main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek','nifti');
cd /home/anna.skrzatek
addpath '/network/lustre/iss02/cenir/software/irm/spm12/toolbox/conn/'

e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
for i = 1: length(e_REMINARY)

    e_REMINARY(i).addSerie('RS$','run_RS',1);
    e_REMINARY(i).getSerie('run_RS').addVolume('^w','wT1c',1);
    e_REMINARY(i).addSerie('model_meica$','model',1);
    e_REMINARY(i).getSerie('model').addVolume('SPM.mat','SPM',1);

    %e_REMINARYS1 = exam(main_dir,'REMINARY_\w{2}_.*1$');
    %e_REMINARYS2 = exam(main_dir,'REMINARY_\w{2}_.*2$');
    e_REMINARY.explore

    ffunc = e_REMINARY(i).getSerie('run_RS').getVolume('wT1c');
    fspm = e_REMINARY(i).getSerie('model').getVolume('SPM');


    %% batch preprocessing for single-subject single-session data 

    % Selects functional / anatomical volumes
    STRUCTURAL_FILE=cellstr(ffunc.path);
    if isempty(STRUCTURAL_FILE),return;end
    SPM_FILE=cellstr(fspm.path);
    if isempty(SPM_FILE),return;end
    % Selects TR (seconds)
    %TR=inputdlg('Enter Repetition-Time (in seconds)','TR',1,{num2str('2')});
    %TR=str2num(TR{1});
    TR=1.6;

    %% CONN Setup
    fname = char(addsuffixtofilenames(addsuffixtofilenames('conn_singlesubject_test_',string(i)),'.mat'));
    batch.filename=fullfile(pwd,fname);
    if ~isempty(dir(batch.filename)), 
        Ransw=questdlg('conn_singlesubject_nopreprocesing01 project already exists, Overwrite?','warning','Yes','No','No');
        if ~strcmp(Ransw,'Yes'), return; end; 
    end
    batch.Setup.spmfiles=SPM_FILE;
    batch.Setup.structurals=STRUCTURAL_FILE;
    batch.Setup.nsubjects=1;
    batch.Setup.RT=TR;
    batch.Setup.rois.names={'atlas'};
    batch.Setup.rois.files{1}=fullfile(fileparts(which('conn')),'rois','atlas.nii');
    batch.Setup.isnew=1;
    batch.Setup.done=1;
    %% CONN Denoising
    batch.Denoising.filter=[0.01, 0.1];          % frequency filter (band-pass values, in Hz)
    batch.Denoising.done=1;
    %% CONN Analysis
    batch.Analysis.analysis_number=1;       % Sequential number identifying each set of independent first-level analyses
    batch.Analysis.measure=1;               % connectivity measure used {1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 'regression (bivariate)', 4 = 'regression (multivariate)';
    batch.Analysis.weight=2;                % within-condition weight used {1 = 'none', 2 = 'hrf', 3 = 'hanning';
    batch.Analysis.sources={};              % (defaults to all ROIs)
    batch.Analysis.done=1;

    conn_batch(batch);
%% CONN Display
% launches conn gui to explore results
conn
conn('load',fullfile(pwd,'conn_singlesubject_test.mat'));
conn gui_results


end