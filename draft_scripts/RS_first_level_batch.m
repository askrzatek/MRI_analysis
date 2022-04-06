function RS_first_level_batch (par)
    clear all
    clc
    %% PART TO REMOVE AND REPLACE BY A JOB-FUNCTION : job_rs (par, fanat, ffunc)
    if ~exist('par','var')
        par = '';
    end

%% defpar
    defpar.anat_file_reg = '^w.*.nii';
    defpar.rs_file_reg = '^wms.*.nii';
    defpar.sge = 0;
    defpar.redo = 0;
    defpar.fake = 0;
    defpar.verbose = 2;
    defpar.jobname = 'job_rs';

    par = complet_struct(par, defpar);
    
%% CREATE A FOLDER conn_RS_data at the level of all subjects (/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti/)
    main_dir = fullfile('/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek','nifti');
    working_dir = char(r_mkdir(main_dir,'conn_RS_data'));
    cd /home/anna.skrzatek
    addpath '/network/lustre/iss02/cenir/software/irm/spm12/toolbox/conn/'
    %addpath '/network/lustre/iss02/cenir/software/irm/matvol/'
    
%% DEFINE FOLDERS & PATHS FOR ANAT & FUNCTIONAL FILES (part not to remove) & result to pass as arg in job_rs

%% Session 1
    e_REMINARY_S1 = exam(main_dir,'REMINARY_\w{2}_.*.1$');
    e_REMINARY_S2 = exam(main_dir,'REMINARY_\w{2}_.*.2$');

    e_REMINARY_S1.addSerie('RS$','run_RS',1);
    e_REMINARY_S1.getSerie('run_RS').addVolume('^w','wT1c',1);
    e_REMINARY_S1.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
    e_REMINARY_S1.getSerie('anat_T1').addVolume('^wms.*p2.nii','s',1)
    
    S1_dir = char(r_mkdir(working_dir,'REMINARY_S1')); % create a directory for S1 data
    
    fanat1 = e_REMINARY_S1.getSerie('anat_T1').getVolume('s');
    ffunc1 = e_REMINARY_S1.getSerie('run_RS').getVolume('wT1c');
    for subj = 1 : length(fanat1)
        job_subj = sprintf('#################### [%s] SUBJECT %d/%d for %s #################### \n', e_REMINARY_S1(subj).name, subj, length(fanat1), par.jobname); % initialize
        
        S1_subj = char(r_mkdir(S1_dir,e_REMINARY_S1(subj).name)); % create a directory for S1 data for subject "subj"
        
        S1_anat = char(r_mkdir(S1_subj,'anat')); % create a directory for S1 anat data for subject "subj"
        [ ~ , job_tmp ] = r_movefile(fanat1(subj).path, S1_anat, 'linkn', par);
        
        job_subj = [job_subj char(job_tmp) sprintf('\n')];
        
        S1_RS = char(r_mkdir(S1_subj,'RS')); % create a directory for S1 anat data for subject "subj"
        [ ~ , job_tmp ] = r_movefile(ffunc1(subj).path, S1_RS, 'linkn', par);
        
        job_subj = [job_subj char(job_tmp) sprintf('\n')];
        
        job {subj} = job_subj;
    end
    
    % Prepare Cluster job optimization
    if par.sge
        if par.nrCPU == 0
            par.nrCPU = 7; % on the cluster, each node have 28 cores and 128Go of RAM
        end
        par.sge_nb_coeur = par.nrCPU;
        par.mem          = 2000*(par.sge_nb_coeur+1) ;
        par.walltime     = sprintf('%0.2d',nrRun); % roughtly 1h per run, in case of slow convergeance
    end

    % Run CPU, run !
    job = do_cmd_sge(job, par);

%% Session 2
    e_REMINARY_S2.addSerie('RS$','run_RS',1);
    e_REMINARY_S2.getSerie('run_RS').addVolume('^w','wT1c',1);
    e_REMINARY_S2.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
    e_REMINARY_S2.getSerie('anat_T1').addVolume('^wms.*p2.nii','s',1)
    
    S2_dir = char(r_mkdir(working_dir,'REMINARY_S2')); % create a directory for S2 data

    fanat2 = e_REMINARY_S1.getSerie('anat_T1').getVolume('s');
    ffunc2 = e_REMINARY_S1.getSerie('run_RS').getVolume('wT1c');
    
    for subj = 1 : length(fanat2)
        job_subj = sprintf('#################### [%s] SUBJECT %d/%d for %s #################### \n', e_REMINARY_S2(subj).name, subj, length(fanat2), par.jobname); % initialize
        
        S2_subj = char(r_mkdir(S2_dir,e_REMINARY_S2(subj).name)); % create a directory for S2 data for subject "subj"
        
        S2_anat = char(r_mkdir(S2_subj,'anat')); % create a directory for S2 anat data for subject "subj"
        [ ~ , job_tmp ] = r_movefile(fanat2(subj).path, S2_anat, 'linkn', par);
        
        job_subj = [job_subj char(job_tmp) sprintf('\n')];
        
        S2_RS = char(r_mkdir(S2_subj,'RS')); % create a directory for S2 anat data for subject "subj"
        [ ~ , job_tmp ] = r_movefile(ffunc2(subj).path, S2_RS, 'linkn', par);
        
        job_subj = [job_subj char(job_tmp) sprintf('\n')];
        
        job {subj} = job_subj;
    end
    
    % Prepare Cluster job optimization
    if par.sge
        if par.nrCPU == 0
            par.nrCPU = 7; % on the cluster, each node have 28 cores and 128Go of RAM
        end
        par.sge_nb_coeur = par.nrCPU;
        par.mem          = 2000*(par.sge_nb_coeur+1) ;
        par.walltime     = sprintf('%0.2d',nrRun); % roughtly 1h per run, in case of slow convergeance
    end

    % Run CPU, run !
    job = do_cmd_sge(job, par);


    
%% END OF PART TO REPLACE BY A JOB-FUNCTION
    %% FIND functional/structural files in the newly created folder
    % note: this will look for all data in these folders, irrespestive of the specific download subsets entered as command-line arguments
    cwd=working_dir;
    
    assert( length(ffunc1) == length(fanat1), 'dir_func & dir_anat must be the same length' )
    assert( length(ffunc2) == length(fanat2), 'dir_func & dir_anat must be the same length' )
    
    %data = {'REMINARY_S1', 'REMINARY_S2'};
    session_dir=dir(fullfile(cwd,'REMINARY_S*'));
    if isempty(session_dir), return; end
    session_dir={session_dir.name};
    %session_dir=regexp(session_dir,'^REMINARY_S([12]).*','tokens','once');
    %data=unique([session_dir{:}]);
    %job = cell(data,1); %troublesome dunno what is the problem and not
    %sure whether we really need it

    fprintf('\n')
%%
    
    
    NSUBJECTS = length(fanat1);
    
    FUNCTIONAL_FILE={};
    STRUCTURAL_FILE={};
    
    for n=1:numel(session_dir),
        sessionpath = [working_dir,session_dir{n}];
        subj_dir = dir(fullfile(sessionpath,'*REMINARY_*'));
        for i=1:numel(subj_dir)
            subjpath = fullfile(sessionpath, subj_dir(i).name);
            tFUNCTIONAL_FILE=cellstr(conn_dir(fullfile(subjpath,'RS',   'wS05_RS_T1c_medn_nat.nii')));
            tSTRUCTURAL_FILE=cellstr(conn_dir(fullfile(subjpath,'anat', 'wms_S03_t1mpr_S256_0_8iso_p2.nii')));
            FUNCTIONAL_FILE=[FUNCTIONAL_FILE;tFUNCTIONAL_FILE(:)];
            STRUCTURAL_FILE=[STRUCTURAL_FILE;tSTRUCTURAL_FILE(:)];
        end
    end
    if ~NSUBJECTS, NSUBJECTS=length(STRUCTURAL_FILE); end
    if rem(length(FUNCTIONAL_FILE),NSUBJECTS),error('mismatch number of functional files %n', length(FUNCTIONAL_FILE));end
    if rem(length(STRUCTURAL_FILE),NSUBJECTS),error('mismatch number of anatomical files %n', length(FUNCTIONAL_FILE));end
    nsessions=length(FUNCTIONAL_FILE)/NSUBJECTS;
    FUNCTIONAL_FILE=reshape(FUNCTIONAL_FILE,[NSUBJECTS,nsessions]);
    STRUCTURAL_FILE={STRUCTURAL_FILE{1:NSUBJECTS}};
    disp([num2str(size(FUNCTIONAL_FILE,1)),' subjects']);
    disp([num2str(size(FUNCTIONAL_FILE,2)),' sessions']);
    TR=1.6; % Repetition time = 2 seconds
    

    %% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
    %% Prepares batch structure
    clear batch;
    batch.filename=fullfile(cwd,'conn_RS_REMINARY.mat');            % New conn_*.mat experiment name

    %% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
    % CONN Setup                                            % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options 
    % CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
    batch.Setup.isnew=1;
    batch.Setup.nsubjects=NSUBJECTS;
    batch.Setup.RT=TR;                                        % TR (seconds)
    batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);       % Point to functional volumes for each subject/session
    for nsub=1:NSUBJECTS,
        for nses=1:nsessions,
            batch.Setup.functionals{nsub}{nses}{1}=FUNCTIONAL_FILE{nsub,nses}; 
        end
    end %note: each subject's data is defined by three sessions and one single (4d) file per session
    batch.Setup.structurals=STRUCTURAL_FILE;                  % Point to anatomical volumes for each subject
    nconditions=nsessions;                                  % treats each session as a different condition (comment the following three lines and lines 84-86 below if you do not wish to analyze between-session differences)
    if nconditions==1
        batch.Setup.conditions.names={'rest'};
        for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
    else
        batch.Setup.conditions.names=[{'rest'}, arrayfun(@(n)sprintf('Session%d',n),1:nconditions,'uni',0)];
        for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
        for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=1:nsessions,  batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=[];batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=[]; end;end;end
        for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=ncond,        batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=inf;end;end;end % session-specific conditions
    end
    batch.Setup.preprocessing.steps={'functional_center', 'structural_center'};
    %batch.Setup.preprocessing.sliceorder='ascending';
    batch.Setup.done=0;
    batch.Setup.overwrite='Yes';                            

    % uncomment the following 3 lines if you prefer to run one step at a time:
%     conn_batch(batch); % runs Preprocessing and Setup steps only
%     clear batch;
%     batch.filename=fullfile(cwd,'conn_RS_REMINARY.mat');            % Existing conn_*.mat experiment name

%
    %% DENOISING step
    % CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options 
    batch.Denoising.filter=[0.01, 0.1];                 % frequency filter (band-pass values, in Hz)
    batch.Denoising.done=0;
    batch.Denoising.overwrite='Yes';

    % uncomment the following 3 lines if you prefer to run one step at a time:
    % conn_batch(batch); % runs Denoising step only
    % clear batch;
    % batch.filename=fullfile(cwd,'conn_NYU.mat');            % Existing conn_*.mat experiment name

    %% FIRST-LEVEL ANALYSIS step
    % CONN Analysis                                     % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options 
    batch.Analysis.done=1;
    batch.Analysis.overwrite='Yes';
    
    %PERFORMS FIRST-LEVEL ANALYSES (dynamic connectivity) %!
    batch.dynAnalysis.done=1;
    batch.dynAnalysis.overwrite='Yes';
    
    %PERFORMS SECOND-LEVEL ANALYSES (ROI-to-ROI and Seed-to-Voxel analyses) %!
    batch.Results.done=1;
    batch.Results.overwrite='Yes';

    %% Run all analyses
    conn_batch(batch);

    %% CONN Display
    % launches conn gui to explore results
    conn
    conn('load',fullfile(cwd,'conn_RS_REMINARY.mat'));
    conn gui_results
end
