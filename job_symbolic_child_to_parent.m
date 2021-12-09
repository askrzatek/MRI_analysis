function [job] = job_symbolic_child_to_parent(dir_func, par)
%%
% function making a symbolic link of a file contained in a subdir of the
% run dir
%%
if ~exist('par','var')
    par = ''; % for defpar
end

%% Check dir_func architecture % from job_meica_afni.m
% Transform dir_func into a multi-level cell
% Transform dir_anat into a cellstr (because only 1 anat per subj)

if ischar(dir_func)
    dir_func = cellstr(dir_func);
end

if ischar(dir_func{1})
%     if size(dir_func{1},1)>1
    if size(dir_func,1)>1 % needed for RS processing because of a different dir architecture
        dir_func = {cellstr(dir_func{1})};
    else
        dir_func = {dir_func};
    end
end

%% par structure

defpar.subdir        = 'tedana_vtd_mle';
defpar.warp_file_reg = 'wdn';

par = complet_struct(par,defpar);

%% MATLAB par structure

defpar.pct           = 0; % Parallel Computing Toolbox, will execute in parallel all the subjects
defpar.sge           = 0; % for ICM cluster, run the jobs in paralle
defpar.redo          = 0; % overwrite previous files
defpar.fake          = 0; % do everything exept running
defpar.verbose       = 2; % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all
defpar.jobname       = 'job_symbolic_child_to_parent';

par = complet_struct(par,defpar);

%% Setup that allows this scipt to prepare the commands only, no execution

parsge  = par.sge;
par.sge = -1; % only prepare commands

parverbose  = par.verbose;
par.verbose = 0; % don't print anything yet


%% for now just a copy of job_meica_afni creation of symbolic links
nrSubject = length (dir_func);
    
for subj = 1 : nrSubject

    % Extract subject name, and print it
    subjectName = get_parent_path(dir_func{subj}(1));
    
    % Echo in terminal & initialize job_subj
    fprintf('[%s]: Preparing JOB %d/%d for %s \n', mfilename, subj, nrSubject, subjectName{1});
    job_subj = sprintf('#################### [%s] JOB %d/%d for %s #################### \n', mfilename, subj, nrSubject, dir_func{subj}{1}); % initialize
    
    nrRun = length(dir_func{subj});
    
    % nrEchoAllRuns = zeros(nrRun,1);
    
    % Create the working dir
    working_dir = char(subjectName);
    tedana_dir  = char(get_subdir_regex(dir_func{subj}, par.subdir));
    
    %-Anat
    %======================================================================
    
    % Make symbolic link of tha anat in the working directory
    assert( exist(char(dir_func{subj}),'dir')==7 , 'not a dir : %s', char(dir_func{subj}) )
    A_src = cellstr(char(get_subdir_regex_files( tedana_dir, par.warp_file_reg, struct('verbose',0))));
    A_src = char(A_src{1}); % keep the first volume, with the shorter name
    
    job_subj = [job_subj sprintf('### Func @ %s \n', char(dir_func{subj})) ]; %#ok<*AGROW>
    
    % File extension ?
    if strcmp(A_src(end-6:end),'.nii.gz')
        ext_func = '.nii.gz';
    elseif strcmp(A_src(end-3:end),'.nii')
        ext_func = '.nii';
    else
        error('WTF ? supported files are .nii and .nii.gz')
    end
    
    [~, func_name, ~] = fileparts(A_src(1:end-length(ext_func))); % remove extension to parse the file name
    func_filename = sprintf('%s%s',func_name,ext_func);
    
    A_dst = fullfile(char(dir_func{subj}),func_filename);
    [ ~ , job_tmp ] = r_movefile(A_src, A_dst, 'linkn', par);
    job_subj = [job_subj char(job_tmp) sprintf('\n')];
    
    % Save job_subj in the job (job container)
    job{subj} = job_subj;
end     

% jobs are prepared

%% Remove skipable jobs (to be corrected in the future - not operational yet)

skip = false(length(job),1);
for j = 1 : length(job)
    has_ln    = ~isempty( strfind(job{j}, 'ln -sf ') );
    if ~has_ln
        skip(j) = true;
    end
end

%% Run the jobs

% Fetch origial parameters, because all jobs are prepared
par.sge     = parsge;
par.verbose = parverbose;

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

end