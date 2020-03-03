function [job] = job_symbolic_child_copy(dir_func, par)
%%
%%
if empty (par)
    par = '';
end

defpar.subdir = 'tedana_vtd_mle';

%% for now just a copy of job_meica_afni creation of symbolic links

    for subj = 1 : nrSubject
    
    % Extract subject name, and print it
    subjectName = get_parent_path(dir_func{subj}(1));
    
    % Echo in terminal & initialize job_subj
    fprintf('[%s]: Preparing JOB %d/%d for %s \n', mfilename, subj, nrSubject, subjectName{1});
    job_subj = sprintf('#################### [%s] JOB %d/%d for %s #################### \n', mfilename, subj, nrSubject, dir_func{subj}{1}); % initialize
    
    nrRun = length(dir_func{subj});
    
    % nrEchoAllRuns = zeros(nrRun,1);
    
    % Create the working dir
    working_dir = char(r_mkdir(subjectName,par.subdir));
    
    %-Anat
    %======================================================================
    
    % Make symbolic link of tha anat in the working directory
    assert( exist(dir_anat{subj},'dir')==7 , 'not a dir : %s', dir_anat{subj} )
    A_src = cellstr(char(get_subdir_regex_files( dir_anat{subj}, par.anat_file_reg, struct('verbose',0))));
    A_src = char(A_src{1}); % keep the first volume, with the shorter name
    
    job_subj = [job_subj sprintf('### Anat @ %s \n', dir_anat{subj}) ]; %#ok<*AGROW>
    
    % File extension ?
    if strcmp(A_src(end-6:end),'.nii.gz')
        ext_anat = '.nii.gz';
    elseif strcmp(A_src(end-3:end),'.nii')
        ext_anat = '.nii';
    else
        error('WTF ? supported files are .nii and .nii.gz')
    end
    
    [~, anat_name, ~] = fileparts(A_src(1:end-length(ext_anat))); % remove extension to parse the file name
    anat_filename = sprintf('%s%s',anat_name,ext_anat);
    
    A_dst = fullfile(working_dir,anat_filename);
    [ ~ , job_tmp ] = r_movefile(A_src, A_dst, 'linkn', par);
    job_subj = [job_subj char(job_tmp) sprintf('\n')];
    
    ext_anat = '.nii.gz'; % force this : AFNI only generates this .nii.gz volumes
    
    %-All echos
    %======================================================================
    
    for run = 1 : nrRun
        
        % Check if dir exist
        run_path = dir_func{subj}{run} ;
        if isempty(run_path), continue, end % empty string
        assert( exist(run_path,'dir')==7 , 'not a dir : %s', run_path )
        fprintf('In run dir %s ', run_path);
        [~, serie_name] = get_parent_path(run_path);
        
        job_subj = [job_subj sprintf('### Run %d/%d @ %s \n', run, nrRun, dir_func{subj}{run}) ];
        
        prefix = serie_name;
        
        if par.redo
            % pass
        elseif exist(fullfile(working_dir,[prefix '_ctab.txt']),'file') == 2
            fprintf('[%s]: skiping %s because %s exist \n',mfilename,run_path,'ctab.txt')
            continue
        end
        
        % Fetch json dics
        jsons = get_subdir_regex_files(run_path,'^dic.*json',struct('verbose',0));
        assert(~isempty(jsons), 'no ^dic.*json file detected in : %s', run_path)
        
        % % Verify the number of echos
        % nrEchoAllRuns(run) = size(jsons{1},1);
        % assert( all( nrEchoAllRuns(1) == nrEchoAllRuns(run) ) , 'all dir_func does not have the same number of echos' )
        
        % Fetch all TE and reorder them
        res = get_string_from_json(cellstr(jsons{1}),'EchoTime','numeric');
        allTE = cell2mat([res{:}]);
        [sortedTE,order] = sort(allTE);
        fprintf(['TEs are : ' repmat('%g ',[1,length(allTE)]) ], allTE)
        
        % Fetch volume corrsponding to the echo
        allEchos = cell(length(order),1);
        for echo = 1 : length(order)
            if order(echo) == 1
                allEchos(echo) = get_subdir_regex_files(run_path, ['^f\d+_' serie_name '.nii'], 1);
            else
                allEchos(echo) = get_subdir_regex_files(run_path, ['^f\d+_' serie_name '_' sprintf('V%.3d',order(echo)) '.nii'], 1);
            end
        end % echo
        fprintf(['sorted as : ' repmat('%g ',[1,length(sortedTE)]) 'ms \n'], sortedTE)
        
        % Make symbolic link of the echo in the working directory
        E_src = cell(length(allEchos),1);
        E_dst = cell(length(allEchos),1);
        for echo = 1 : length(allEchos)
            
            E_src{echo} = allEchos{echo};
            
            % File extension ?
            if strcmp(E_src{echo}(end-6:end),'.nii.gz')
                ext_echo = '.nii.gz';
            elseif strcmp(E_src{echo}(end-3:end),'.nii')
                ext_echo = '.nii';
            else
                error('WTF ? supported files are .nii and .nii.gz')
            end
            
            filename = sprintf('%s_e%.3d%s',prefix,echo,ext_echo);
            
            E_dst{echo} = fullfile(working_dir,filename);
            [ ~ , job_tmp ] = r_movefile(E_src{echo}, E_dst{echo}, 'linkn', par);
            job_subj = [job_subj char(job_tmp)];
            
            E_dst{echo} = filename;
            
            ext_echo = '.nii.gz'; % force this : AFNI only generates this .nii.gz volumes
            
        end % echo
        

end