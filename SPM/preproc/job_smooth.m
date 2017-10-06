function  jobs = job_smooth(img,par)
% JOB_SMOOTH - SPM:Spatial:Smooth
%
% To build the image list easily, use get_subdir_regex & get_subdir_regex_files
%
% See also get_subdir_regex get_subdir_regex_files


%% Check input arguments

if ~exist('par','var')
    par = ''; % for defpar
end


%% defpar

defpar.smooth   = [8 8 8];
defpar.prefix   = 's';

defpar.sge      = 0;
defpar.jobname  = 'spm_smooth';
defpar.walltime = '00:30:00';

defpar.redo     = 0;
defpar.run      = 0;
defpar.display  = 0;

par = complet_struct(par,defpar);


%% SPM:Spatial:Smooth

skip = [];

for subj = 1 : length(img)
    % Test if exist
    folast = addprefixtofilenames(cellstr(char(img(subj))),par.prefix);
    if ~par.redo   &&  exist(folast{end},'file')
        skip = [skip subj];
        fprintf('[%s]: skiping subj %d because %s exist \n',mfilename,subj,folast{1});
    end
    jobs{subj}.spm.spatial.smooth.data = cellstr(char(img(subj))); %#ok<*AGROW>
    jobs{subj}.spm.spatial.smooth.fwhm = par.smooth;
    jobs{subj}.spm.spatial.smooth.dtype = 0;
    jobs{subj}.spm.spatial.smooth.prefix = par.prefix;
end


%% Other routines

[ jobs ] = job_ending_rountines( jobs, skip, par );


end % function
