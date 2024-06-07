function jobs = job_RS_doublerun_model_specify(dirFonc, dirOutput, par)

skip = [];
if ~exist ('par','var')
   par = '';
end
    
%% defpar
defpar.file_reg = '^s.*nii';
defpar.rp       = 0;
defpar.rp_regex = '^rp.*txt';

% Masking
defpar.mask_thr = 0.8; % spm default option
defpar.mask     =  {}; % cell(char) of the path for the mask of EACH model : N models means N paths

defpar.cvi      = 'AR(1)'; % 'AR(1)' / 'FAST' / 'none'

% Regressors
%-----------
% multilevel_cells(struct) for used defined regressors : they will NOT be convolved
defpar.user_regressor = {}; 
% multilevel_cells(char  ) for used defined regressors : they will NOT be convolved
% The regressors in the file will be concatenated with rp_*.txt

defpar.jobname  = 'spm_glm';
defpar.walltime = '04:00:00';

defpar.sge      = 0;
defpar.run      = 0;
defpar.display  = 0;
defpar.redo     = 0;
defpar.nruns    = 1;

par = complet_struct(par,defpar);

%%
n = 1;
for ip = 1 : size(dirFonc(:,1))
%     for ir = 1 : par.nruns
        jobs{n}.spm.stats.fmri_spec.dir = dirOutput(ip);
        jobs{n}.spm.stats.fmri_spec.timing.units = 'secs';
        jobs{n}.spm.stats.fmri_spec.timing.RT = 1.6;
        jobs{n}.spm.stats.fmri_spec.timing.fmri_t = 16;
        jobs{n}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        %%
%         subjectRun{ir} = get_subdir_regex_files(dirFonc{ip}{1}, par.file_reg);
%         unzip_volume(subjectRun1);
%         subjectRun{ir} = get_subdir_regex_files(dirFonc{ip}{1}, par.file_reg, struct('verbose',0));
        subjectRun = get_subdir_regex_files(dirFonc{ip}, par.file_reg, struct('verbose',0));

        clear allVolumes
        allVolumes = spm_select('expand',subjectRun);
        jobs{n}.spm.stats.fmri_spec.sess(1).scans = allVolumes; %

        %%
        jobs{n}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});

        % jobs{n}.spm.stats.fmri_spec.sess(1).multi = {'/home/anna.skrzatek/data/behav/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/20200313T141827_PARKGAMEII_P044_CK_S01_13_03_2020_MRI_run01_SPM.mat'};
        jobs{n}.spm.stats.fmri_spec.sess(1).multi = {''};
        jobs{n}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});

        % jobs{n}.spm.stats.fmri_spec.sess(1).multi_reg = {'/home/anna.skrzatek/data/nifti_test/2020_03_13_PARKGAMEII_044_CK_13_03_2020_V1_c/S12_ACTIVATION/wts_multiple_regressors.txt'};
        if par.rp && isempty(par.rp_file) && isempty(par.user_regressor)
            jobs{n}.spm.stats.fmri_spec.sess(1).multi_reg = get_subdir_regex_files(dirFonc{ip}{1},par.rp_regex);
        elseif par.rp && ~isempty(par.rp_file) && isempty(par.user_regressor)
            jobs{n}.spm.stats.fmri_spec.sess(1).multi_reg = par.rp_file{ip};
        elseif par.rp && ~isempty(par.rp_file) && ~isempty(par.user_regressor)
            jobs{n}.spm.stats.fmri_spec.sess(1).multi_reg = [par.user_regressor(ip) ; par.rp_file(ip)];
        elseif par.rp == 0 && ~isempty(par.user_regressor)
            jobs{n}.spm.stats.fmri_spec.sess(1).multi_reg = par.user_regressor(ip);
        end
%         if ~isempty(par.user_regressor)
%            jobs{n}.spm.stats.fmri_spec.sess(1).multi_reg = par.user_regressor(ip);
%         end
%         
        jobs{n}.spm.stats.fmri_spec.sess(1).hpf = 128;
        
        %%
        jobs{n}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        jobs{n}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        jobs{n}.spm.stats.fmri_spec.volt = 1;
        jobs{n}.spm.stats.fmri_spec.global = 'None';
        jobs{n}.spm.stats.fmri_spec.mthresh = par.mask_thr;
        
        if ~isempty(par.mask)
            if par.mask{ip} == ''
                jobs{n}.spm.stats.fmri_spec.mask = unzip_volume(get_subdir_regex_files(dirFonc{ip},'wbet.*mask'));
            else
                jobs{n}.spm.stats.fmri_spec.mask = unzip_volume(par.mask);
            end
        else
            jobs{n}.spm.stats.fmri_spec.mask = {''};
        end
        
        jobs{n}.spm.stats.fmri_spec.cvi = 'AR(1)';

        spm('defaults','FMRI')

    %    job = do_cmd_sge(jobs{n},par);
        n = n+1;
%     end
end

[ jobs ] = job_ending_rountines( jobs, skip, par );

    
end