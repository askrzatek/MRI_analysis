clearvars;
clc;


% suj_path = '/network/lustre/iss02/cenir/analyse/irm/users/asya.ekmen/AMEDYST/IRM/analyse';

suj_path = '/network/lustre/iss02/cenir/analyse/irm/users/anna.skrzatek/nifti_test/'

suj = get_subdir_regex(suj_path);
%e = exam(suj_path,'PARKGAME');
%suj = gpath(e);

[pp, name] = get_parent_path(suj);
t1_format = '^v_.*p2.nii$';
t2_format = 'nothing';

number_of_threads = 4;
jobs = {};
par.verbose = 0;

for i = 1 : length(suj)
    
    temp = regexp(name{i}, '_', 'split');
    sujID = [temp{end-1}, '_', temp{end}];
    
    suit_path = fullfile(suj_path,name{i},'suit');
%     t1_path = fullfile(suj_path,'T1');
%     t2_path = fullfile(suj_path,'T2');
    t1_oriented = fullfile(suit_path,[sujID,'_t1_LPI.nii.gz']);     
    t1_biascor = fullfile(suit_path,[sujID,'_t1corr.nii.gz']);
    
    t1_fold = get_subdir_regex(suj{i}, 'Anat');
    if ~isempty(t1_fold)
        if ~exist(suit_path, 'dir'), mkdir(suit_path), end
%         if ~exist(t1_path, 'dir'), mkdir(t1_path), end
        
        par = '';
        par.verbose=0;
        par.keep = 'last';
        par.file_type = 'T1';
        t1_fold = keep_one_file(t1_fold, name{i}, par);
        t1_file = char(get_subdir_regex_files(t1_fold,t1_format,par));
%         r_movefile(t1_file,t1_path,'link');
        fprintf('--%s t1: ', sujID);
        pre_cmd = '';
        pre_cmd = sprintf('%s\ncp -r %s %s\n', pre_cmd, t1_file, t1_oriented);
        pre_cmd = sprintf('%s\nfslswapdim %s -x y z %s\n', pre_cmd, t1_oriented, t1_oriented);
        pre_cmd = sprintf('%s\nfslorient -forceneurological %s\n', pre_cmd, t1_oriented);
        unix(pre_cmd);
        
        t1_file = t1_oriented;
        t1_file = remove_gz(t1_file);
        spm_image('Display', t1_file);
        fprintf('==== Set origin of t1 to anterior commisure ===\n');
        pause;
       
        cmd = '';
        cmd = sprintf('%s\nmodule load ANTs FSL\n\n', cmd);
        cmd = sprintf('%s\nexport ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=%d\n', cmd, number_of_threads);
        cmd = sprintf('%s\nN4BiasFieldCorrection -d 3 -s 3 -i %s -o %s \n', cmd, t1_file, t1_biascor);
        cmd = sprintf('%s\nfslcpgeom %s %s\n', cmd, t1_file, t1_biascor);
        
        t2_fold = get_subdir_regex(suj{i}, t2_format);
        if ~isempty(t2_fold)
%             if ~exist(t2_path, 'dir'), mkdir(t2_path), end
            par = '';
            par.keep = 'last';
            par.verbose = 0;
            par.file_type = 'T2';
            t2_fold = keep_one_file(t2_fold, name{i}, par);
            t2_file = char(get_subdir_regex_files(t2_fold,'^s.*nii.gz',par));
            
            t2_oriented = fullfile(suit_path,[sujID,'_t2_LPI.nii.gz']);
            t2_reg = fullfile(suit_path,[sujID,'_t2reg']);
            
            fprintf('--%s t2: ', sujID);
            pre_cmd = '';
            pre_cmd = sprintf('%s\ncp -r %s %s\n', pre_cmd, t2_file, t2_oriented);
            pre_cmd = sprintf('%s\nfslswapdim %s -x y z %s\n', pre_cmd, t2_oriented, t2_oriented);
            pre_cmd = sprintf('%s\nfslorient -forceneurological %s\n', pre_cmd, t2_oriented);
            unix(pre_cmd);
            
            t2_file = t2_oriented;
            
            cmd = sprintf('%s\nantsRegistrationSyN.sh -d 3 -f %s -m %s -o %s -t r -n %d\n', cmd, t1_biascor, t2_file, t2_reg, number_of_threads);
            cmd = sprintf('%s\nfslcpgeom %s %sWarped.nii.gz\n', cmd, t1_biascor, t2_reg);
        end
        
        qsm_path = fullfile(suj_path,sujID, 'qsm_analysis');
        qsm_p2 = char(get_subdir_regex_files(qsm_path, 'p2.nii.gz$',par));
        qsm_reconstructed = char(get_subdir_regex_files(fullfile(qsm_path, 'qsm_nifti'), '.nii.gz$', par));
        t1_toqsm = fullfile(suit_path, [sujID, '_t12qsm']);
        qsm_reg = fullfile(suit_path, [sujID, '_qsmreg.nii.gz']);
        
        if ~isempty(qsm_p2) || ~isempty(qsm_reconstructed)
            qsm_p2_oriented = fullfile(suit_path,[sujID,'_qsm_LPI.nii.gz']);
            fprintf('--%s qsm_p2: ', sujID);
            pre_cmd = '';
            pre_cmd = sprintf('%s\ncp -r %s %s\n', pre_cmd, qsm_p2, qsm_p2_oriented);
            pre_cmd = sprintf('%s\nfslswapdim %s -x y z %s\n', pre_cmd, qsm_p2_oriented, qsm_p2_oriented);
            pre_cmd = sprintf('%s\nfslorient -forceneurological %s\n', pre_cmd, qsm_p2_oriented);
            unix(pre_cmd);
            
            qsm_reconstructed_oriented = fullfile(suit_path,[sujID,'_qsm_e5_LPI.nii.gz']);
            fprintf('--%s qsm_reconstructed :', sujID);
            pre_cmd = '';
            pre_cmd = sprintf('%s\ncp -r %s %s\n', pre_cmd, qsm_reconstructed, qsm_reconstructed_oriented);
            pre_cmd = sprintf('%s\nfslswapdim %s -x y z %s\n', pre_cmd, qsm_reconstructed_oriented, qsm_reconstructed_oriented);
            pre_cmd = sprintf('%s\nfslorient -swaporient %s\n', pre_cmd, qsm_reconstructed_oriented);
            unix(pre_cmd);
            
            
            qsm_p2 = qsm_p2_oriented;
            qsm_reconstructed = qsm_reconstructed_oriented;
            
            cmd = sprintf('%s\nantsRegistrationSyN.sh -d 3 -f %s -m %s -o %s -t r -n %d\n', cmd, qsm_p2, t1_biascor, t1_toqsm, number_of_threads);
            cmd = sprintf('%s\nantsApplyTransforms -d 3 -i %s -r %s -o %s -t [%s0GenericAffine.mat,1]\n', cmd, qsm_reconstructed, t1_biascor, qsm_reg, t1_toqsm);
            cmd = sprintf('%s\nfslcpgeom %s %s\n', cmd, t1_biascor, qsm_reg);
        end
        
        jobs{end+1} = cmd;
        
    end
    
    if i == length(suj)
        fprintf('====all done enjoy!!====\n');
    end
end


jobdir = '/network/lustre/iss02/cenir/analyse/irm/users/cecile.gallea/ASYA/jobs_cluster';
if ~exist(jobdir, 'dir'); mkdir(jobdir); end
cd(jobdir);

par.sge_queu = 'normal'; %{CENIR: normal, express,bigmen}, {CMRR: all.q, long.q, verylong.q, short.q, veryshort.q}
par.jobname = 'suit_preproc_encore';
par.mem = 8000; %use this format in icm '48G';
par.walltime = '01:00:00';
par.sge_nb_coeur = number_of_threads;

do_cmd_sge(jobs,par);


