% preproc_tedana

clear all
clc

%% Initialisation
par.pct = 0;
main_dir = fullfile('/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek','nifti');

e_PARKGAME = exam(main_dir,'PARKGAME');
%e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
%e = e_REMINARY(1:3);
%[e,ei] = removeIncomplete(e_PARKGAME); 
e = e_PARKGAME; % (3:length(e_PARKGAME)); % choose one
e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)

e.getSerie('run').addVolume('^f\d{3}','f',3)

e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^s.*p2.nii','s',1)

e.reorderSeries('path');
e.unzipVolume(par);

dir_func = e.getSerie('run') .toJob;
dir_anat = e.getSerie('anat').toJob(0);


%% quantity control
% functional ACTI
p = e.getSerie('run_ACTIVATION').getVolume('^f').print;
p = cellstr(p);
res =regexp(p,'f630');
res = cellfun(@isempty,res);
bad_volumes = char(p(res));

if ~isempty(bad_volumes)
    disp(bad_volumes)
    %error('wrong number of volumes')
else
    fprintf('all echos seem to have 630 volumes \n')
end

% functional RS
p = e.getSerie('run_RS').getVolume('^f').print;
p = cellstr(p);
res =regexp(p,'f300');
res = cellfun(@isempty,res);
bad_volumes = char(p(res));

if ~isempty(bad_volumes)
    disp(bad_volumes)
    %error('wrong number of volumes')
else
    fprintf('all echos seem to have 300 volumes \n')
end

%%
clear par

par.fname        = 'meinfo'; % name of the .mat file that will be saved
par.sge          = 0;
par.redo         = 0;
par.run          = 1;
par.fake         = 0;

run_list = e.getSerie('run');

[meinfo, jobs] = job_sort_echos( run_list , par )

meinfo.volume = e.getSerie('run').getVolume('^e')
meinfo.anat   = e.getSerie('anat').getVolume('^s')

%% AFNI
clear par
par.blocks   = {'despike','tshift','volreg'}; % now codded : despike, tshift, align, volreg
par.seperate = 1;                             % each volume is treated seperatly : useful when runs have different orientations
par.execute  = 1;                             % execute afni_proc.py generated tcsh script file immidatly after the generation

par.run      = 1;
par.sge      = 0;

par.jobname  = 'job_afni_proc_multi_echo';
par.subdir   = 'afni';

par.pct      = 0;
par.redo     = 0;
par.fake     = 0;
par.verbose  = 1;

tic
job_afni_proc_multi_echo(meinfo, par)
toc

%% get pp from afni

e.getSerie('run').addVolume('^vtde1.nii','vtde1',1)
e.getSerie('run').addVolume('^vtde2.nii','vtde2',1)
e.getSerie('run').addVolume('^vtde3.nii','vtde3',1)

[ec_vtd, ei_vtd] = e.removeIncomplete;
e = ec_vtd;

%% temporal mean

clear par
par.run   = 1;
par.fake  = 0;
par.sge   = 0;
par.fsl_output_format = 'NIFTI';

par.jobname = 'fslmean';
vtde_in  = e.getSerie('run').getVolume('^vtde1').removeEmpty.toJob;
vtde_out = addprefixtofilenames(vtde_in,'mean_');
do_fsl_mean(vtde_in,vtde_out,par);

e.getSerie('run').addVolume('^mean_vtde1.nii','mean_vtde1',1)


%% bet (brain extraction)

clear par
par.run   = 1;
par.fake  = 0;
par.sge   = 0;
par.fsl_output_format = 'NIFTI';

run_dir = e.getSerie('run').removeEmpty.toJob(0);

job = cell(0);
for iJob = 1 : length(run_dir)
    job{iJob,1} = sprintf('export FSLOUTPUTTYPE=NIFTI; cd %s; /network/lustre/iss01/cenir/software/irm/fsl5/bin/bet %s %s -R -m -f 0.3',...
        run_dir{iJob},...
        'mean_vtde1.nii',...
        'bet_mean_vtde1.nii');
end
par.jobname = 'fslbet';
do_cmd_sge(job,par);

e.getSerie('run').addVolume('^bet_mean_vtde1_mask.nii','mask_vtde1',1)
e.getSerie('run').addVolume('^bet_mean_vtde1.nii'     , 'bet_vtde1',1)


%% Tedana

% tedana.py main arguments
par.tedpca     = 'mle'; % mle, kundu, kundu-stabilize (tedana default = mle)
par.maxrestart = []; % (tedana default = 10)
par.maxit      = []; % (tedana default = 500)
par.png        = 1;  % tedana will make some PNG files of the components Beta map, for quality checks

% tedana.py other arguments
par.cmd_arg = ''; % Allows you to use all addition arguments not scripted in this job_tedana.m file

% matvol classic options
par.pct      = 0; % Parallel Computing Toolbox, will execute in parallel all the subjects
par.redo     = 0; % overwrite previous files
par.fake     = 0; % do everything exept running
par.verbose  = 2; % 0 : print nothing, 1 : print 2 first and 2 last messages, 2 : print all
par.run      = 1;

% Cluster
if par.sge               % for ICM cluster, run the jobs in paralle
    par.jobname  = 'job_tedana';
    par.walltime = '08:00:00';      % HH:MM:SS
    par.mem      = 16000;           % MB
    par.sge_queu = 'normal,bigmem'; % use both
end


job_tedana( meinfo, 'vtd', 'tedana_vtd_mle', 'bet_mean_vtde1_mask.nii', par );


%return %not sure if necessary

%% tedana report %how to use it, when meinfo is in main_dir

%tedana_report (main_dir);

%% Add new preprocessed data to exam object : dn_ts_OC.nii 

e.addSerie('ACTIVATION','tedana_vtd','tedana_ACTIVATION',1);
e.addSerie('RS','tedana_vtd','tedana_RS',1);
e.getSerie('tedana').addVolume('^dn.*.nii$','dn',1);

[ec_dn, ei_dn] = e.removeIncomplete;
e = ec_dn;

% if all in run directory instead of tedana_vtd_mle
%e(1).getSerie('run').addVolume('^dn.*.nii$','dn',1);

%% to delete if the lower part works
%% Prepare CAT12 segmentation
%%anat segment
%fanat = e.getSerie('anat').getVolume('^s').getPath;
%clear par

% global cat; cat_defaults; cat.extopts.subfolders=0; cat.extopts.expertgui=1;clear defaults; spm_jobman('initcfg');
% 
% par.run     = 1;
% par.display = 0;
% par.redo    = 1;
% par.sge     = 0;
% 
% par.jobname      = 'spm_segmentCAT12';
% par.mem          = 8000;
% par.sge_nb_coeur = 2
% 
% par.walltime  = '08:00:00';
% par.subfolder = 0;
% expert_mode   = 1;
% par.cmd_prepend = sprintf('global cat; cat_defaults; cat.extopts.subfolders=%d; cat.extopts.expertgui=%d;clear defaults; spm_jobman(''initcfg'');',...
%     par.subfolder,expert_mode);
% par.matlab_opt = ' -nodesktop ';
% 
% 
% jobs = cell(1,length(fanat));
% for j = 1 : length(fanat)
%    
%     jobs{j}.spm.tools.cat.estwrite.data = fanat(j);
%     jobs{j}.spm.tools.cat.estwrite.nproc = 0;
%     jobs{j}.spm.tools.cat.estwrite.opts.tpm = {'/network/lustre/iss01/cenir/software/irm/spm12/tpm/TPM.nii'};
%     jobs{j}.spm.tools.cat.estwrite.opts.affreg = 'mni';
%     jobs{j}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
%     jobs{j}.spm.tools.cat.estwrite.opts.samp = 3;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 0.5;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.WMHCstr = 0.5;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 1;
%     jobs{j}.spm.tools.cat.estwrite.extopts.segmentation.restypes.best = [0.5 0.3];
%     jobs{j}.spm.tools.cat.estwrite.extopts.registration.darteltpm = {'/network/lustre/iss01/cenir/software/irm/spm12/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii'};
%     jobs{j}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {'/network/lustre/iss01/cenir/software/irm/spm12/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'};
%     jobs{j}.spm.tools.cat.estwrite.extopts.registration.regstr = 0;
%     jobs{j}.spm.tools.cat.estwrite.extopts.vox = 1.5;
%     jobs{j}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
%     jobs{j}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
%     jobs{j}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
%     jobs{j}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0;
%     jobs{j}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 0;
%     jobs{j}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
%     jobs{j}.spm.tools.cat.estwrite.extopts.admin.print = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.surface = 12;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.GM.native = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.GM.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.GM.mod = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.GM.dartel = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.WM.native = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.WM.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.WM.mod = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.WM.dartel = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.CSF.native = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.CSF.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.CSF.mod = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.CSF.dartel = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.WMH.native = 0;
%     jobs{j}.spm.tools.cat.estwrite.output.WMH.warped = 0;
%     jobs{j}.spm.tools.cat.estwrite.output.WMH.mod = 0;
%     jobs{j}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
%     jobs{j}.spm.tools.cat.estwrite.output.label.native = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.label.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.label.dartel = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.bias.native = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.bias.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.bias.dartel = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.las.native = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.las.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.las.dartel = 2;
%     jobs{j}.spm.tools.cat.estwrite.output.jacobian.warped = 1;
%     jobs{j}.spm.tools.cat.estwrite.output.warps = [1 1];
%    
% end
% 
% 
% 
% job_ending_rountines(jobs,[],par)

%e.getSerie('UNI').addVolume('^mcs' ,'mcs' ,1)
%e.getSerie('UNI').addVolume('^wmcs','wmcs',1)

%% to delete if the upper part works
%% Run CAT12 segmentation
clear par
%anat segment
fanat = e.getSerie('anat').getVolume('^s').getPath;

% Retrocompatibility for SPM:Spatial:Segment options
par.GM        = [1 1 1 1]; % warped_space_Unmodulated(wp*) / warped_space_modulated(mwp*) / native_space(p*) / native_space_dartel_import(rp*)
par.WM        = [1 1 1 1];
par.CSF       = [1 1 1 1];
par.bias      = [1 1 0] ;  % native normalize dartel     [0 1]; % bias field / bias corrected image
par.warp      = [1 1]; % warp field native->template / warp field native<-template

par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;
par.doSurface = 0;
par.subfolder = 0; % all results in the same subfolder

par.run     = 1;
par.display = 0;
par.pct     = 0;
par.redo    = 0;

j_segment = job_do_tedana_segmentCAT12(fanat,par);


e.getSerie('anat').addVolume('^y'  ,'y' );
e.getSerie('anat').addVolume('^p0' ,'p0' );
%e.getSerie('anat').addVolume('^ms' ,'ms');
%e.getSerie('anat').addVolume('^wms','wms');

% %% bet @ anat_T1_UNI - mcs
%
% clear par
% par.run   = 0;
% par.fake  = 0;
% par.sge   = 1;
% par.fsl_output_format = 'NIFTI';
%
% anat_in  = gfile(e.getSerie('anat_T1_UNI').getPath,'^p0');
% anat_out = addprefixtofilenames(anat_in,'bet_');
%
% job = cell(0);
% for iJob = 1 : length(anat_in)
%     job{iJob,1} = sprintf('export FSLOUTPUTTYPE=NIFTI; /network/lustre/iss01/cenir/software/irm/fsl5/bin/bet %s %s -v -R -m',...
%         anat_in {iJob} ,...
%         anat_out{iJob} );
% end
% par.jobname = 'fslbet';
% do_cmd_sge(job,par);


%%

% e.getSerie('anat_T1_UNI').addVolume('^bet_mcs_.*_fatnav.nii'     , 'bet_mcs',1)
% e.getSerie('anat_T1_UNI').addVolume('^bet_mcs_.*_fatnav_mask.nii','mask_mcs',1)


%% Coregister func @ anat
% correction of any former transformations of ^dn

    par.ask = 0;
    par.redo = 1;
    par.pct = 0;
    
    if par.redo==1
        
        e.getSerie('tedana').addVolume('^ts','ts',1);
        
        %% for ACTIVATION vols
        origin = e.getSerie('tedana_ACTIVATION').getVolume('^ts').toJob(0);
        influencer = e.getSerie('run_ACTIVATION').getVolume('^bet_vtde1').toJob(0);
        follower = e.getSerie('tedana_ACTIVATION').getVolume('^dn').toJob(0);
        %e.getSerie('tedana_ACTIVATION').addVolume('.*.mat$','mat',1);
        %transf = e.getSerie('tedana_ACTIVATION').getVolume('mat');
        %if ~isempty(transf)
    %         for i = 1 : length(transf)
    %             if exist(transf(i).path, 'file') && par.redo==1
    %                 par.ask = 0;
    %                 do_fsl_copy_header(influencer{i}, origin{i}, par); % clear the former transformations if any
    %                 do_fsl_copy_header(follower{i}, influencer{i}, par); % clear the former transformations if any
    %             end
    %         end
            
            do_fsl_copy_header(influencer, origin,      par);
            do_fsl_copy_header(follower,   influencer,  par);
            clear par
            par.type   = 'estimate';
            par.interp = 1;
            par.prefix = 'r';
            par.sge    = 0;
            par.redo   = 0;
            par.run    = 1;
            par.display= 0;
        %end
        ref = e.getSerie('anat').getVolume('^p0').toJob(0);
        src = e.getSerie('run_ACTIVATION').getVolume('^bet_vtde1').toJob(0);
        oth = e.getSerie('tedana_ACTIVATION').getVolume('^dn').toJob(0);

        % include the skip option

        job_coregister(char(src),char(ref),char(oth),par)

    %% for RS vols
        origin = e.getSerie('tedana_RS').getVolume('^ts').toJob(0);
        influencer = e.getSerie('run_RS').getVolume('^bet_vtde1').toJob(0);
        follower = e.getSerie('tedana_RS').getVolume('^dn').toJob(0);
        %e.getSerie('tedana_ACTIVATION').addVolume('.*.mat$','mat',1);
        %transf = e.getSerie('tedana_ACTIVATION').getVolume('mat');
        %if ~isempty(transf)
    %         for i = 1 : length(transf)
    %             if exist(transf(i).path, 'file') && par.redo==1
    %                 par.ask = 0;
    %                 do_fsl_copy_header(change{i}, origin{i}, par); % clear the former transformations if any
    %             end
    %         end
            par.ask = 0;
            par.pct = 0;
            do_fsl_copy_header(influencer, origin,      par);
            do_fsl_copy_header(follower,   influencer,  par);
            clear par
            par.type   = 'estimate';
            par.interp = 1;
            par.prefix = 'r';
            par.sge    = 0;
            par.redo   = 0;
            par.run    = 1;
            par.display= 0;
        %end
        ref = e.getSerie('anat').getVolume('^p0').toJob(0);
        src = e.getSerie('run_RS').getVolume('^bet_vtde1').toJob(0);
        oth = e.getSerie('tedana_RS').getVolume('^dn').toJob(0);
        
        % include the skip option

        job_coregister(char(src),char(ref),char(oth),par)
    
    end
    

    
    %% Normalize

    clear par
    %par.preserve = 0;
    %par.bb       = [NaN NaN NaN ; NaN NaN NaN];
    %par.vox      = [2.5 2.5 2.5];
    %par.interp   = 4;
    %par.wrap     = [0 0 0];
    %par.prefix   = 'w';

    par.redo    = 1;
    par.sge     = 0;
    par.run     = 1;
    par.display = 0;
    par.jobname = 'spm_apply_norm';

    warp_field = e.getSerie('anat').getVolume('^y');
    img = e.getSerie('tedana').getVolume('^dn');

    job_apply_normalize(warp_field,img, par)

    e.explore

    save('e','e')
    
    %%

    %e.getSerie('tedana_ACTIVATION').addVolume('^wdn_ts_OC.nii','wdn_ts_OC',1)


%     %% Smooth
% 
% 
%     clear par
% 
%     par.sge      = 0;
%     par.redo     = 0;
%     par.run      = 1;
%     par.display  = 0;
% 
%     img = e.getSerie('tedana').getVolume('^wdn_ts_OC').removeEmpty;
% 
%     par.smooth   = [5 5 5];
%     par.prefix   = 's5';
%     job_smooth(img,par)
% 
%     par.smooth   = [8 8 8];
%     par.prefix   = 's8';
%     job_smooth(img,par)
% 
% end
%%end of run loop

%% apply normalization on fmri tedana volume, using SPM warp field from the segementation
% 
% % Prepare
% 
% ffunc = e.getSerie('tedana').getVolume('^dn');
% fy    = e.getSerie('anat').getVolume('^y');
% 
% par.display= 0;
% par.run    = 1;
% par.pct    = 1;
% 
% 
% %% Run normalize
% 
% par.auto_add_obj=0;
% job_apply_normalize(fy,ffunc,par);
% e.getSerie('tedana').addVolume('^w.*nii','wdn',1)
% 
%% Add the warped data to exam object : (?) wdn_ts_OC(_dn?).nii


