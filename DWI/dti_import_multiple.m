function dti_import_multiple(dti_spm_dir,outdir,par)



if ~exist('par')
    par='';
end

defpar.skip_vol='';
defpar.sge=1;
defpar.data4D= '4D_dwi';
defpar.imgregex = '^[fs]';
defpar.swap=''; %for swap dim for sag par.swap='-y z -x' for coro par.swap='x z -y '
defpar.skip_if_exist = 1;
defpar.make_even_number_of_slice=1;
defpar.include_all = 0;
defpar.force_eddy=1;
defpar.eddy_add_cmd=' --data_is_shelled';

par = complet_struct(par,defpar);
choose_sge=par.sge;

%multiple dti subdir to import
if iscell(dti_spm_dir{1})
    for k=1:length(dti_spm_dir)
        dti_import_multiple(dti_spm_dir{k},outdir{k},par);
    end
    return
end

%case of on dti subdir but multiple subjectc
if iscell(outdir)
    for k=1:length(outdir)
        dti_import_multiple(dti_spm_dir(k),outdir{k},par);
    end
    return
end



if exist(fullfile(outdir,'bvecs'),'file') & par.skip_if_exist
    fprintf('\nsiking import because file %s  exist\n',fullfile(outdir,'bvecs'))
    return
end

dti_files=get_subdir_regex_images(dti_spm_dir,par.imgregex);
dti_files=cellstr(char(dti_files));
for k=1:length(dti_spm_dir)
    try
        bval_f(k) = get_subdir_regex_files(dti_spm_dir(k),'bvals',1);
        bvec_f(k) = get_subdir_regex_files(dti_spm_dir(k),'bvecs',1);
    catch
        fprintf('WARNING no bvals bvecs found in %s\n',dti_spm_dir{k})
        fprintf('supposing B0 acquisition\n')
        %check if dti file start with s which mean only one volume
        fftest = get_subdir_regex_files(dti_spm_dir{k},'^s.*nii')
        if isempty(fftest), error(sprintf('missing bval or bvec and no s file \n so check %s', dti_spm_dir{k})); end
        
        dd=get_parent_path(which('dti_import_multiple.m'));
        bval_f(k) = get_subdir_regex_files(dd,'bvals_b0',1);
        bvec_f(k) = get_subdir_regex_files(dd,'bvecs_b0',1);
    end
end

bval=[];bvec=[];
for k=1:length(bval_f)
    aa = load(deblank(bval_f{k}));    bb = load(deblank(bvec_f{k}));
    if par.skip_vol
        aa(par.skip_vol) = []; bb(:,par.skip_vol)=[];
    end
    bval = [bval aa];    bvec = [bvec,bb];
end
%make B0 real B0
bval(bval<50) = 0;

%copy the json, since b0 are now in the outputdir
dicjson=get_subdir_regex_files(dti_spm_dir,'^dic.*json');
if isempty(dicjson)
    %try to get in in /export/dataCENIR/dicom/nifti_raw
    dicparam = get_subdir_regex_files(dti_spm_dir,'dicom_info.mat',1);
    for nn=1:length(dicparam)
        l = load(dicparam{nn});
        [ed pd sd ] = get_description_from_dicom_hdr(l.hh);
        dirdicom = get_subdir_regex('/export/dataCENIR/dicom/nifti_raw',ed,pd,sd);
        dicjson(nn)  = get_subdir_regex_files(dirdicom,'^dic.*json',1);
    end
    
end


if iscell(outdir), outdir = char(outdir);end
if ~exist(outdir,'dir')
    mkdir(outdir);
end

dicjson = r_movefile(dicjson,outdir,'copy');

%check if eddy may be apply after topup
%acqp=topup_param_from_nifti_cenir(dti_files);
[acqp,session]=topup_param_from_json_cenir(dti_files,'',dicjson,1);
[B bi bj ]=unique(acqp,'rows');

if par.force_eddy, do_eddy=1; else do_eddy=0; end

if length(bi)==2
    if length(find(bj==1)) == length(find(bj==2)) %this mean there is the same number of AP PA acquisition so do eddy
        opo_ind=[];
        do_eddy=1;
    else %just one B0 in oposite file you do not include it in the 4D
        if length(find(bj==1)) > length(find(bj==2))
            opo_ind = find(bj==2);
        else
            opo_ind = find(bj==1);
        end
        if par.force_eddy  %in case of missing data, let's give a chance
            opo_ind=[];        do_eddy=1;
        end
    end
    
elseif length(bi)==1
    opo_ind =[];
else %there is more than 2 phase orientation
    fprintf('WARNING more than 2 phase direction for %s ',dti_spm_dir{1})
    maxdir=0;
    for kk=1:max(bj)
        if length(find(bj==kk))>maxdir, maxdir=length(find(bj==kk)); inddir=kk;end
    end
    opo_ind = find(bj~=inddir);
end

if (par.include_all) %added for case of bad acquisition for instance 1 serie AP 1 RL
    opo_ind=[]
end

ind=find(bval>50);
%should change to write in outdir with serie name
par.output_dir = outdir;
[ppp fff] = get_parent_path(dti_files);
for kkk=1:length(fff)
    fffo{kkk} = fullfile(outdir,['meanB0_' fff{kkk}]);
end
par.B0_name = fffo;
par.bval = bval_f;par.bvec = bvec_f;

fb0 = transform_4D_to_oneB0(dti_files,par);

%remove the BO in oposite phase
ind_series_toremove= unique(session(opo_ind));

dti_files(ind_series_toremove) = [];
bvec(:,opo_ind)=[];  bval(:,opo_ind)=[];


%DO MERGE
fodti=fullfile(outdir,par.data4D);
do_fsl_merge(dti_files,fodti,struct('checkorient',1));

if ~isempty(par.swap)
    fff = get_subdir_regex_files(outdir,'.*gz$');
    fff=cellstr(char(fff));
    
    do_fsl_swapdim(fff,par);
    
    vi = nifti_spm_vol(dti_files{1});
    mat=vi(1).mat(1:3,1:3);    vox = sqrt(diag(mat'*mat));      e=eye(3) ;e(1,1)=vox(1);e(2,2)=vox(2);e(3,3)=vox(3);
    rot0=mat/e;
    
    v =  nifti_spm_vol(fff{1});
    mat=v(1).mat(1:3,1:3);    vox = sqrt(diag(mat'*mat));      e=eye(3) ;e(1,1)=vox(1);e(2,2)=vox(2);e(3,3)=vox(3);
    rot=mat/e;
    bvec = inv(rot)*rot0*bvec;
    clear v vi rot mat e vox
end

%Writing bvals and bvec

fid = fopen(fullfile(outdir,'bvals'),'w');
fprintf(fid,'%d ',bval);  fprintf(fid,'\n');  fclose(fid);

fid = fopen(fullfile(outdir,'bvecs'),'w');
for kk=1:3
    fprintf(fid,'%f ',bvec(kk,:));
    fprintf(fid,'\n');
end
fclose(fid);

if par.make_even_number_of_slice
    ff=get_subdir_regex_files(outdir,'mean');
    v = nifti_spm_vol(ff{1}(1,:));
    if mod(v(1).dim(3),2)>0 %then add a slice
        ff=get_subdir_regex_files(outdir,'gz$')
        pppar.fsl_output_format ='NIFTI_GZ';    pppar.prefix = '';
        do_fsl_add_one_slice(ff,pppar);
    end
    
end

if size(B,1) == 1
    if  par.force_eddy
        %let's try eddy without topup
        fid = fopen(fullfile(outdir,'index.txt'),'w');
        fid3 = fopen(fullfile(outdir,'acqp.txt'),'w');
        fprintf(fid3,'0 1 0 0.05');
        ssess = ones(1,length(bval)); fprintf(fid,'%d ',ssess);
        fclose(fid);    fclose(fid3);

        par.sge=-1;
        [job par.mask] = do_fsl_bet({fodti},par);
        par.sge=choose_sge; par.topup='';
        do_fsl_eddy({fodti},par,job)

        
    else        
        fprintf('Sorry but there is a unique phase direction for all acquisitions, I can not do topup\n');
        topup_param_from_json_cenir(fb0,{outdir},dicjson,1); %just to write acqp.txt
        [job foDTIeddycor] = do_fsl_dtieddycor({fodti})
    end
else
    
    topup=r_mkdir({outdir},'topup');
    topup_param_from_json_cenir(fb0,topup,dicjson,1);
    fo=fullfile(topup{1},'4D_B0');
    
    fme=cellstr(char(fb0));
    for nbb=2:length(fme)
        if compare_orientation(fme(1),fme(nbb)) == 0
            fprintf('WARNING reslicing mean image because DIFFERENT ORIENTATION \n');
            bb= do_fsl_reslice( fme(nbb),fme(1));
            fme(nbb) = bb;
        end
    end
    do_fsl_merge(fme,fo,struct('checkorient',1));
    
    ffname = fullfile(topup{1},'remove_one_slice');
    if exist(ffname,'file') %then remove the last slice
        error('WARNING you need to change the number of slice should not happend !!! romain !!!\n')        
    end
    
    if do_eddy
        %then perform topup on raw B0 file (because it does the realignment and distortion correction at the same time
        
        par.sge=-1;
        job = do_fsl_topup({fo},par);
        
        par.topup_dir = topup;
        par.topup = '4D_B0_topup';
        
        %eddy
        %this write the index file for eddy each DTI refer the
        fid = fopen(fullfile(outdir,'index.txt'),'w');
        fid2 = fopen(fullfile(outdir,'session.txt'),'w');
        ll=load(fullfile(topup{1},'session.txt'));
        ind=find(bval==0);
        for k=1:length(session)
            aa  = find(ind<=k);
            fprintf(fid,'%d ',aa(end));
            fprintf(fid2,'%d ',ll(aa(end)));
            if any(acqp(k,:)-acqp(ind(aa(end)),:))
                error('ERROR bad index topup to eddy association');
            end
        end
        fclose(fid);    fclose(fid2);
        [job par.mask] = do_fsl_bet({fodti},par,job);
        par.sge=choose_sge;
        do_fsl_eddy({fodti},par,job)
        
    else %only topup so
        if par.include_all
            par.index=1;
        else
            %find out the first no opposite phase BO (first volume in 4D)
            kkk=1;
            while opo_ind(kkk)==kkk
                kkk=kkk+1;
                if kkk>length(opo_ind)
                    break
                end
            end
            %index is the first non opoB0
            par.index=kkk;
        end
        
        %so firt you need to perform eddy_correct on the 4DB0lt
        par.sge=-1; par.refnumber = kkk-1;
        [job foeddycor] = do_fsl_dtieddycor({fo},par)
        par.refnumber = 0;
        [job foDTIeddycor] = do_fsl_dtieddycor({fodti},par,job)
        
       
        job = do_fsl_topup(foeddycor,par,job);
        
        par.topup_dir = topup;
        [pp fB0name] = get_parent_path(foeddycor);
        par.topup = addsuffixtofilenames(fB0name{1},'_topup');
                        
        par.fsl_output_format = 'NIFTI_GZ';
        par.sge=choose_sge;        par.walltime='12:00:00';
        do_fsl_apply_topup(foDTIeddycor,addsuffixtofilenames(topup,['/' par.topup]),par,job)
    end
    
end

