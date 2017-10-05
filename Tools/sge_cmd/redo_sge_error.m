function job = redo_sge_error(jobdir,par)

if ~exist('par'),par ='';end

defpar.grep_error=1; %if ==0 then it will look into log file
defpar.msg='illed';
defpar.jobname='redoo';
defpar.logdir = 'bad_log';
defpar.sge=1;

par = complet_struct(par,defpar);

grep_error = par.grep_error;
error_msg = par.msg;

if ~exist('jobdir','var'), jobdir=pwd;end

if ischar(jobdir)
    jobdir={jobdir};
end

%erf = get_subdir_regex_files(jobdir,[jobname{1} '.e']);
erf = get_subdir_regex_files(jobdir,'err');
erf = cellstr(char(erf));
logf =  get_subdir_regex_files(jobdir,'log');logf = cellstr(char(logf));


redodir = r_mkdir(jobdir,par.logdir);
jobagain={};
for k=1:length(erf)
    if grep_error,    thefile=erf{k};
    else thefile = logf{k}; end
    
    s=dir(thefile );
    if s.bytes % only non empty files
        if grep_error
            cmd = sprintf('cat %s |grep "%s"',thefile,error_msg);
        else
            cmd = sprintf('cat %s |grep "%s"',thefile,error_msg);
        end
        
        [a b] = unix(cmd);
        if ~isempty(b)
            %find job_number
            [pp, efile] = fileparts(erf{k});
            [pp jobname] = fileparts(pp);
            
            ii = strfind(efile,'_');
            numjob = str2num(efile(ii+1:end));
            
            r_movefile(erf(k),redodir{1},'move');
            r_movefile(logf(k),redodir{1},'move');
            
            jobfile = sprintf('j%.2d_%s',numjob,jobname);
            fid = fopen(jobfile);
            l = fgetl(fid);cmd=''; %skip first line #bash
            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                cmd = sprintf('%s%s\n',cmd,tline);
            end
            fclose(fid);
            
            jobagain{end+1} = cmd;

        else
            if grep_error
            fprintf(' Error file %s not empty but not error message %s\n',erf{k},error_msg);
            end
        end
    end
end


do_cmd_sge(jobagain,par);

if 0 %delete freesurfer sujdir
    a=str2num(strjob)
    for k=1:length(a)
        jname=sprintf('j%.2d_freesurfer_reconall',a(k))
        l=readtext(jname);
        aa=split(l{3})
        cmd{k} = sprintf('rm -rf %s/%s',aa{7},aa{5})
        par.sge=0;
        do_cmd_sge(cmd,par)
        
    end
end

