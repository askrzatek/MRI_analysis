% preproc physIO tapas
clc
clear all
%load e

main_dir = fullfile('/network/iss/cenir/analyse/irm/users/anna.skrzatek','nifti_test');
e_PARKGAME = exam(main_dir, 'PARKGAME');
%e_REMINARY = exam(main_dir, 'REMINARY');

e = e_PARKGAME; % e_PARKGAME + e_REMINARY;

e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)

e.getSerie('anat').addVolume('^s.*p2.nii$','s'  ,1)
e.getSerie('anat').addVolume('^y'     ,'y'  ,1)
e.getSerie('anat').addVolume('^p0'    ,'p0' ,1)
e.getSerie('anat').addVolume('^wp0'   ,'wp' ,1)

e.addSerie('ACTIVATION$','run_ACTIVATION',1);
e.addSerie(        'RS$','run_RS'        ,1);
e.addSerie('ACTIVATION','tedana_vtd','tedana_ACTIVATION',1)
e.addSerie('RS'        ,'tedana_vtd','tedana_RS'        ,1)

e.getSerie('tedana').addVolume('^ts.*.nii$'  ,'ts'  ,1)
e.getSerie('run').addVolume('^bet.*vtde1.nii','bet' ,1)
e.getSerie('run').addVolume('^bet.*mask.nii' ,'mask',1)


%% Coregistration of both runs

% ACTIVATION
    clear par
    par.type   = 'estimate';
    par.interp = 1;
    par.prefix = 'r';
    par.sge    = 0;
    par.redo   = 0;
    par.run    = 1;
    par.display= 0;

    ref = e.getSerie('anat').getVolume('^p0').toJob(0);
    src = e.getSerie('run_ACTIVATION').getVolume('^bet').toJob(0);
    oth = e.getSerie('tedana_ACTIVATION').getVolume('^ts').toJob(0);

        % include the skip option (?)


    
% RS
    ref_RS = e.getSerie('anat').getVolume('^p0').toJob(0);
    src_RS = e.getSerie('run_RS').getVolume('^bet').toJob(0);
    oth_RS = e.getSerie('tedana_RS').getVolume('^ts').toJob(0);

    % include the skip option

    job_coregister(char(src),char(ref),char(oth),par)
    
    job_coregister(char(src_RS),char(ref_RS),char(oth_RS),par)
    

