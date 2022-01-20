function jobs = job_cat_TIV_estimate(catreport_xml,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% catreport_xml : cell string file paths to cat12 report data
%                 {'/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/2018_07_18_PARKGAMEII_001_NB_18_07_2018_V1_a/S03_t1mpr_S256_0_8iso_p2/cat_v_PARKGAMEII_001_NB_18_07_2018_V1_S3_t1mpr_S256_0_8iso_p2.xml'};
% par.fname : string output file name
%         'TIV'
% par.foutdir : string outdir path for the file to be created
%           '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    skip = [];
    if ~exist ('par','var')
        par = '';
    end
    
%% par structure
    defpar.jobname  = 'cat12_TIV_estimate';
    defpar.walltime = '04:00:00';
    defpar.run = 1;
    defpar.fname = 'TIV';
    defpar.foutdir = '/network/lustre/iss01/cenir/analyse/irm/users/anna.skrzatek/nifti_test/';

    par = complet_struct(par,defpar);

%% batch
    jobs{1}.spm.tools.cat.tools.calcvol.data_xml = catreport_xml(:);

    jobs{1}.spm.tools.cat.tools.calcvol.calcvol_TIV = 1;

    cd (par.foutdir)
    fname = sprintf('%s.txt',par.fname);
    jobs{1}.spm.tools.cat.tools.calcvol.calcvol_name = fname;

    [ jobs ] = job_ending_rountines (jobs, skip, par);
end