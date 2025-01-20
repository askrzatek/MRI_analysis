%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A suite for a SUIT Preproc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
%
%% load data and dependencies
addpath /network/iss/cenir/software/irm/spm12b/
addpath /network/iss/cenir/analyse/irm/users/cecile.gallea/ASYA/asyasuit/
addpath /home/anna.skrzatek/MRI_analysis/

main_dir = '/home/anna.skrzatek/data/nifti_test/';
cd (main_dir)

% %if e.mat not existing
e_PARKGAME = exam(main_dir,'PARKGAME');
%e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
e = e_PARKGAME; %+ e_REMINARY; % (3:length(e_PARKGAME)); % choose specific
%
e.addSerie('ACTIVATION$','run_ACTIVATION',1)
e.addSerie(        'RS$','run_RS'        ,1)
e.addSerie('ACTIVATION','tedana009.*_vtd','tedana_ACTIVATION',1);
e.addSerie('RS','tedana009.*_vtd','tedana_RS',1);
e.getSerie('tedana_RS').addVolume('^ts.*nii','ts',1)

e.addSerie('t1mpr_S256_0_8iso_p2$','anat_T1',1)
e.getSerie('anat_T1').addVolume('^v.*p2.nii','s',1)
e.getSerie('anat').addVolume('^p0','p0',1)

% %otherwise
%load('e','e');

%% create cell variables with paths to your images
t1 = e.getSerie('anat').getVolume('s').path;
p0 = e.getSerie('anat').getVolume('p0').path;

%% go to spm
spm fmri
%% step 1: segment & isolate the cerebellum
suit_isolate_seg({t1,p0});

% add & get outputs to vars
e.getSerie('anat').addVolume('^c_v.*pcereb.nii','pcereb',1)
e.getSerie('anat').addVolume('^v.*seg1.nii','gray',1)
e.getSerie('anat').addVolume('^v.*seg2.nii','white',1)

%% step 2: normalize 
%structure per subject getting each tissue mask
job.subjND.gray = {e.gser('anat').gvol('gray').path};
job.subjND.white = {e.gser('anat').gvol('white').path};
job.subjND.isolation = {e.gser('anat').gvol('pcereb').path};

suit_normalize_dartel(job)
e.gser('anat').addVolume('^Affine.*seg1.mat','Af1',1)
e.gser('anat').addVolume('^u_a.*seg1.nii','ff',1)

clear job

%% step 3: reslice (??)
job.subj.affineTr = {e.gser('anat').gvol('Af').path};
job.subj.flowfield = {e.gser('anat').gvol('ff').path};
job.subj.resample = {e.gser('tedana_RS').gvol('ts').path};
job.subj.mask = {e.gser('anat').gvol('pcereb').path};

suit_reslice_dartel(job) %% V(i)=spm_vol(S.resample{i}) Subscripted assignment dimension mismatch.  112 V(i)=spm_vol(S.resample{i}); 

%%