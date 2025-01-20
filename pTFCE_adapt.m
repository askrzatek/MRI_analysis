%% Script for automatic pTFCE contrast correction in SPM fmri data adapted from pTFCE toolbox and script
% A.SKRZATEK
% January 2022
%
%% INPUTS explanation and examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPMpath   : char variable indicating the path to SPM.mat file
%           SPMpath = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/full_resliced_multiple_regression/full_RS_clinic/AXIAL/Caudate_L/SPM.mat'
% con       : char variable indicating the path to contrast file spmT_000#.nii in SPM.mat directory
%           con = '/network/iss/cenir/analyse/irm/users/anna.skrzatek/nifti_test/full_resliced_multiple_regression/full_RS_clinic/AXIAL/Caudate_L/spmT_0001.nii'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pTFCE_Z, pTFCE_p] = pTFCE_adapt(SPMpath, con)

addpath /home/anna.skrzatek/matvol
addpath /home/anna.skrzatek/MRI_analysis

load(SPMpath)
rD = SPM.xVol.R(4);
V = SPM.xVol.S;

%get original contrast volume
V1 = spm_vol(con);
vol = spm_read_vols(V1);

parts = strsplit(SPMpath, '/');
SPMparent = sprintf('/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/',parts{2:end-1});
mask = fullfile(SPMparent,'mask.nii');

imgZ = niftiread(con);
[pTFCE_Z, pTFCE_p] = pTFCE(imgZ, mask, rD, V);

%V1.private.dat

[p nm e v] = spm_fileparts(V1.fname);

% vol1 = pTFCE_Z.*vol;
% V1.fname = [p filesep 'multivolZpTFCE_' nm e];
% spm_write_vol(V1,vol1);
% 
% y1 = vol+vol1;
% V1.fname = [p filesep 'y1_ZpTFCE_' nm e];
% spm_write_vol(V1,y1);
% 
% vol2 = vol + pTFCE_Z;
% V1.fname = [p filesep 'addvolZpTFCE_' nm e];
% spm_write_vol(V1,vol2);
% 
% vol3 = vol - pTFCE_Z;
% V1.fname = [p filesep 'diffvolZpTFCE_' nm e];
% spm_write_vol(V1,vol3);

V1.fname = [p filesep 'ZpTFCE_' nm e];
spm_write_vol(V1,pTFCE_Z);

%% p-value pTFCE file // right output image comes from its transformation
% I only need to find how the maths should sum up to get it
% I'm looking for the formula to get the right intensity - clusters are
% correct

V1.fname = [p filesep 'pval_pTFCE_' nm e];
spm_write_vol(V1,pTFCE_p);

ipval = vol.*pTFCE_p;
V1.fname = [p filesep 'multivolpval_pTFCE_' nm e];
spm_write_vol(V1,ipval);

y = vol-ipval;
V1.fname = [p filesep 'y_pTFCE_' nm e];
spm_write_vol(V1,y); % we miss .10 ^14 

% 
% ipval100 = vol.*(pTFCE_p*100);
% V1.fname = [p filesep 'multi100volpval_pTFCE_' nm e];
% spm_write_vol(V1,ipval100);
% 
% y100 = vol-ipval100;
% V1.fname = [p filesep 'y100_pTFCE_' nm e];
% spm_write_vol(V1,y100);

end