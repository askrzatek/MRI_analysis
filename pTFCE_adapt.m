function [pTFCE_Z, pTFCE_p] = pTFCE_adapt(SPMpath, con)

load(SPMpath)
rD = SPM.xVol.R(4);
V = SPM.xVol.S;
vol = spm_read_vols(V);

[p nm e v] = spm_fileparts(V.fname);

parts = strsplit(SPMpath, '/');
SPMparent = sprintf('/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/',parts{2:end-1});
mask = fullfile(SPMparent,'mask.nii');

imgZ = niftiread(con);
[pTFCE_Z, pTFCE_p] = pTFCE(imgZ, mask, rD, V);

end