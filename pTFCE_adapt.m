function [pTFCE_Z, pTFCE_p] = pTFCE_adapt(SPMpath, con)

load(SPMpath)
rD = SPM.xVol.R(4);
V = SPM.xVol.S;

parts = strsplit(SPMpath, '/');
SPMparent = sprintf('/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/',parts{2:end-1});
mask = fullfile(SPMparent,'mask.nii');

[pTFCE_Z, pTFCE_p] = pTFCE(con, mask, rD, V);

end