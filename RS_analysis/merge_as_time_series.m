function merge_as_time_series (root, foldregex, targetregex, outdir, outname)
    %root = {'/home/anna.skrzatek/data_orig/nifti/'};
    suj = get_subdir_regex(root,{'^2'});
    %foldregex = 'p2$';
    targetfold = get_subdir_regex(suj,{foldregex});
    %targetregex = '^ws';
    target = get_subdir_regex_files(targetfold,{targetregex});

    %outdir = {'/home/anna.skrzatek/data_orig/nifti'};
    cd(outdir{1});
    %outname = 'merged_t1'
    do_fsl_merge(target,outname);
end