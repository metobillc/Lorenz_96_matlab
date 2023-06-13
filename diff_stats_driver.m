%% Driver for diff_stats.m function

%% Nature run and DA output
da_spinup = 100; % 2 weeks
[truth, prior, posterior, trueparms, H, R, yobs] =...
    get_run_data();

[AmB,AmBsmooth,OmB,OmBsmooth,OmA,OmAsmooth,AmT,AmTsmooth,BmT,BmTsmooth,OmT,OmTsmooth] = diff_stats(posterior,prior,yobs,truth,H,da_spinup);