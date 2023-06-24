%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biascor test, no model run
% Bill Campbell
% Last modified 6/12/2023

%% 1) General parameters
infile = fullfile('D:','Lorenz_96_model','Nature_Runs',...
    'L05M3_N960_K32_F15.00_I12_b10.00_c2.50_tf0.05_spinup1000_tsteps11000_seed5129331.mat');
outfolder = fullfile('D:','Lorenz_96_model','biascor_tests');

%% 2) Nature run
% Load nature run discarding spinup
allvars = load(infile); % full nature run
Xt = allvars.Xt.';
spinup = 1000;
Ncycles = 500;
Xt = Xt(:, spinup+1:spinup+Ncycles); % discard nature run spinup

%% 3) Observation parameters
% obs fields: seed,first,skip,err_true,err_assumed,bias,biasfac
obs = get_obs_parms();

%% 4) Forward model and observations drawn from file or constructed
Nx = size(Xt,1);
[H, R, Rhat, obs_locs_post] = get_forward_model(Nx, obs);
% Note that the true oberr (R) is used for the obs, and Rhat will be
% used in the DA step

%% 5) Draw obs and apply current bias
yt = observe_truth(Xt, H, R, obs.seed);  % Nobs x Ncycles
y0 = apply_bias(yt, obs);

%% 6) Get and smooth statistical bias corrections from previous run
% biascor: OmHB_prior,fname_OmHB_prior
biascor.obs = true;
biascor.obs_locs_post = obs_locs_post;
biascor = get_bias_corrections(Nx, outfolder, biascor);
% Apply smoothers to innovations used for bias correction
biascor = apply_sgolay(biascor);

%% 7) Parameters diagnostic print and save
fprintf('Parms values:\n');
display(obs);
display(biascor);
parmfile = fullfile(outfolder,'biascor_test_parms.mat');
save(parmfile,'obs','biascor');

%% 8) Fake DA run
biascor = fakeDA(y0,Xt,H,biascor,outfolder);

%% Fake DA routines
function biascor = fakeDA(y0,Xt,H,biascor,outfolder)
    y = fakeKSQ(y0,biascor);
    % Time mean ensemble means of differences
    OmHB = fake_diff_stats(y,Xt,H);
    % Always save the following stats
    statsfile = fullfile(outfolder,'biascor_test_output_diffstats.mat');
    save(statsfile,'OmHB','-v7.3');
end
function y = fakeKSQ(y0, biascor)
    % Apply obs bias correction directly to obs
    if biascor.obs  % Subtract mean of O - H*B from another run
        [~,K] = size(y0);
        ycor = repmat(biascor.OmHB_sg_post,1,K);
        y = y0 - ycor;
    end
end
function OmHB = fake_diff_stats(y,Xt,H)
    if ~isequal(size(y), size(H*Xt))
        error("y and H*Xt must have same size, aborting...");
    end
    Xtbar = squeeze(mean(Xt,2)); % Nx x 1
    OmHB = mean(y - H*Xtbar,2); % Nobs x 1
end

%% Observation parameters input
function obs = get_obs_parms()
    name='Obs Parameters';
    numlines=[1 80];
    opts='on';
    prompt={'Seed','First var observed','Grid skip',...
            'True Ob error','Assumed Ob error','Ob bias','Ob bias factor'};
    default={'4249687','1','1',...
             '1.0','1.0','0.3','1.0'};
    answer=inputdlg(prompt,name,numlines,default,opts);i=1;
    % Seed for rng to draw obs
    obs.seed = str2double(answer{i});i=i+1;
    % First observed variable
    obs.first = str2double(answer{i});i=i+1;
    % Variables skipped (unobserved)
    obs.skip = str2double(answer{i});i=i+1;
    % observation error variance
    obs.err_true = str2double(answer{i});i=i+1;
    obs.err_assumed = str2double(answer{i});i=i+1;
    % observation additive bias
    obs.bias = str2double(answer{i});i=i+1;
    % observation multiplicative bias factor
    obs.biasfac = str2double(answer{i});
end

%% Forward model parameters
function [H, R, Rhat, obs_locs] = get_forward_model(Nx, obs)
% Construct forward operator H, true diagonal obs error covariance R,
% and assumed obs error covariance Rhat
    oblist = zeros(1,Nx);
    rlist_true = oblist;
    rlist_assumed = oblist;
    % Observed gridpoints
    gmask = obs.first:obs.skip:Nx;
    oblist(gmask) = 1;
    rlist_true(gmask) = obs.err_true;
    rlist_assumed(gmask) = obs.err_assumed;
    % Set up observation operator and ob error covariance matrix
    % H will be Nobs x Nx, R and Rhat will be Nobs x Nobs, Nobs <= Nx
    [~, H, R, obs_locs] = forward(oblist, rlist_true);
    [~, ~, Rhat] = forward(oblist, rlist_assumed);
end

%% Observe the true state
function y0 = observe_truth(Xt, H, R, seed)
    % Create simulated obs with error covariance R and forward model H
    Rsqrt = chol(R, 'lower'); % R need not be diagonal
    Ncycles = size(Xt,2);
    Nobs = size(R,1);
    rng(seed,'twister');
    y0 = H * Xt + Rsqrt * randn(Nobs, Ncycles);
end

%% Bias the true obs
function ybias = apply_bias(yt, obs)
    % Multiplicative, then additive
    ybias = yt .* obs.biasfac;
    ybias = ybias + obs.bias;
end

%% Bias corrections input
function biascor = get_bias_corrections(Nx, outfolder,biascor)
    if biascor.obs
%         fname = fullfile(outfolder,...
%             'diffstats_alpha_1.640_tf_0.05_R_1.00_loc_n10_nc_1000_1st_1_skip_1_Nx_960_Kp_32_seed_29418_GKSQ.mat');
        fname = fullfile(outfolder,'skip5_with_bc',...
            'biascor_test_output_diffstats.mat');
        biascor.fname_OmHB_prior = fname;
        % Load bias corrections for obs
        allvars = load(fname);
        biascor = interpolate(Nx, allvars, biascor);
    end
end

%% Get OmHB_prior, OmB_prior, and obs_locs
function biascor = interpolate(Nx, allvars, biascor)
    % First look for full-sized OmB
    ombprior = false;
    if isfield(allvars, 'OmB')
        ombprior = true;
        biascor.OmB_prior = allvars.OmB(:);
        Nx_prior = length(biascor.OmB_prior);
        if Nx_prior ~= Nx
            error('OmB length is %d, should be %d, aborting.\n',...
                  Nx_prior,Nx);
        end
    end
    % Next look for OmHB_prior
    if isfield(allvars, 'OmHB_prior')
        biascor.OmHB_prior = allvars.OmHB_prior(:);
    elseif isfield(allvars,'OmHB')
        biascor.OmHB_prior = allvars.OmHB(:);
    elseif ombprior
        biascor.OmHB_prior = biascor.OmB_prior(:);
    else  % Default to zero bias correction
        print('None of OmHB_prior, OmHB, and OmB found, defaulting to 0.\n')
        biascor.OmHB_prior = zeros(Nx,1);
        biascor.OmB_prior = zeros(Nx,1);
        biascor.obs_locs_prior = (1:Nx);
    end
    Nobs = length(biascor.OmHB_prior);

    % Find index of obs locations (stored or inferred)
    if isfield(allvars, 'obs_locs')
        Nlocs = length(allvars.obs_locs);
        if Nobs ~= Nlocs
            error('OmHB_prior length=%d incompatible with obs_locs length=%d, aborting.\n',...
                   Nobs,Nlocs);
        end
        biascor.obs_locs_prior = allvars.obs_locs(:)';
    else  % Nobs is length(allvars.OmHB)
        if Nobs < Nx
            first = biascor.obs_locs_post(1);
            biascor.obs_locs_prior = infer_obs_locs(Nobs, Nx, first);
        else
            biascor.obs_locs_prior = (1:Nx);
        end
        Nlocs = length(biascor.obs_locs_prior);
    end
    % Now get biascor.OmB using obs_locs
    if Nlocs ~= Nx
        % Periodic cubic spline fit on the obs_locs_prior grid
        % csape requires Curve Fitting Toolbox
        pp = csape(biascor.obs_locs_prior, biascor.OmHB_prior, 'periodic');
        % Evaluate spline on full grid
        biascor.OmB = fnval(pp, (1:Nx)');
    else
        biascor.OmB = biascor.OmHB_prior;
    end
end

%% Apply smoothers to bias correction mean increments, innoavtions
function biascor = apply_sgolay(biascor)
    biascor.innov_smoother = 'sg';
    % Apply periodic savitzky-golay filter
    order = 9;
    framelen = 21;
    biascor.innov_smoother_parms = [order,framelen];
    biascor.OmHB_sg_prior = circ_sg_filt(biascor.OmHB_prior, order, framelen);
    % Downscale to full grid by periodic cubic spline fit to the
    % sg-filtered state on the coarse grid
    pp = csape(biascor.obs_locs_prior, biascor.OmHB_sg_prior, 'periodic');
    % Evaluate spline on full grid
    biascor.OmHB_sg_post = fnval(pp, biascor.obs_locs_post');
end
