%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parLorGKSQ with prior inflation
% Using Lorenz 96 Model I, II or III
% Global Kalman Square Root filter
% Bill Campbell
% Last modified 6/24/2023

start = datestr(clock);
fprintf('Started: %s\n',start);

%% 1) General parameters
% run: expname,Ncycles,printcycle,use_obs_file,save_state,ring_movie,reserve
mainfolder = 'D:\Lorenz_96_model';
[infolder, run] = get_run_parms(mainfolder);
% Guard against overwriting
outfolder = create_outfolder(infolder, run);
run.outfolder = outfolder;

%% 2) Nature run
% Optional compatible observations from file,
% otherwise will draw obs as needed
% nature: ftruth,truepath,fobs,obspath,
%         K,F,I,b,c,timestep,spinup,seed,abstol,reltol
nature = get_nature_run_info(infolder, run.use_obs_file);
% Load nature run discarding spinup
[Xt, nature.abstol, nature.reltol] =...
    load_truth(nature.truepath, nature.ftruth);
Xt = Xt(:, nature.spinup+1:end); % discard nature run spinup

%% 3) Model parameters
% Default is the same as the nature run
% Deviation from default introduces model bias
% model: K,F,I,b,c,timestep,abstol,reltol,type
% Allow experiments with parameters that differ from the nature run parameters
model = get_model_parms(nature);

%% 4) DA parameters
% da: cycle_skip,K,ci,spinup,alpha,loctype,locstr,locrad
da = get_da_parms();
% Make output folder for this ensemble size if needed
ensout = fullfile(run.outfolder,['K',num2str(da.K,'%d')]);
[status, msg] = mkdir(ensout);

%% 5) Generate initial ensemble from nature run climatology
% Equally spaced through time (could randomize this later)
ics = floor(linspace(1,size(Xt,2),da.K+1)); % K x 1
XIC = Xt(:,ics(2:end)); % Nx x K
% Now can discard unneeded cycles
Xt = Xt(:,1:run.Ncycles);  % Nx x Ncycles

%% 6)  Observation parameters
% obs: seed,first,skip,err_true,err_assumed,bias,biasfac
obs = get_obs_parms();
if run.use_obs_file
    % What if nature.fobs has different obs params? Should they be loaded?
    obs.path = nature.obspath;
    obs.file = nature.fobs;
end

%% 7) Forward model
Nx = size(Xt,1);
[H, R, Rhat, obs_locs_post] = get_forward_model(Nx, obs);
% NB obs_locs_post is a column vector, want a row vec in structs

%% 8) Load or (draw and apply bias to) observations
if run.use_obs_file
    y = load_obs(obs.path, obs.file);  % observations consistent with nature run
    y = y(:, nature.spinup+1:end); % must discard same spinup period
    if (run.Ncycles > size(y,2))
        error('Error -- obs file not long enough: Ncycles=%d\n > len(obs)',...
              run.Ncycles,size(y,2));
    else % Only use needed portion of obs
        y = y(:,1:run.Ncycles);  % Nobs x Ncycles
    end
else
    yt = observe_truth(Xt, H, R, obs.seed); % Nx x Ncycles
    % Construct imperfect obs
    y = apply_obs_bias(yt, obs); % Nx x Ncycles
end

%% 9) Bias correction parameters and corrections from stats files
% biascor: obs,model,simobs_only,innov_smoother,innov_smoother_parms,
%         incr_smoother,incr_smoother_parms
biascor = get_bias_correction_parms();

%% Bias corrections from stats files
% biascor: AmB_raw,OmB_raw,fname_AmB_raw,fname_OmB_raw
biascor = get_bias_corrections(outfolder, biascor);

%% Apply smoothers to raw increments, innovations used for bias correction
% biascor: AmB, OmB
biascor = apply_smoothers(biascor);

%% Parameters diagnostic print and save
fprintf('Parms values:\n');
display(run);
display(nature);
display(model);
display(da);
display(obs);
display(biascor);
parmfile = fullfile(outfolder,...
                    ['K',num2str(da.K,'%d')],...
                    [run.expname,'_parms_GKSQ.mat']);
save(parmfile,'run','nature','model','da','obs','biascor');

%% Cycling DA run
% Execute up to Ncycles-1 DA cycles
parallelize(run.reserve);
tic;
[ensmean,post_stmse,post_stvarmse,bmse_pct,bvar_pct] = ...
    parDA_GKSQ(XIC,Xt,y,outfolder,H,Rhat,biascor,da,model,nature,run);
toc
fprintf('Posterior MSE %7.5f, %4.1f%% due to squared bias\n',...
    round(post_stmse,3,'significant'),round(bmse_pct,3,'significant'));
fprintf('(Posterior varMSE %7.5f, %4.1f%% due to squared bias)\n',...
     round(post_stvarmse,3,'significant'),round(bvar_pct,3,'significant'));

%% Plot time series of truth and ensemble mean
if run.ring_movie
    tskip = 4; % One ringplot per day
    plot_ensmean_truth(Xt,ensmean,tskip)
end

%% Output file listing
fprintf('Run output is in %s:\n',ensout);
dir(fullfile(ensout,[run.expname,'*']))

finish = datestr(clock);
fprintf('Started: %s\nFinished: %s\n',start,finish);


%% General run input
function [infolder, run] = get_run_parms(mainfolder)
    % Get general parameters for the run
    infolder = uigetdir(mainfolder,...
        'Choose output folder with Experiments subfolder:');
    name='Run Input';
    numlines=[1 120];
    opts='on';
    prompt={'Experiment Name','Cycles (4/dy)','Cycles per print','Use obs file',...
            'Save state','Ring movie','Reserved procs'};
    default={'Refactor_test','416','25','0',...
             '0','1','1'};
    answer=inputdlg(prompt,name,numlines,default,opts);i=1;
    % Unique experiment name
    run.expname = answer{i};i=i+1;
    % Number of cycles to run
    run.Ncycles = str2double(answer{i});i=i+1;
    % Cycles between print/plot of run
    run.printcycle = str2double(answer{i});i=i+1;
    % Use pre-drawn obs from file
    run.use_obs_file = logical(str2double(answer{i}));i=i+1;
    % Save full time history of state vector
    run.save_state = logical(str2double(answer{i}));i=i+1;
    % Display ring movie at end of run
    run.ring_movie = logical(str2double(answer{i}));i=i+1;
    % Processors to reserve for non-Matlab use
    run.reserve = str2double(answer{i});
end

%% Create output folder, warning if already exists
function outfolder = create_outfolder(infolder, run)
    outfolder = fullfile(infolder,'Experiments',run.expname);
    if exist(outfolder, 'dir')
        reply = input(sprintf('Directory %s exists, OK to overwrite? Y/N: ',...
            strrep(outfolder,'\','\\')),'s');
        if lower(reply(1:1)) ~= 'y'
            error('Aborting run to avoid potential overwrite of prior experiment %s in %s',...
                  run.expname,outfolder);
        end
    else
        success = mkdir(outfolder);
        if ~success
            error('Failed to create %s folder, aborting...\n',outfolder);
        end
    end
end

%% Nature run, optional observations compatible with nature run
function nature = get_nature_run_info(infolder, use_obs_file)
% Nature run file must have a filename with the following structure:
% L05M3_N960_K32_F15.00_I12_b10.00_c2.50_tf0.05_spinup100_tsteps10000_seed51422
    if exist(fullfile(infolder,'Nature_runs'),'dir')
        infolder = fullfile(infolder,'Nature_runs');
    end    
    fname = fullfile(infolder,'*L05*');
    [nature.ftruth, nature.truepath] =...
        uigetfile(fname,'Choose truth file:'); % truth trajectory
    nature.fobs = 'Unused';
    nature.obspath = 'Unused';
    if use_obs_file
        compatible = false;
        while ~compatible
            fname = fullfile(infolder,'obs*');
            [nature.fobs, nature.obspath] =...
                uigetfile(fname,'Choose obs file:'); % obs trajectory
            % Check compatibility of truth and obs
            compatible = new_check_obs_compatibility(nature.ftruth, nature.fobs);
            if ~compatible
                fprintf('Obs file\n%s\nincompatible with truth file\n%s\n, retrying...\n',...
                    nature.fobs, nature.ftruth);
            end
        end
    end
    % Extract Lorenz 96 Model I, II or III parameters from nature run filename
    [~,fname,~] = fileparts(nature.ftruth);
    string = strsplit(fname,'_');
    nature.K = str2double(string{3}(2:end));
    nature.F = str2double(string{4}(2:end));
    nature.I = str2double(string{5}(2:end));
    nature.b = str2double(string{6}(2:end));
    nature.c = str2double(string{7}(2:end));
    nature.timestep = str2double(string{8}(3:end));
    nature.spinup = str2double(string{9}(7:end));
    nature.seed = str2double(string{11}(5:end));
end

%% Model parameters
function model = get_model_parms(nature)
    % Allow experiments with parameters that differ from
    % the nature run parameters
    name='Model Parameters';
    numlines=[1 80];
    opts='on';
    prompt={'K(parm)','I(parm)',...
            'F(orcing)','b (damping)',...
            'c (coupling)','timestep',...
            'abstol','reltol'};
    default={num2str(nature.K),num2str(nature.I),...
        num2str(nature.F),num2str(nature.b),...
        num2str(nature.c),num2str(nature.timestep),...
        num2str(nature.abstol),num2str(nature.reltol)};
    answer=inputdlg(prompt,name,numlines,default,opts);i=1;
    model.K = str2double(answer{i});i=i+1;
    model.I = str2double(answer{i});i=i+1;
    model.F = str2double(answer{i});i=i+1;
    model.b = str2double(answer{i});i=i+1;
    model.c = str2double(answer{i});i=i+1;
    model.timestep = str2double(answer{i});i=i+1;
    model.abstol = str2double(answer{i});i=i+1;
    model.reltol = str2double(answer{i});
    model.type = detect_lorenz2005(model);
end

%% DA parameters input
function da = get_da_parms()
    name='DA Parameters';
    numlines=[1 80];
    opts='on';
    prompt={'Cycle skip','Ensemble size','Confidence level','Spinup',...
            'Prior inflation','Localization type','Localization radius'};
    default={'1','128','0.95','56',...
             '1.64','gc','160'};
    answer=inputdlg(prompt,name,numlines,default,opts);i=1;
    % DA cycle skipping (no obs assimilated)
    da.cycle_skip = str2double(answer{i});i=i+1;
    % number of ensemble members
    da.K = str2double(answer{i});i=i+1;
    % Confidence level (0.95 is typical)
    da.ci = str2double(answer{i});i=i+1;
    da.ci = abs(da.ci);
    if (da.ci > 1.0) % test for pctg rather than decimal
        da.ci = 100 * da.ci;
    end
    if (da.ci > 1.0)
        da.ci = 0.95; % default to 95%
    end
    % DA spinup
    da.spinup = str2double(answer{i});i=i+1;
    % prior inflation
    da.alpha = str2double(answer{i});i=i+1;
    % Localization type
    switch lower(answer{i})
        case {'gaspari','gasparicohn','g','gc','gaspari-cohn','g-c'}
            da.loctype = 1;
            da.locstr = 'g';
        case {'boxcar','box','flat','b'}
            da.loctype = 2;
            da.locstr = 'b';
        case {'triangular','tri','t','tent'}
            da.loctype = 3;
            da.locstr = 't';
        otherwise % No localization
            da.loctype = 0;
            da.locstr = 'n';
    end
    i=i+1;
    % Localization radius (in gridpoints)
    da.locrad = abs(round(str2double(answer{i})));
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

%% Bias correction parameters input
function biascor = get_bias_correction_parms()
    name='Bias Correction Input';
    numlines=[1 80];
    opts='on';
    prompt={'Apply BC to obs','Apply BC to model AND simobs',...
            'Apply BC ONLY to simobs'};
    default={'1','1',...
             '0'};
    answer=inputdlg(prompt,name,numlines,default,opts);i=1;
    % apply bias correction to obs
    biascor.apply_to_obs = logical(str2double(answer{i}));i=i+1;
    % apply bias correction to model (and therefore simulated obs)
    biascor.apply_to_model = logical(str2double(answer{i}));i=i+1;
    % apply bias correction ONLY to simulated obs, not model
    biascor.apply_to_simobs_only = logical(str2double(answer{i}));
end

%% Bias corrections input
function biascor = get_bias_corrections(outfolder, biascor)
    fname = fullfile(outfolder,'K*','diffstat*');
    % Load bias corrections for model and/or for simulated obs
    if (biascor.model || biascor.simobs_only)
        % Can e.g. add long-term mean of analysis minus background here
        [fn,pn] = uigetfile(fname,'Grab A-B');
        load([pn,fn],'AmB');
        biascor.AmB_raw = AmB;
        biascor.fname_AmB_raw = fullfile(pn,fn);
    end
    % Load bias corrections for obs
    if (biascor.obs)
        % Subtract from ob the mean of O-B from previous obs bias run with
        % perfect(i.e. identical model params e.g. forcing) model and no bias correction
        [fn,pn] = uigetfile(fname,'Grab O-B');
        load([pn,fn],'OmHB');
        biascor.OmB_raw = OmHB;
        biascor.fname_OmB_raw = fullfile(pn,fn);
    end
end

%% Parallelization options and initialization
function parallelize(reserve)
    NPROCS = feature('numcores');
    parprocs = max(1, NPROCS - abs(floor(reserve))); % Use all processes available
    if parprocs > 1 && isempty(gcp)
        parpool(parprocs);
    end
end

%% Plot time series of truth and ensemble mean
function plot_ensmean_truth(Xt,ensmean,tskip)
    myfilt = 'GKSQ';
    [Nx,Ncycles] = size(Xt);
    figure;
    plot(Xt(:,1),'k-');
    hold on; grid on;
    plot(ensmean(:,1),'r-');
    ylim([-15,20]);
    xlim([1,Nx]);
    xtt = floor(Nx/16);
    xticks([1 xtt:xtt:Nx]);
    legend('Truth','Ensmean');
    title(sprintf('%s time=1',myfilt));
    for t=tskip+1:tskip:Ncycles
        pause(0.5);
        clf;
        plot(Xt(:,t),'k-');
        hold on; grid on;
        plot(ensmean(:,t),'r-');
        ylim([-15,20]);
        xlim([1,Nx]);
        xtt = floor(Nx/16);
        xticks([1 xtt:xtt:Nx]);
        legend('Truth','Ensmean');
        title(sprintf('%s time=%d',myfilt,t));
    end % for t=2:Ncycles
end

%% Observe the true state
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
    % H will be Nobs x Nx, R will be Nobs x Nobs, Nobs <= Nx
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

%% Observation biases
function y = apply_obs_bias(y0, obs)
   % Multiplicative and additive obs biases
   y = y0.*obs.biasfac + obs.bias;
end

%% Apply smoothers to bias correction mean increments, innoavtions
function biascor = apply_smoothers(biascor)
    % Apply smoother to increments
    if (biascor.apply_to_model || biascor.apply_to_simobs_only)
        name='Increment Smoother';
        numlines=[1 80];
        opts='on';
        prompt={'Smoother type','Smoother parmlist',};
        default={'Savitzky-Golay','[9, 21]'};
        answer=inputdlg(prompt,name,numlines,default,opts);i=1;
        switch lower(answer{i})
            case {'savitzky-golay','savitzky_golay','savgolay','sg'}
                biascor.incr_smoother = answer{i};i=i+1;
                biascor.incr_smoother_parms = str2num(answer{i});
                % Apply periodic savitzky-golay filter
                order = biascor.incr_smoother_parms(1);
                framelen = biascor.incr_smoother_parms(2);
                biascor.AmB = circ_sg_filt(biascor.AmB_raw, order, framelen);
            otherwise
                biascor.incr_smoother = 'none';
                biascor.incr_smoother_parms = [];
                biascor.AmB = biascor.AmB_raw;
        end
    end
    if (biascor.apply_to_obs)
        % Apply smoother to innovationss
        name='Innnovation Smoother';
        numlines=[1 80];
        opts='on';
        prompt={'Smoother type','Smoother parmlist',};
        default={'Savitzky-Golay','[9, 21]'};
        answer=inputdlg(prompt,name,numlines,default,opts);i=1;
        switch lower(answer{i})
            case {'savitzky-golay','savitzky_golay','sgolay','sg'}
                biascor.innov_smoother = answer{i};i=i+1;
                biascor.innov_smoother_parms = str2num(answer{i});
                % Apply periodic savitzky-golay filter
                order = biascor.innov_smoother_parms(1);
                framelen = biascor.innov_smoother_parms(2);
                biascor.OmHB = circ_sg_filt(biascor.OmHB_raw, order, framelen);
            otherwise
                biascor.innov_smoother = 'none';
                biascor.innov_smoother_parms = [];
                biascor.OmHB = biascor.OmHB_raw;
        end
    end
end

%% Periodic Savitzky-Golay filter
function Xsg = circ_sg_filt(X, order, framelen)
    % Framelen must be odd
    framelen = abs(floor(framelen));
    if rem(framelen,2)==0
        error('Framelen must be odd for Savitzky-Golay filter');
    end
    % Order must be < framelen
    order = abs(floor(order));
    if order >= framelen
        error('Order must be less than framelen for Savitzky-Golay filter');
    end
    % Vector length must exceed framelen
    if isvector(X)
        nrow = length(X);
    else
        error('X must be a vector for this Savitzky-Golay filter');
    end
    if nrow < framelen
        error('Number of rows must be >= framelen for Savitzky-Golay filter');
    end
    row_output = isrow(X);
    % Convert to row vector
    if ~isrow(X)
        X = X.';
    end
    % Pad X to apply periodic boundary conditions
    pad = (framelen - 1)/2;
    Xpad = [X(end-pad+1:end) X X(1:pad)];
    Xpadsg = sgolayfilt(Xpad, order, framelen);
    Xsg = Xpadsg(pad+1:end-pad);
    if ~row_output
        Xsg = Xsg(:);
    end
end
