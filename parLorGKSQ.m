%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parLorGKSQ with prior inflation
% Using Lorenz 96 Model I, II or III
% Global Kalman Square Root filter
% Bill Campbell
% Last modified 9/17/2022

%% General parameters
[outfolder, use_obs_file, run_parms] = get_run_parms;

%% Nature run
% Optional compatible observations from file, otherwise will draw obs
% as needed
nature_parms = get_nature_run_info(outfolder, use_obs_file);

%% Model parameters
% Default is the same as the nature run
model_parms = get_model_parms(nature_parms);

%% Observation parameters
obs_parms = get_obs_parms(use_obs_file);
if use_obs_file
    obs_parms.obspath = nature_parms.obspath;
    nature_parms = rmfield(nature_parms, 'obspath');
    obs_parms.fobs = nature_parms.fobs;
    nature_parms = rmfield(nature_parms, 'fob');
end

%% DA parameters
da_parms = get_da_parms(outfolder);

%% Bias correction parameters
biascor_parms = get_bias_correction_parms(outfolder);

%% Startup
parprocs = parallelize(run_parms.reserve);

%% Load nature run discarding spinup
[Xt, model_parms.abstol, model_parms.reltol] = load_truth(nature_parms); % full nature run
Xt = Xt(:,nature_parms.spinup+1:end); % discard nature run spinup
Nx = size(Xt,1);

%% Load obs file if desired, discarding spinup
if use_obs_file
    y = load_obs(obs_parms);  % observations consistent with nature run
    y = y(:,nature_parms.spinup+1:end); % must discard same spinup period
    if (run_parms.Ncycles > size(y,2))
        error('Error -- obs file not long enough: Ncycles=%d\n > len(obs)',...
              run_parms.Ncycles,size(y,2));
    else % Only use needed portion of nature run and obs
        y = y(:,1:run_parms.Ncycles);  % Nx x Ncycles
    end
else
    y = NaN(Nx, run_parms.Ncycles);
end

%% Generate initial ensemble from nature run climatology
% Equally spaced through time
ics = floor(linspace(1,size(Xt,2),da_parms.K+1)); % K x 1
XIC = Xt(:,ics(2:end)); % Nx x K
% Now can discard unneeded cycles
Xt = Xt(:,1:run_parms.Ncycles);  % Nx x Ncycles

%% Construct forward operator H and diagonal obs error covariance R
oblist = zeros(1,Nx);
rlist = oblist;
% Observed gridpoints
gmask = obs_parms.first:obs_parms.skip:Nx;
oblist(gmask)=1;
rlist(gmask)=obs_parms.oberr;
% Set up observation operator and ob error covariance matrix
% H will be Nobs x Nx, R will be Nobs x Nobs, Nobs <= Nx
[~,H,R] = forward(oblist,rlist);
% Compress the observation vector to match compressed R
y = y(gmask,:);  % Nobs x Ncycles

%% Cycling DA run
tstart=tic;
pa
% Execute up to Ncycles-1 DA cycles
[tot_err,tot_avar] = ...
    parDA_GKSQ(XIC,Xt,y,outfolder,use_obs_file,H,R,...
    da_parms,model_parms,nature_parms,obs_parms);
terr = tot_err(ialpha,:);
toc

%% Plot time series of truth and ensemble mean
if run_parms.ring_movie
    plot_ensmean_truth(outfolder,Xt,...
        da_parms,model_parms,obs_parms,run_parms)
end
finish = datestr(clock);
fprintf('Started: %s\nFinished: %s\n',start,finish);


%% General run input
function [outfolder, use_obs_file, run_parms] = get_run_parms()
    % Get general parameters for the run
    outfolder = uigetdir('','Choose output folder:');
    name='Run Input';
    numlines=1;
    prompt={'Cycles','Cycles per print','Seed','Use obs file',...
            'Save state','Ring movie','Reserved procs'};
    default={'1000','100','24918','1','1','0','1'};
    answer=inputdlg(prompt,name,numlines,default);i=1;    
    % Number of cycles to run
    run_parms.Ncycles = str2double(answer{i});i=i+1;
    % Cycles between print/plot of run
    run_parms.printcycle = str2double(answer{i});i=i+1;
    % Seed for rng in ensemble DA
    run_parms.seed = str2double(answer{i});i=i+1;
    % Use pre-drawn obs from file
    use_obs_file = logical(str2double(answer{i}));i=i+1;
    % Save full time history of state vector
    run_parms.save_state = logical(str2double(answer{i}));i=i+1;
    % Display ring movie at end of run
    run_parms.ring_movie = logical(str2double(answer{i}));i=i+1;
    % Processors to reserve for non-Matlab use
    run_parms.reserve = str2double(answer{i});
end

%% DA parameters input
function da_parms = get_da_parms(outfolder)
    name='DA Parameters';
    numlines=1;
    prompt={'Time skip','Ensemble size','Confidence level','Prior inflation',...
            'Localization type','Localization radius'};
    default={'1','1','1','500','0.95','1.64','gc','160'};
    answer=inputdlg(prompt,name,numlines,default);i=1;
    % Time steps skipped (completely unobserved)
    da_parms.tskip = str2double(answer{i});i=i+1;
    % number of ensemble members
    da_parms.K = str2double(answer{i});i=i+1;
    % Make output folder for this ensemble size
    if ~exist([outfolder,'K',num2str(da_parms.K,'%d\n')], 'dir')
        mkdir([outfolder,'K',num2str(da_parms.K,'%d\n')]);
    end
    % Confidence level (0.95 is typical)
    da_parms.ci = str2double(answer{i});i=i+1;
    % prior inflation
    da_parms.alpha = str2double(answer{i});i=i+1;
    % Localization type
    switch lower(answer{i})
        case {'gaspari','gasparicohn','g','gc','gaspari-cohn','g-c'}
            da_parms.loctype = 1;
            da_parms.locstr = 'g';
        case {'boxcar','box','flat','b'}
            da_parms.loctype = 2;
            da_parms.locstr = 'b';
        case {'triangular','tri','t','tent'}
            da_parms.loctype = 3;
            da_parms.locstr = 't';
        otherwise % No localization
            da_parms.loctype = 0;
            da_parms.locstr = 'n';
    end
    i=i+1;
    % Localization radius (in gridpoints)
    da_parms.locrad = abs(round(str2double(answer{i})));
end

%% Observation parameters input
function obs_parms = get_obs_parms()
    name='Obs Parameters';
    numlines=1;
    prompt={'First var observed','Grid skip',...
        'Ob error','Ob bias','Ob bias factor',};
    default={'1','1','1.0','0.0','1.0'};
    answer=inputdlg(prompt,name,numlines,default);i=1;
    % First observed variable
    obs_parms.first = str2double(answer{i});i=i+1;
    % Variables skipped (unobserved)
    obs_parms.skip = str2double(answer{i});i=i+1;
    % observation error variance
    obs_parms.oberr = str2double(answer{i});i=i+1;
    % observation additive bias
    obs_parms.obbias = str2double(answer{i});i=i+1;
    % observation multiplicative bias factor
    obs_parms.obbiasfac = str2double(answer{i});
end

%% Bias correction parameters input
function biascor_parms = get_bias_correction_parms(outfolder)
    name='Bias Correction Input';
    numlines=1;
    prompt={'Apply BC to obs','Apply BC to model',...
            'Apply BC ONLY to simobs'};
    default={'0','0','0'};
    answer=inputdlg(prompt,name,numlines,default);i=1;
    % apply bias correction to obs
    biascor_parms.obs = logical(str2double(answer{i}));i=i+1;
    % apply bias correction to model (and therefore simulated obs)
    biascor_parms.model = logical(str2double(answer{i}));i=i+1;
    % apply bias correction ONLY to simulated obs, not model
    biascor_parms.simobs_only = logical(str2double(answer{i}));
    
    % Load bias corrections for model and/or for simulated obs
    if (biascor_parms.model || biascor_parms.simobs_only)
        % Can e.g. add long-term mean of analysis minus background here
        [fn,pn,~] = uigetfile([outfolder,'K*/diffstat*'],'Grab A-B');
        biascor_parms.AmB = load([pn,fn],'AmB');
    end

    % Load bias corrections for obs
    if (biascor_parms.obs)
        % Subtract from ob the mean of O-B from previous obs bias run with
        % perfect(i.e. identical model params e.g. forcing) model and no bias correction
        [fn,pn,~] = uigetfile([outfolder,'K*/diffstat*'],'Grab O-B');
        biascor_parms.OmB = load([pn,fn],'OmB');
    end
end

%% Nature run, optional observations compatible with nature run
function nature_parms = get_nature_run_info(outfolder, use_obs_file)
% Nature run file must have a filename with the following structure:
% L05M3_N960_K32_F15.00_I12_b10.00_c2.50_tf0.05_spinup100_tsteps10000_seed51422
    [nature_parms.ftruth, nature_parms.truepath] =...
        uigetfile([outfolder,'*L05*'],'Choose truth file:'); % truth trajectory
    nature_parms.fobs = '';
    nature_parms.obspath = '';
    if use_obs_file
        compatible = false;
        while ~compatible
            [nature_parms.fobs, nature_parms.obspath] =...
                uigetfile([outfolder,'obs*'],'Choose obs file:'); % obs trajectory
            % Check compatibility of truth and obs
            compatible = new_check_obs_compatibility(nature_parms.ftruth, nature_parms.fobs);
            if ~compatible
                fprintf('Obs file\n%s\nincompatible with truth file\n%s\n, retrying...\n',...
                    nature_parms.fobs, nature_parms.ftruth);
            end
        end
    end
    % Extract Lorenz 96 Model I, II or III parameters from nature run filename
    string = strsplit(nature_parms.ftruth,'_');
    nature_parms.K = str2double(string{3}(2:end));
    nature_parms.F = str2double(string{4}(2:end));
    nature_parms.I = str2double(string{5}(2:end));
    nature_parms.b = str2double(string{6}(2:end));
    nature_parms.c = str2double(string{7}(2:end));
    nature_parms.timestep = str2double(string{8}(3:end));
    nature_parms.spinup = str2double(string{9}(7:end));
    nature_parms.seed = str2double(string{11}(5:end));
end

%% Model parameters
function model_parms = get_model_parms(nature_parms)
    % Allow experiments with parameters that differ from
    % the nature run parameters
    name='Model Parameters';
    numlines=1;
    prompt={'K(parm)','I(parm)','F(orcing)','b (damping)',...
            'c (coupling)','timestep'};
    default={num2str(nature_parms.K),num2str(nature_parms.I),...
        num2str(nature_parms.F),num2str(nature_parms.b),...
        num2str(nature_parms.c),num2str(nature_parms.timestep)};
    answer=inputdlg(prompt,name,numlines,default);i=1;
    model_parms.K = str2double(answer{i});i=i+1;
    model_parms.I = str2double(answer{i});i=i+1;
    model_parms.F = str2double(answer{i});i=i+1;
    model_parms.b = str2double(answer{i});i=i+1;
    model_parms.c = str2double(answer{i});i=i+1;
    model_parms.timestep = str2double(answer{i});
end

%% Parallelization options and initialization
function parprocs = parallelize(reserve)
    NPROCS = feature('numcores');
    parprocs = max(1, NPROCS - abs(int(reserve))); % Use all processes available
    if parprocs > 1 && isempty(gcp)
        parpool(parprocs);
    end
    fprintf('Started: %s\n',datestr(clock));
end

%% Plot time series of truth and ensemble mean
    function plot_ensmean_truth(outfolder,Xt,K,alpha,tF,R,locstr,...
                            locrad,Ncycles,first,skip,parms,myseed)
    myfilt = 'GKSQ';
    Nx = size(Xt,1);
    posterior = [outfolder,'K',num2str(K,'%d\n'),'\',...
        'posterior_alpha_',num2str(alpha,'%5.3f\n'),...
        '_tf_',num2str(tF,'%4.2f\n'),...
        '_R_',num2str(R(2,2),'%4.2f\n'),...
        '_loc_',locstr,num2str(locrad,'%d\n'),...
        '_nc_',num2str(Ncycles,'%d\n'),...
        '_1st_',num2str(first,'%d\n'),...
        '_skip_',num2str(skip,'%d\n'),...
        '_Nx_',num2str(Nx,'%d\n'),...
        '_Kp_',num2str(parms.K,'%d\n'),...
        '_seed_',num2str(myseed,'%d\n'),...
        '_',myfilt,'.mat'];
    Xa = load(posterior);
    fn = fieldnames(Xa);
    Xa = Xa.(fn{1});
    figure;
    plot(Xt(:,1),'k-');
    hold on; grid on;
    plot(mean(Xa(:,:,1),2),'r-');
    ylim([-15,20]);
    xlim([1,Nx]);
    xtt = floor(Nx/16);
    xticks([1 xtt:xtt:Nx]);
    legend('Truth','Ensmean');
    title(sprintf('%s time=1',myfilt));
    tskip = 4; % Plot once/day
    for t=tskip+1:tskip:Ncycles
        pause(0.5);
        clf;
        plot(Xt(:,t),'k-');
        hold on; grid on;
        plot(mean(Xa(:,:,t),2),'r-');
        ylim([-15,20]);
        xlim([1,Nx]);
        xtt = floor(Nx/16);
        xticks([1 xtt:xtt:Nx]);
        legend('Truth','Ensmean');
        title(sprintf('%s time=%d',myfilt,t));
    end % for t=2:Ncycles
end

%% Observation biases
function y = apply_obs_bias(obs_parms,y)
   % Multiplicative and additive obs biases
   y = y.*obs_parms.obbiasfac + obs_parms.obbias;
end
