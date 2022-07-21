%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parLorGKSQ with prior inflation
% Using Lorenz 96 Model I, II or III
% Global Kalman Square Root filter
% Bill Campbell
% Last modified 7/21/2022

%% Parameter input
[Ncycles, first, skip, oberr, Kvec, ci, alphavec, tFvec, fstr, seedlist,...
 printcycle, loctype, locstr, locradlist, savestate, plotit] =...
input_params;

%% Nature run and observations of the truth
[outfolder, ftruth, truepath, fobs, obspath, nature, spinup] =...
get_nature_run_info;

%% Parameter adjustment
parms = adjust_parameters(nature);

%% Output directory creation
make_ensemble_dirs(Kvec, outfolder)

%% Parallelization options and initialization
[fname, fid, tot_err, tot_var, start, parprocs] = parallelize(alphavec);

%% Main loop
% Forward model time steps, NOT a time loop through Ncycles
for itf = 1:length(tFvec)
    tF = tFvec(itf);
    % Load nature run and observations for this time step
    ftruth = regexprep(ftruth,'tf_....',['tf_',num2str(tF,'%4.2f\n')]);
    fobs = regexprep(fobs,'tf_....',['tf_',num2str(tF,'%4.2f\n')]);
    [Xt,abstol,reltol] = load_truth(truepath,ftruth); % full nature run
    Xt = Xt(:,spinup+1:end); % discard nature run spinup
    parms.abstol = abstol;
    parms.reltol = reltol;
    y = load_obs(obspath,fobs);  % observations consistent with nature run
    y = y(:,spinup+1:end); % must discard same spinup perio
    % Compute climatological covariance
    % Take time mean at each N point
    tru_t_ave = mean(Xt,2); % Nx x 1
    mean_diff = (bsxfun(@minus,y,tru_t_ave)).^2;  % Nx x Ncycles
    climvar = mean(mean_diff,'all'); % scalar
    fprintf('Climatoligical error variance is %8.6f\n',climvar);
    % Generate initial ensemble from nature run climatology
    % Equally spaced through time
    ics = floor(linspace(1,size(Xt,2),max(Kvec)+1)); % Kmax+1 x 1
    XIC = Xt(:,ics(2:end)); % Nx x Kmax
    if (Ncycles > size(y,2))
        fprintf('Error -- obs file not long enough: Ncycles=%d\n > len(obs)',...
            Ncycles,size(y,2));
        break;
    else % Only use needed portion of nature run and obs
        Xt = Xt(:,1:Ncycles);  % Nx x Ncycles
        y = y(:,1:Ncycles);  % Nx x Ncycles
    end
    % Set oblist and rlist here
    Nx = size(y,1);
    oblist = zeros(1,Nx);
    rlist = oblist;
    % Observed gridpoints
    gmask = first:skip:Nx;
    oblist(gmask)=1;
    rlist(gmask)=oberr;
    % Set up observation operator and ob error covariance matrix
    % H will be Nobs x Nx, R will be Nobs x Nobs, Nobs <= Nx
    [~,H,R] = forward(oblist,rlist);
    % Observed time steps -- implement later
    %tmask = first:tskip:Ncycles; this does not work
    % Compress the observation vector to match compressed R
    y = y(gmask,:);  % Nobs x Ncycles
    % Ensemble sizes
    for iK = 1:length(Kvec)
        K = Kvec(iK);
        % Random seeds
        for iseed = 1:length(seedlist)
            myseed = seedlist(iseed);
            % Localization radii
            for irad = 1:length(locradlist)
                locrad = locradlist(irad);
                for ialpha = 1:length(alphavec) % Must set filenames outside parfor loop if myprocs < length(alphavec)
                    alpha = alphavec(ialpha);
                    fname{ialpha,:} = [outfolder,'K',num2str(K,'%d\n'),...
                        '\alpha_',num2str(alpha,'%4.2f\n'),...
                        '_tF_',num2str(tF,'%4.2f\n'),...
                        '_loc_',locstr,num2str(locrad,'%d\n'),...
                        '_nc_',num2str(Ncycles,'%d\n'),...
                        '_1st_',num2str(first,'%d\n'),...
                        '_skip_',num2str(skip,'%d\n'),...
                        '_',num2str(myseed,'%d\n'),'.txt'];
                end
                if parprocs > 1
                    % Multiplicative inflation values
                    parfor (ialpha = 1:length(alphavec),parprocs) % Prior inflation factors
                        palpha = alphavec(ialpha);
                        fid(ialpha) = fopen(fname{ialpha,:},'at');
                        fprintf(fid(ialpha),'K=%d, tF=%4.2f, loc %s, rad %d, %d cycles\t%s\t%s\t%s\t%s\t%s\n',...
                            K,tF,locstr,locrad,Ncycles,fstr); %#ok<*PFBNS>
                        tstart=tic;
                        % Execute up to Ncycles-1 DA cycles
                        [tot_err(ialpha,:),tot_avar(ialpha,:)] = ...
                            parDA_GKSQ(XIC,Xt,y,outfolder,H,R,...
                            first,skip,tF,ci,palpha,K,myseed,fstr,savestate,...
                            printcycle,loctype,locstr,locrad,parms,nature,...
                            climvar);
                        terr = tot_err(ialpha,:);
                        elapsed(tstart,palpha,ialpha,terr)
                    end % alpha
                else
                    % Multiplicative inflation values
                    for ialpha = 1:length(alphavec) % Prior inflation factors
                        alpha = alphavec(ialpha);
                        fid(ialpha) = fopen(fname{ialpha,:},'at');
                        fprintf(fid(ialpha),'K=%d, tF=%4.2f, loc %s, rad %d, %d cycles\t%s\t%s\t%s\t%s\t%s\n',...
                            K,tF,locstr,locrad,Ncycles,fstr);
                        tstart=tic;
                        % Execute up to Ncycles-1 DA cycles
                        [tot_err(ialpha,:),tot_avar(ialpha,:)] = ...
                            parDA_GKSQ(XIC,Xt,y,outfolder,H,R,...
                            first,skip,tF,ci,alpha,K,myseed,fstr,savestate,...
                            printcycle,loctype,locstr,locrad,parms,nature,...
                            climvar);
                        terr = tot_err(ialpha,:);
                        elapsed(tstart,alpha,ialpha,terr)
                    end % alpha
                    % Plot time series of truth and ensemble mean
                    if plotit
                        plot_ensmean_truth(outfolder,Xt,K,alpha,tF,R,...
                            locrad,Ncycles,first,skip,myseed)                  
                    end % if plotit
                end % myprocs > 0
            end % localization radii
        end % seeds
    end % ensemble sizes
end % different time steps, NOT a time loop through Ncycles
% if ~isempty(gcp('nocreate')),
%     delete(gcp);
% end
finish = datestr(clock);
fprintf('Started: %s, Finished %s\n',start,finish);
disp([prompt',answer]);
fclose('all');

%% Parameter input
function [Ncycles, first, skip, oberr, Kvec, ci, alphavec,...
          tFvec, fstr, seedlist, printcycle, loctype, locstr,...
          locradlist, savestate, plotit] =...
          input_params
    prompt={'Cycles','First var observed','Grid Skip','Time Skip',...
        'Ob error','Ensemble size(s)','Confidence level(s)','Prior inflation(s)',...
        'Cycling interval(s0','Seed(s)','Cycles per print','Localization Type',...
        'Localization radi(us,i)','Save state','Plot'};
    name='Filter Input';
    numlines=1;
    default={'1000','1','1','1',...
        '1.0','[500]','[0.95]','[1.64]',...
        '0.05','29418','100','none',...
        '[10]','1','0'};
    answer=inputdlg(prompt,name,numlines,default);i=1;
    % Number of cycles to run
    Ncycles = str2num(answer{i});i=i+1;
    % First observed variable
    first = str2num(answer{i});i=i+1;
    % Variables skipped (unobserved)
    skip = str2num(answer{i});i=i+1;
    % Time steps skipped (completely unobserved)
    tskip = str2num(answer{i});i=i+1;
    % observation error variance
    oberr = str2num(answer{i});i=i+1;
    % number of ensemble members --- 32, 64, 128, 256, 512, 1024
    Kvec = str2num(answer{i});i=i+1;
    % Confidence level (0.95 is typical)
    ci = str2num(answer{i});i=i+1;
    % inflation parameters (tune these for each experiment)
    alphavec = str2num(answer{i});i=i+1;
    % cycling interval -- should get this from nature run
    tFvec = str2num(answer{i});i=i+1;
    fstr = 'Global Kalman Square Root';
    seedlist = str2num(answer{i});i=i+1;
    printcycle = str2num(answer{i});i=i+1;
    % Localization type
    switch lower(answer{i})
        case {'gaspari','gasparicohn','g','gc','gaspari-cohn','g-c'}
            loctype = 1;
            locstr = 'g';
        case {'boxcar','box','flat','b'}
            loctype = 2;
            locstr = 'b';
        case {'triangular','tri','t','tent'}
            loctype = 3;
            locstr = 't';
        otherwise % No localization
            loctype = 0;
            locstr = 'n';
    end
    i=i+1;
    % Localization radius (in gridpoints)
    locradlist = abs(round(str2num(answer{i})));i=i+1;
    savestate = logical(str2num(answer{i}));i=i+1;
    plotit = logical(str2num(answer{i}));i=i+1;
    disp([prompt',answer]);
end

%% Nature run and observations of the truth
function [outfolder, ftruth, truepath, fobs, obspath, nature, spinup] =...
          get_nature_run_info
% Truth file e.g. L05M3_N960_K32_F15.00_I12_b10.00_c2.50_tf0.05_spinup100_tsteps10000_seed51422
    % Obs file e.g. obs_tf0.05_nc10000_R1.00_N960_K32_F15.00_I12_b10.00_c2.50_tsteps960_tseed51422_oseed73033
    % Nature run location
    outfolder = 'C:\Users\campbell\Documents\MATLAB\Lorenz_96_model\'; % Local hard drive
    % Load truth and compatible observations
    compatible = false;
    while ~compatible
        [ftruth,truepath]=uigetfile([outfolder,'*L05*'],'Choose truth file:'); % truth trajectory
        [fobs,obspath]=uigetfile([outfolder,'obs*'],'Choose obs file:'); % obs trajectory
        % Check compatibility of truth and obs
        compatible = new_check_obs_compatibility(ftruth,fobs);
    end
    % Extract Lorenz 96 Model I, II or III parameters from nature run filename
    string = strsplit(ftruth,'_');
    nature.K = str2double(string{3}(2:end));
    nature.F = str2double(string{4}(2:end));
    nature.I = str2double(string{5}(2:end));
    nature.b = str2double(string{6}(2:end));
    nature.c = str2double(string{7}(2:end));
    % Also need nature run spinup
    spinup = str2double(string{9}(7:end));
end

%% Parameter adjustment
function parms = adjust_parameters(nature)
    % Allow experiments with parameters that differ from
    % the nature run parameters
    prompt={'K(parm)','I(parm)','F(orcing)','b (damping)','c (coupling)'};
    name='Parameter adjustment';
    numlines=1;
    default={num2str(nature.K),num2str(nature.I),num2str(nature.F),...
        num2str(nature.b),num2str(nature.c)};
    answer=inputdlg(prompt,name,numlines,default);i=1;
    parms.K = str2double(answer{i});i=i+1;
    parms.I = str2double(answer{i});i=i+1;
    parms.F = str2double(answer{i});i=i+1;
    parms.b = str2double(answer{i});i=i+1;
    parms.c = str2double(answer{i});
end

%% Parallelization options and initialization
function [fname, fid, tot_err, tot_var, start, parprocs] =...
          parallelize(alphavec)
    NPROC = feature('numcores'); % This machine has 4 processors
    reserve = 1;
    allprocs = NPROC - reserve; % Use all processes available
    % If have enough processors, parallelize over both the alpha parameter and
    % number of ensemble members
    parprocs = min([allprocs,length(alphavec)]); % Devote max resources to parameters, devote remainder to ensemble members
    % To run with only one processor, uncomment the two lines below
    % parprocs = 0;
    % myprocs = 1;
    myprocs = parprocs; % Can be less then allprocs for small ensembles
    if myprocs > 1 && isempty(gcp)
        parpool(myprocs);
    end
    fname = cell(length(alphavec),1);
    fid = zeros(length(alphavec),1);
    tot_err = zeros(length(alphavec),1);
    tot_var = tot_err;
    start = datestr(clock);
    fprintf('Started: %s\n',start);
end

%% Output directory creation
function make_ensemble_dirs(Kvec, outfolder)
    for iK = 1:length(Kvec) % Ensemble sizes
        K = Kvec(iK);
        if ~exist([outfolder,'K',num2str(K,'%d\n')], 'dir')
            mkdir([outfolder,'K',num2str(K,'%d\n')]);
        end
    end
end

%% Plot time series of truth and ensemble mean
function plot_ensmean_truth(outfolder,Xt,K,alpha,tF,R,...
                            locrad,Ncycles,first,skip,myseed)
    myfilt = 'GKSQ';
    posterior = [outfolder,'K',num2str(K,'%d\n'),'\',...
        'posterior_alpha_',num2str(alpha,'%4.2f\n'),...
        '_tf_',num2str(tF,'%4.2f\n'),...
        '_R_',num2str(R(2,2),'%4.2f\n'),...
        '_loc_',locstr,num2str(locrad,'%d\n'),...
        '_nc_',num2str(Ncycles,'%d\n'),...
        '_1st_',num2str(first,'%d\n'),...
        '_skip_',num2str(skip,'%d\n'),...
        '_seed_',num2str(myseed,'%d\n'),...
        '_',myfilt,'.mat'];
    Xa = load(posterior);
    fn = fieldnames(Xa);
    Xa = Xa.(fn{1});
    figure;
    plot(Xt(:,1),'k-');
    hold on; grid on;
    plot(mean(Xa(:,:,1),2),'r-');
    legend('Truth','Ensmean');
    title(sprintf('%s time=1',myfilt));
    for t=2:Ncycles
        pause(0.5);
        clf;
        plot(Xt(:,t),'k-');
        hold on; grid on;
        plot(mean(Xa(:,:,t),2),'r-');
        legend('Truth','Ensmean');
        title(sprintf('%s time=%d',myfilt,t));
    end % for t=2:Ncycles
    pause;
end

%% Timing and cleanup
function elapsed(tstart,alpha,ialpha,terr)
    tend=toc(tstart);
    thour=floor(tend/3600);
    tminute=floor((tend-thour*3600)/60);
    tsecond=round(mod(tend,60));
    fprintf('Elapsed time %0.2d:%0.2d:%0.2d\n',thour,tminute,tsecond);
    fprintf('%4.2f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
        alpha,terr);
    fprintf(fid(ialpha),'%4.2f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
        alpha,terr);
    fclose(fid(ialpha));
end