function [ensmean,post_stmse,post_stvarmse,varargout] =...
    parDA_GKSQ(XIC,Xt,y,outfolder,H,Rhat,biascor,da,model,nature,run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bill Campbell Last Modified 9/21/2022
% Global-solve Kalman Square Root filter
% Added covariance localization
% Added ensemble variance for comparison with squared error
% Choose initial ensemble from sample of truth run
% Parallelize over ensemble members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Our standard choice: very good precision is enforced
options=odeset('RelTol', model.reltol, 'AbsTol', model.abstol);

% define temporal grid for the solver
% Can do this outside the ncycle loop, as the lorenz models have no
% explicit time dependence
Nt = 3;
tsteps = linspace(0, model.timestep, Nt);

[Nx,Ncycles]=size(Xt);

% Define covariance localization matrix CL
periodic = true;
CL = covloc(Nx, da.loctype, da.locrad, periodic);

% Allocate memory
ZZ = zeros(Nx,da.K);
% Initial conditions
XX = XIC(:,1:da.K);
ZKSQ = zeros(Nx,da.K,Ncycles);
ZKSQ(:,:,1) = XX;
XKSQ = zeros(Nx,da.K,Ncycles);
XKSQ(:,:,1) = XX;
scoreKSQ = zeros(Nx*2,Ncycles);
ensvarKSQ = scoreKSQ;
errKSQ = scoreKSQ;

fprintf('Cycle = 1 of %d: da.K=%d,model.timestep=%4.2f,da.alpha=%5.3f,seed=%d\n',...
    Ncycles,da.K,model.timestep,da.alpha,nature.seed);
tstart=tic;
first = true;
loranon = @(t, x) circ_lorenz2005(t, x, model);
for ncycle = 1:Ncycles-1
    %%%%%%%%% Global Kalman Square Root
    %%%%%%%% build ensemble members
    parfor kk = 1:da.K
        % integrate the equations with one of the available integrators, in this
        % case the Runga-Kutta 4,5 method (good for simple, non-stiff systems).
        M0 = XX(:, kk); % posterior ensemble member from previous time step
        [~,M] = ode45(loranon, tsteps, M0, options);
        ZZ(:, kk) = M(end, :).'; % forecast (background) ensemble member, which is input as the prior to the DA routine
    end
    ZKSQ(:, :, ncycle+1) = ZZ; % forecast (background) ensemble, which is input as the prior to the DA routine

    % Do DA every da.cycle_skip time steps
    if ~mod(ncycle,da.cycle_skip)
        [XX,scoreKSQ(1:2*Nx,ncycle+1),...
            ensvarKSQ(1:2*Nx,ncycle+1),...
            errKSQ(1:2*Nx,ncycle+1)] = ...
            KSQ(ZZ,da.alpha,da.K,H,Rhat,y(:,ncycle+1),Xt(:,ncycle+1),...
                CL,biascor);
    else
        XX = ZZ; % No DA, set posterior equal to prior
    end
    XKSQ(:, :, ncycle+1) = XX; % posterior

    % Diagnostic output
    if ~mod(ncycle,run.printcycle) || ncycle==Ncycles-1
        % Compute mean value so far
        % Will want to add calculations for ensvar also?
        [time_avg,time_stdev] =...
            mynanstats(scoreKSQ(:,ncycle-run.printcycle+2:ncycle+1).'); % 1 x 2*Nx
        [time_vavg,~] =...
            mynanstats(ensvarKSQ(:,ncycle-run.printcycle+2:ncycle+1).');
        if isscalar(time_stdev)
            time_stdev = repmat(time_stdev,1,2*Nx);
        end
        % Normalize by climvar for plotting (climvar = 1 for now)
        climvar = 1;
        prior_space_mse_norm   = mean(time_avg(1:Nx))/climvar; % scalar
        post_space_mse_norm   = mean(time_avg(Nx+1:2*Nx))/climvar; % scalar
		time_var = time_stdev.^2;
        prior_space_varse_norm = mean(time_var(1:Nx))/climvar; % scalar
        post_space_varse_norm = mean(time_var(Nx+1:2*Nx))/climvar; % scalar
        % Estimate timing
        telapsed=toc(tstart);
        fprintf('Cycle = %d of %d: da.K=%d,model.timestep=%4.2f,da.alpha=%5.3f,seed=%d,time=%5.2f seconds\n',...
            ncycle,Ncycles,da.K,model.timestep,da.alpha,nature.seed,telapsed);
        display_results('GKSQ',time_avg.',time_vavg.');
        % First fix plot_results, then add ensvar
        plot_results(first,ncycle,prior_space_mse_norm,prior_space_varse_norm,...
                     post_space_mse_norm,post_space_varse_norm,...
                     Ncycles,run.printcycle,da.ci,da.K);
        if first
            t_estimate = Ncycles./run.printcycle.*telapsed;
            thour=floor(t_estimate/3600);
            tminute=floor((t_estimate-thour*3600)/60);
            tsecond=round(mod(t_estimate,60));
            fprintf('Estimated time to completion is %d hours %d minutes %d seconds\n',...
                thour,tminute,tsecond);
            first = false;
        end
    end % diagnostic output
end % end number of cycles

fprintf('Mean Posterior Error after %d-cycle spinup\n',da.spinup);
[tmse,tbiassq,stmse,stbiassq] =...
    compute_stats(errKSQ, da.spinup);  % tmse,tbiassq are 1 x 2*Nx
post_stmse = stmse(2); % ring mean of time mean of posterior squared error
if nargout >= 4
    varargout{1} = 100*stbiassq(2)./stmse(2); % Percent due to squared bias
end
[tvarmse,tvarbiassq,stvarmse,stvarbiassq] =...
    compute_stats(ensvarKSQ, da.spinup);  % tvarmse,tvarbiassq are 1 x 2*Nx
post_stvarmse = stvarmse(2);  % ring mean of time mean of posterior ensemble variance
if nargout >= 5
    varargout{2} = 100*stvarbiassq(2)./stvarmse(2); % Percent due to squared bias
end
display_results('GKSQ',tmse(Nx+1:2*Nx),tvarmse(Nx+1:2*Nx));

% Time mean ensemble means of differences
[AmB,OmB,OmA,AmT,BmT,OmT] = diff_stats(XKSQ,ZKSQ,y,Xt,da.spinup);
% Print spatial means of time means
fprintf('Mean A-B is %d, mean O-B is %d, mean O-A is %d after %d-cycle spinup\n',...
    mean(AmB), mean(OmB), mean(OmA), da.spinup);
fprintf('Mean A-T is %d, mean B-T is %d, mean O-T is %d after %d-cycle spinup\n',...
    mean(AmT), mean(BmT), mean(OmT), da.spinup);

% Always save the following stats
fprintf('Saving output stats...\n');
generic = fullfile(outfolder,...
        ['K',num2str(da.K,'%d')],...
        [run.expname,'_generic_GKSQ.mat']);
errstats = strrep(generic,'generic','errstats');
save(errstats,'errKSQ','tmse','tbiassq','stmse','stbiassq');
ensvarstats = strrep(generic,'generic','ensvarstats');
save(ensvarstats,'ensvarKSQ','tvarmse','tvarbiassq','stvarmse','stvarbiassq');
statsfile = strrep(generic,'generic','diffstats');
save(statsfile,'AmB','OmB','OmA','AmT','BmT','OmT','-v7.3');
% Optional save obs
if ~run.use_obs_file
    fprintf('Saving obs (includes spinup)...\n');
    obsfile = strrep(generic,'generic','obs');
    save(obsfile,'y');
end
% Optional save full prior, posterior
if (run.save_state)    % Create filenames, and open files
    fprintf('Saving prior and posterior full states (includes spinup)...\n');
    prior = strrep(generic,'generic','prior');
    save(prior,'ZKSQ','-v7.3');
    posterior = strrep(generic,'generic','posterior');
    save(posterior,'XKSQ','-v7.3');
end

% Compute ensemble mean of posterior for plotting
ensmean = squeeze(mean(XKSQ,2));

end

function [tmse,tbiassq,stmse,stbiassq] =...
        compute_stats(space_time, spinup)
    Nx = size(space_time,1) / 2; % prior and posterior
    [bias,~,rmse] = mynanstats(space_time(:,spinup+1:end).');  % 1 x 2*Nx
    tbiassq = bias.^2;  % 1 x 2*Nx
    prior_stbiassq = mean(tbiassq(1:Nx));
    post_stbiassq = mean(tbiassq(Nx+1:2*Nx));
    stbiassq = [prior_stbiassq post_stbiassq];  % 1 x 2
    tmse = rmse.^2;  % 1 x 2*Nx
    prior_stmse = mean(tmse(1:Nx));
    post_stmse = mean(tmse(Nx+1:2*Nx));
    stmse = [prior_stmse post_stmse];  % 1 x 2
end

function display_results(tstr,aerr,varargin)
    Nx=length(aerr)/2;
    nmid = floor(Nx/2);
    fprintf('%s: Prior squared error %s = %s + %s + ... + %s + %s + ...\n',...
        tstr,num2str(sum(aerr(1:Nx))),num2str(aerr(1)),num2str(aerr(2)),num2str(aerr(nmid)),num2str(aerr(nmid+1)));
    if nargin==3
        avar = varargin{1};
        fprintf('%s: Prior error variance %s = %s + %s + ... + %s + %s + ...\n',...
            tstr,num2str(sum(avar(1:Nx))),num2str(avar(1)),num2str(avar(2)),num2str(avar(nmid)),num2str(avar(nmid+1)));
    end
    fprintf('%s: Posterior squared error %s = %s + %s + ... + %s + %s + ...\n',...
        tstr,num2str(sum(aerr(Nx+1:2*Nx))),num2str(aerr(Nx+1)),num2str(aerr(Nx+2)),num2str(aerr(Nx+nmid)),num2str(aerr(Nx+nmid+1)));
    if nargin==3
        fprintf('%s: Posterior error variance %s = %s + %s + ... + %s + %s + ...\n',...
            tstr,num2str(sum(avar(Nx+1:2*Nx))),num2str(avar(Nx+1)),num2str(avar(Nx+2)),num2str(avar(Nx+nmid)),num2str(avar(Nx+nmid+1)));
    end
end

