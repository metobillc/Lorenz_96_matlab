function [tot_err,tot_avar,varargout] = parDA_GKSQ(XIC,Xt,y,outfolder,...
          H,R,g1st,skip,tF,ci,alpha,K,myseed,fstr,savestate,printcycle,...
          loctype,locstr,locrad,parms,trueparms,obs_parms,climvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bill Campbell Last Modified 9/4/2022
% Global-solve Kalman Square Root filter
% Added covariance localization
% Added ensemble variance for comparison with squared error
% Choose initial ensemble from sample of truth run
% Parallelize over ensemble members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Skipping DA cycles -- nskip = 1 turns it off
nskip = 1;

% Our standard choice: very good precision is enforced
% Try to use the same tolerances for DA as in the nature run
if isfield(parms, 'reltol') && isfield(parms, 'abstol')
  options=odeset('RelTol',parms.reltol,'AbsTol',parms.abstol);
else
  options=' ';
end

% define temporal grid for the solver
% Can do this outside the ncycle loop, as the lorenz models have no
% explicit time dependence
Nt = 3;
tsteps = linspace(0,tF,Nt);

% initialize random number generators
rng(myseed);

[Nx,Ncycles]=size(Xt);

% Define covariance localization matrix CL
periodic = true;
CL = covloc(Nx,loctype,locrad,periodic);


% Allocate memory
ZZ = zeros(Nx,K);
% Initial conditions
XX = XIC(:,1:K);
ZKSQ = zeros(Nx,K,Ncycles);
ZKSQ(:,:,1) = XX;
XKSQ = zeros(Nx,K,Ncycles);
XKSQ(:,:,1) = XX;
scoreKSQ = zeros(Nx*2,Ncycles);
ensvarKSQ = scoreKSQ;
errKSQ = scoreKSQ;

fprintf('Cycle = 1 of %d: K=%d,tF=%4.2f,alpha=%5.3f,seed=%d\n',...
    Ncycles,K,tF,alpha,myseed);
first = true;
ci = abs(ci);
if (ci > 1.0) % test for pctg rather than decimal
    ci = 100*ci;
end
if (ci > 1.0)
    ci = 0.95; % default to 95%
end
save_printcycle = printcycle;
printcycle = min(printcycle,1000);
detect_lorenz2005(parms);

tstart=tic;
loranon = @(t, x) circ_lorenz2005(t, x, parms);
for ncycle = 1:Ncycles-1
    %%%%%%%%% Global Kalman Square Root
    %%%%%%%% build ensemble members
    parfor k=1:K
        % integrate the equations with one of the available integrators, in this
        % case the Runga-Kutta 4,5 method (good for simple, non-stiff systems).
        M0 = XX(:,k); % posterior ensemble member from previous time step
        [~,M] = ode45(loranon,tsteps,M0,options);
        ZZ(:,k) = M(end,:).'; % forecast (background) ensemble member, which is input as the prior to the DA routine
    end
    ZKSQ(:,:,ncycle+1) = ZZ; % forecast (background) ensemble, which is input as the prior to the DA routine

    % function [Xa,score,ensvar] = KSQ(Zfull,alpha,K,H,R,y,Xt,varargin)

    if ~mod(ncycle,nskip) % Do DA every nskip time steps
        [XX,scoreKSQ(1:2*Nx,ncycle+1),...
            ensvarKSQ(1:2*Nx,ncycle+1),...
            errKSQ(1:2*Nx,ncycle+1)] = ...
            KSQ(ZZ,alpha,K,H,R,y(:,ncycle+1),Xt(:,ncycle+1),CL);
    else
        XX = ZZ; % No DA, set posterior equal to prior
    end
    XKSQ(:,:,ncycle+1) = XX; % posterior
    
    % Diagnostic output
    if ~mod(ncycle,printcycle) || ncycle==Ncycles-1
        % Compute mean value so far
        % Will want to add calculations for ensvar also?
        [time_avg,time_stdev] =...
            mynanstats(scoreKSQ(:,ncycle-printcycle+2:ncycle+1).'); % 1 x 2*Nx
        [time_vavg,~] =...
            mynanstats(ensvarKSQ(:,ncycle-printcycle+2:ncycle+1).');
        if isscalar(time_stdev)
            time_stdev = repmat(time_stdev,1,2*Nx);
        end
        % Normalize by climvar for plotting
        prior_space_mse_norm   = mean(time_avg(1:Nx))/climvar; % scalar
        post_space_mse_norm   = mean(time_avg(Nx+1:2*Nx))/climvar; % scalar
		time_var = time_stdev.^2;
        prior_space_varse_norm = mean(time_var(1:Nx))/climvar; % scalar
        post_space_varse_norm = mean(time_var(Nx+1:2*Nx))/climvar; % scalar
        % Estimate timing
        telapsed=toc(tstart);
        fprintf('Cycle = %d of %d: K=%d,tF=%4.2f,alpha=%5.3f,seed=%d,time=%5.2f seconds\n',...
            ncycle,Ncycles,K,tF,alpha,myseed,telapsed);
        display_results(fstr,time_avg.',time_vavg.');
        % First fix plot_results, then add ensvar
        plot_results(first,ncycle,prior_space_mse_norm,prior_space_varse_norm,...
                     post_space_mse_norm,post_space_varse_norm,...
                     Ncycles,printcycle,ci,K);
        if first
            t_estimate = Ncycles./printcycle.*telapsed;
            thour=floor(t_estimate/3600);
            tminute=floor((t_estimate-thour*3600)/60);
            tsecond=round(mod(t_estimate,60));
            fprintf('Estimated time to completion is %d hours %d minutes %d seconds\n',...
                thour,tminute,tsecond);
            printcycle = save_printcycle;
            first = false;
        end
    end % diagnostic output
    
end % end number of cycles

spinup = 100;
fprintf('Mean Posterior Error after %d-cycle spinup\n',spinup);
[aerr,stdev] = mynanstats(scoreKSQ(:,spinup+1:end).');
tot_err = sum(aerr(Nx+1:2*Nx));
tot_var = sum(stdev(Nx+1:2*Nx).^2);
[avar,vstdev] = mynanstats(ensvarKSQ(:,spinup+1:end).');
tot_avar = sum(avar(Nx+1:2*Nx));
tot_vvar = sum(vstdev(Nx+1:2*Nx).^2);
display_results(fstr,aerr,avar);

% Save states in .mat files
if (savestate)    % Create filenames, and open files
    tsave=tic;
    fprintf('Saving prior, posterior, parameters...');
    prior = [outfolder,'K',num2str(K,'%d\n'),'\',...
        'prior_alpha_',num2str(alpha,'%5.3f\n'),...
        '_tf_',num2str(tF,'%4.2f\n'),...
        '_R_',num2str(R(2,2),'%4.2f\n'),...
        '_loc_',locstr,num2str(locrad,'%d\n'),...
        '_nc_',num2str(ncycle+1,'%d\n'),...
        '_1st_',num2str(g1st,'%d\n'),...
        '_skip_',num2str(skip,'%d\n'),...
        '_Nx_',num2str(Nx,'%d\n'),...
        '_Kp_',num2str(parms.K,'%d\n'),...
        '_seed_',num2str(myseed,'%d\n')];
    newprior = [prior,'_GKSQ'];
    save([newprior,'.mat'],'ZKSQ','-v7.3');
    posterior = strrep(prior,'prior','posterior');
    newpost  = [posterior,'_GKSQ'];
    % Saving errKSQ instead of scoreKSQ
    save([newpost,'.mat'],'XKSQ','errKSQ','ensvarKSQ','-v7.3');
    parmfile = strrep(prior,'prior','trueparms');
    newparmfile = [parmfile,'_GKSQ'];
    save([newparmfile,'.mat'],'parms','trueparms','obs_parms','-v7.3');
    fprintf('took %5.3f seconds.\n',toc(tsave));
end % if (savestate)



if nargout>=3
    varargout{1} = tot_var;
end
if nargout==4
    varargout{2} = tot_vvar;
end
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
