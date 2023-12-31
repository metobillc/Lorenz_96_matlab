% This script loads posterior .mat files, then calls 
% plot_mse_ring.m to produce a ring plots of prior & posterior

%% Plot options
prompt={'MSE Logscale','MSE Plot Prior', 'Ensvar Logscale', 'Ensvar Plot Prior'};
name='Ring Plot Options';
numlines=1;
default={'0','1','0','1'};
answer=inputdlg(prompt,name,numlines,default);i=1;
% MSE Logscale
mse_logscale = str2double(answer{i});i=i+1;
% MSE Plot Prior
mse_plot_prior = str2double(answer{i});i=i+1;
% Ensvar Logscale
ensvar_logscale = str2double(answer{i});i=i+1;
% Ensvar Plot Prior
ensvar_plot_prior = str2double(answer{i});

%% Load posterior, discarding da_spinup
outfolder = 'C:\Users\campbell\Document\MATLAB\\Lorenz_96\'; % Local hard drive
[fpost,ppath]=uigetfile([outfolder,'K*/posterior*'],'Choose post file:'); % post
posterior = fullfile(ppath,fpost);
fprintf('Loading posterior...\n');
% Even though this file is labeled 'posterior', it contains prior
% information as well
load(posterior,'ensvarKSQ','scoreKSQ','errKSQ');
if exist('errKSQ','var')
    use_errKSQ = 1;
    scoreKSQ = errKSQ.^2;
else
    fprintf('Older code did not save errKSQ, so cannot plot it.\n');
    use_errKSQ = 0;
end
da_spinup = 100;
scoreKSQ = scoreKSQ(:,da_spinup+1:end); % 2*Nx x Nt
ensvarKSQ = ensvarKSQ(:,da_spinup+1:end);

if use_errKSQ==1
    errKSQ = errKSQ(:,da_spinup+1:end);
    %% MSE, Bias, and Variance of the state estimate
    % errKSQ holds differences between ensemble mean and truth,
    % for both prior (1:Nx) and posterior(Nx+1:2*Nx)
    [bias,stdev,rmse] = mynanstats(errKSQ.'); % outputs are 1 x 2*Nx
    tbiassq = bias.^2;
    tvar = stdev.^2;
    % We want to plot the time mean squared error at each spatial point, tmse
    tmse = rmse.^2; % equals tbiassq + tvar
else
    % scoreKSQ already IS the squared error, so we just need its mean
    % We cannot assess the components of squared bias and variance that
    % add up to this MSE
    [tmse,~,~] = mynanstats(scoreKSQ.'); % 1 x 2*Nx
end

%% Upper error bar for MSE
ci = 0.95;
upper = oneside_ci(scoreKSQ.', ci);

%% Plot MSE, confidence interval, and components of MSE
if use_errKSQ==1
    plot_mse_ring(tmse, upper, mse_logscale, mse_plot_prior, tbiassq);
else
    plot_mse_ring(tmse, upper, mse_logscale, mse_plot_prior);
end

%% Ensemble variance
[ens_bias,ens_stdev,ens_rmse] = mynanstats(ensvarKSQ.');
ens_tbiassq = ens_bias.^2;
ens_tvar = ens_stdev.^2; 
ens_tmse = ens_rmse.^2; % equals ens_tbiassq + ens_tvar

%% Confidence interval for ensemble variance
% Need to compute effective N, which is less than N
% because of temporal autocovariance
ens_upper = oneside_ci(ensvarKSQ.', ci);

%% Plot ensemble variance, confidence interval
plot_ensvar_ring(ens_tmse, ens_upper, ensvar_logscale, ensvar_plot_prior);

%% Nature run and observations of the truth
function [ftruth, truepath, fobs, obspath, spinup] =...
          get_nature_run_info(outfolder)
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
    %nature.K = str2double(string{3}(2:end));
    %nature.F = str2double(string{4}(2:end));
    %nature.I = str2double(string{5}(2:end));
    %nature.b = str2double(string{6}(2:end));
    %nature.c = str2double(string{7}(2:end));
    % Also need nature run spinup
    spinup = str2double(string{9}(7:end));
end

function upper = oneside_ci(tseries, ci)
    rho_all = lagged_autocovariance(tseries); % tseries is Nt x Nx
    rho = mean(rho_all); % spatial mean
    factor = (1 - rho) / (1 + rho);
    N = size(tseries,1);
    Neff = N * factor;
    fprintf('Rho = %d, Neff=%d\n',rho,Neff);
    [~,stdev,~] = mynanstats(tseries);
    stderr = stdev ./ sqrt(Neff);
    % Compute one-sided confidence interval
    cifac = norminv(ci);
    upper = cifac.*stderr;
end

function plot_mse_ring(mse, upper, logscale, plot_prior, varargin)
    if nargin==5
        sqbias = varargin{1};
        var = mse - sqbias;
    end
% MSE is 1 x 2*Nx
% First Nx values are the prior; second Nx values are the posterior
    Nx = size(mse,2)/2;
    figure('position',[100 200 1250 500]);
    prior = 1:Nx;
    post = Nx+1:2*Nx;
    ring = 1:Nx;
    if plot_prior==1
        plot(ring,mse(prior),'b-',ring,mse(prior)+upper(prior),'b:');
    end
    hold on
    plot(ring,mse(post),'r-',ring,mse(post)+upper(post),'r:')
    if nargin==5
        if plot_prior==1
            plot(ring,sqbias(prior),'c-',ring,var(prior),'y-');
        end
        plot(ring,sqbias(post),'k-',ring,var(post),'m-')
        gnir = [ring, fliplr(ring)];
        if plot_prior==1
            tween_prior = [var(prior), fliplr(sqbias(prior))];
            fill(gnir, tween_prior, 'y', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.0);
        end
        tween_post = [var(post), fliplr(sqbias(post))];
        fill(gnir, tween_post, 'm', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.0);
    end
    set(gca,'fontsize',12,'fontweight','bold')
    xlim([1,Nx]);
    if logscale==1
        miny = 0.01;
        set(gca,'yscale','log');
        ylim([miny,Inf]);
    else
        ylim([0,Inf]);
    end
    nticks = 15;
    gap = floor(Nx/nticks);
    xticks([1 gap:gap:Nx]);
    grid on
    xlabel('Ring')
    ylabel('Time Mean Squared Error')
    if nargin==5
        if plot_prior==1
            legend({'Prior','','Post','','Prior SqBias','Prior Var',...
                'Post SqBias','Post Var','',''},'location','northwestoutside');
        else
            legend({'Post','','Post SqBias','Post Var',''},...
                'location','northwestoutside');
        end
    else
        if plot_prior==1
            legend({'Prior','','Post',''},'location','northwest');
        else
            legend({'Post',''},'location','northwest');
        end
    end
    smse = mean(mse(post));
    ssqb = mean(sqbias(post));
    bpct = round(100 * ssqb / smse, 1);
    vpct = 100 - bpct;
    ttl1 = sprintf('Posterior Spatial Mean of TMSE: %8.5f',smse);
    if nargin==5
        ttl2 = sprintf('Squared Bias is %4.1f%%, Variance is %4.1f%%',...
            bpct,vpct);
        title({ttl1,ttl2});
    else
        title({ttl1});
    end
end

function plot_ensvar_ring(ensvar, upper, logscale, plot_prior)
% ensvar is 1 x 2*Nx
% First Nx values are the prior; second Nx values are the posterior
    Nx = size(ensvar,2)/2;
    figure('position',[300 600 1250 500]);
    prior = 1:Nx;
    post = Nx+1:2*Nx;
    ring = 1:Nx;
    if plot_prior==1
        plot(ring,ensvar(prior),'b-',ring,ensvar(prior)+upper(prior),'b:');
    end
    hold on
    plot(ring,ensvar(post),'r-',ring,ensvar(post)+upper(post),'r:')
    set(gca,'fontsize',12,'fontweight','bold')
    xlim([1,Nx]);
    if logscale==1
        miny = 0.01;
        set(gca,'yscale','log')
        ylim([miny,Inf]);
    else
        ylim([0,Inf]);
    end
    set(gca,'yscale','log');
    nticks = 15;
    gap = floor(Nx/nticks);
    xticks([1 gap:gap:Nx]);
    grid on
    xlabel('Ring')
    ylabel('Time Mean Ensemble Variance')
    legend({'Prior','','Post',''},'location','northwest')
    ttl1 = sprintf('Posterior Mean Ensemble Variance: %8.5f',...
        mean(ensvar(post)));
    title({ttl1})
end

