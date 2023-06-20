% This script loads posterior .mat files, then calls 
% plot_rmse_ring.m to produce a ring plots of prior & posterior
% As of 6/20/23, this script only calculates and prints RMSE and Bias

%% Plot options
%[rmse_logscale,rmse_plot_prior,ensstdev_logscale,ensstdev_plot_prior] = plot_options;

%% Load posterior, discarding da_spinup
outfolder = 'C:\Users\knisely\Desktop\Lorenz_05\'; % Local hard drive
[fpost,ppath]=uigetfile([outfolder,'K*/posterior*'],'Choose post file:'); % post
posterior = fullfile(ppath,fpost);
fprintf(fpost);
fprintf('\nLoading posterior...\n');
% Even though this file is labeled 'posterior', it contains prior
% information as well
load(posterior,'ensvarKSQ','errKSQ');
if exist('errKSQ','var')
    scoreKSQ = errKSQ.^2;
else
    fprintf('Older code did not save errKSQ. Quitting. \n');
    quit
end
da_spinup = 100;
scoreKSQ = scoreKSQ(:,da_spinup+1:end); % 2*Nx x Nt
ensvarKSQ = ensvarKSQ(:,da_spinup+1:end);
errKSQ = errKSQ(:,da_spinup+1:end);

%% RMSE, Bias, and Std Dev of the state estimate
% errKSQ holds differences between ensemble mean and truth,
% for both prior (1:Nx) and posterior(Nx+1:2*Nx)
[bias,stdev,rmse] = mynanstats(errKSQ.'); % outputs are 1 x 2*Nx
Nx = size(rmse,2)/2;
post = Nx+1:2*Nx;
srmse = mean(rmse(post));
sbias = mean(bias(post));
sstdev = mean(stdev(post));
fprintf('Posterior Spatial Mean of RMSE: %8.3f\n',srmse);
fprintf('Posterior Spatial Mean of Bias: %8.3f\n',sbias);
fprintf('Posterior Spatial Mean of St. Dev.: %8.3f\n',sstdev);

%% Upper error bar for MSE
%ci = 0.95;
%upper = oneside_ci(scoreKSQ.', ci);

%% Plot MSE, confidence interval, and components of MSE
%plot_rmse_ring(rmse, upper, rmse_logscale, rmse_plot_prior, bias);

%% Ensemble variance
%[ens_bias,ens_stdev,ens_rmse] = mynanstats(ensvarKSQ.');

%% Confidence interval for ensemble variance
% Need to compute effective N, which is less than N
% because of temporal autocovariance
%ens_upper = oneside_ci(ensvarKSQ.', ci);

%% Plot ensemble variance, confidence interval
%plot_ensstdev_ring(ens_rmse, ens_upper, ensstdev_logscale, ensstdev_plot_prior);

%% Functions
function [rmse_logscale,rmse_plot_prior,ensstdev_logscale,ensstdev_plot_prior] = plot_options
    prompt={'RMSE Logscale','RMSE Plot Prior', 'Ens Std Dev Logscale', 'Ens Std Dev Plot Prior'};
    name='Ring Plot Options';
    numlines=1;
    default={'0','1','0','1'};,...
    answer=inputdlg(prompt,name,numlines,default);i=1;
    % RMSE Logscale
    rmse_logscale = str2double(answer{i});i=i+1;
    % RMSE Plot Prior
    rmse_plot_prior = str2double(answer{i});i=i+1;
    % Ens Std Dev Logscale
    ensstdev_logscale = str2double(answer{i});i=i+1;
    % Ens Std Dev Plot Prior
    ensstdev_plot_prior = str2double(answer{i});
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

function plot_rmse_ring(rmse, upper, logscale, plot_prior, bias) %Needs fixed
% RMSE is 1 x 2*Nx
% First Nx values are the prior; second Nx values are the posterior
    Nx = size(rmse,2)/2;
    figure('position',[100 200 1250 500]);
    prior = 1:Nx;
    post = Nx+1:2*Nx;
    ring = 1:Nx;
    if plot_prior==1
        plot(ring,rmse(prior),'b-',ring,rmse(prior)+upper(prior),'b:');
    end
    hold on
    plot(ring,rmse(post),'r-',ring,rmse(post)+upper(post),'r:')
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
    smse = mean(rmse(post));
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

function plot_ensstdev_ring(ensvar, upper, logscale, plot_prior) %Needs fixed
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

