% Ring plots of analysis minus background (posterior minus prior) and
% analysis minus truth, with appropriate confidence intervals.
% Bill Campbell Last modified 7/26/2022

%% Nature run and DA output
if ~exist('truth','var')
   da_spinup = 56; % 2 weeks
   % Return arguments are Nt x Nx
   [truth, prior, posterior, trueparms] = get_run_info(da_spinup);
end
[Nt, Nx] = size(truth);

%% Analysis increment statistics
increment = posterior - prior; % Nt x Nx
[inc_time_mean, inc_time_std] =...
    mynanstats(increment); % Nt x Nx => 1 x Nx
inc_rho = lagged_autocovariance(increment); % Nt x Nx => 1 x Nx

inc_ci_pct = 95; % 95% confidence interval
[inc_lower, inc_upper] = ...
    confidence_interval(Nt, inc_time_mean, inc_time_std,...
                        inc_rho, inc_ci_pct);

%% Analysis minus truth statistics
state_error = posterior - truth; % Nt x Nx
[staterr_time_mean, staterr_time_std] =...
    mynanstats(state_error); % Nt x Nx => 1 x Nx
staterr_rho = lagged_autocovariance(state_error); % 1 x Nx

staterr_ci_pct = 95; % 95% confidence interval
[staterr_lower, staterr_upper] = ...
    confidence_interval(Nt, staterr_time_mean, staterr_time_std,...
                        staterr_rho, staterr_ci_pct);

%% Plots
figname = 'Time_mean_increment_plot';
plot_meandiff(inc_time_mean, inc_lower, inc_upper,...
     staterr_time_mean, staterr_lower, staterr_upper,...
     figname)


%% Nature run, prior, posterior, parameters
function [Xt, Xb, Xa, trueparms] = get_run_info(da_spinup)
    if (~exist('da_spinup','var'))
        da_spinup = 0;
    else
        da_spinup = floor(abs(da_spinup));
    end
    % Truth file e.g. L05M3_N960_K32_F15.00_I12_b10.00_c2.50_tf0.05_spinup100_tsteps10000_seed51422
    % Prior file e.g. prior_alpha_1.640_tf_0.05_R_1.00_loc_n10_nc_1000_1st_1_skip_1_Nx_960_Kp_32_seed_29418_GKSQ
    % Nature run location
    outfolder = 'C:\Users\campbell\Documents\MATLAB\Lorenz_96_model\'; % Local hard drive
    % Load truth, prior, derive posterior and parms names, load those
    [ftruth,tpath]=uigetfile([outfolder,'*L05*'],'Choose truth file:'); % truth trajectory
    truth = fullfile(tpath,ftruth);
    % Also need nature run spinup
    string = strsplit(ftruth,'_');
    spinup = str2double(string{9}(7:end));
    [fprior,ppath]=uigetfile([outfolder,'K*/prior*'],'Choose prior file:'); % prior
    prior = fullfile(ppath,fprior);
    posterior = strrep(prior,'prior','posterior');
    tic;

    fprintf('Loading prior...\n');
    load(prior,'ZKSQ');
    Xb = squeeze(mean(ZKSQ,2)).'; % Nt x Nx
    clear ZKSQ
    [pri_time, pri_loc] = size(Xb);

    fprintf('Loading posterior...\n');
    load(posterior,'XKSQ');
    Xa = squeeze(mean(XKSQ,2)).'; % Nt x Nx
    clear XKSQ
    [post_time, post_loc] = size(Xa);

    if any([pri_time, pri_loc]) ~= any([post_time, post_loc])
        error('Prior and posterior have incompatible sizes');
    end

    fprintf('Loading truth...\n');
    load(truth,'Xt');  % Ntfull x Nx
    Xt = Xt(spinup+1:spinup+post_time,:); % Nt x Nx

    fprintf('Loading parameters...\n');
    trueparmsname = strrep(prior,'prior','trueparms');
    load(trueparmsname,'trueparms');

    % Delete DA spinup period unless it is too large
    if (da_spinup < pri_time)
        Xb = Xb(da_spinup+1:end, :);
        Xa = Xa(da_spinup+1:end, :);
        Xt = Xt(da_spinup+1:end, :);
    end
    fprintf('Data read finished after %f seconds.\n',toc)
end

%% Ring plot
function plot_meandiff(tmean1,tlo1,tup1,tmean2,tlo2,tup2,figname)
    figure('position',[100 200 1650 680]);
    Nx = size(tmean1,2);
    plot(1:Nx,tmean1,'b-',1:Nx,tlo1,'b:',1:Nx,tup1,'b:',...
         1:Nx,tmean2,'r-',1:Nx,tlo2,'r:',1:Nx,tup2,'r:')
    set(gca,'fontsize',12,'fontweight','bold')
    xlim([1,Nx]);
    nticks = 15;
    gap = floor(Nx/nticks);
    xticks([1 gap:gap:Nx]);
    grid on
    xlabel('Ring')
    ylabel('Time Mean Difference')
    legend({'Xa-Xb','','','Xa-Xt','',''},'location','northwest')
    title('Time Mean Differences with 95% CI',...
          'fontsize',20,'fontweight','bold')
    mt1 = mean(tmean1);
    text(-gap,mt1,sprintf('%5.3f',mt1),'color','b',...
         'fontsize',14,'fontweight','bold')
    mt2 = mean(tmean2);
    text(-gap,mt2,sprintf('%5.3f',mt2),'color','r',...
         'fontsize',14,'fontweight','bold')
    if exist('figname','var')
        saveas(gcf,figname,'png');
    end
end