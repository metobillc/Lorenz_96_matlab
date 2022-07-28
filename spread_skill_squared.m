% My take on spread skill
% This script calculates spread-skill ratio over time

%% Load prior, posterior, truth
outfolder='C:\Users\Campbell\Documents\MATLAB\Lorenz_96\';
% Prefix is the same:
% prior_alpha_1.640_tf_0.05_R_1.00_loc_n10_nc_10000_1st_1_skip_1_
% Old suffix: seed_29418_GKSQ.mat
% New suffix: Nx_960_Kp_32_seed_29418_GKSQ.mat

da_spinup = 56; % 14 days = 40 cycles
if ~exist('post','var')
    [ftruth,ptruth]=uigetfile([outfolder,'*L05*'],'Truth');
    [fpost,ppost]=uigetfile([outfolder,'K500\posterior*'],'Posterior');
    post=load_ensemble(ppost,fpost);
    post = post(:,:,1+da_spinup:end);
    % Convenient to permute indices
    post = permute(post,[1,3,2]);
    [Nx,Nt,K]=size(post);  % Space x Time x Nmem
    truth=load_truth(ptruth,ftruth);
    string = strsplit(ftruth,'_');
    % Naming convention has 9th field as nature run spinup
    spinup = str2double(string{9}(7:end));
    allspin = spinup + da_spinup;
    truth=truth(:,1+allspin:Nt+allspin); % Space x Time, commensurate with post
end

%% MSE, spread^2, and ratio calculations
% Skill is measured by MSE
% Analysis is the posterior ensemble mean (space x time)
% Truth is truth (space x time)
post_ensmean = squeeze(mean(post,3)); % space x time
SE = (post_ensmean - truth).^2; % space x time
TMSE = mean(SE,2); % time mean is space x 1
SMSE = mean(SE,1); % spatial mean is 1 x time

% Skill is measured by the time mean of the ensemble mean
% of the ensembele variance
ensmean = mean(post,3); % space x time x 1
ensmeanfull = repmat(ensmean,1,1,K); % space x time x member
enssqdif = (post - ensmeanfull).^2; % space x time x member
ensvar = mean(enssqdif,3); % ensemble mean is space x time
tmean_ensvar = mean(ensvar,2); % time mean is space x 1
smean_ensvar = mean(ensvar,1); % spatial mean is 1 x time

% Spread-skill ratios
time_spread_skill_ratio = tmean_ensvar./TMSE; % space x 1
space_spread_skill_ratio = smean_ensvar./SMSE; % 1 x time

%% Ring plot
figure;
ax1 = nexttile;
xx = 1:Nx;
plot(xx,TMSE,'r');
hold on;
plot(xx,tmean_ensvar,'b')
hold on
grid on
plot(xx,time_spread_skill_ratio,'k')
legend('TMSE','Spread','Spread/TMSE','location','best')
xlim([1 Nx])
xtix = xticks;
if xtix(1) > 1
    xtix = [1 xtix];
end
if xtix(end) < Nx
    xtix = [xtix Nx];
end
set(gca,'xtick',xtix,'yscale','log')
xlabel('Ring')
title('Spread-Skill at Each Gridpoint')

%% Time series plot
figure;
ax2 = nexttile;
tt=da_spinup+1:da_spinup+Nt;
plot(tt,SMSE,'r');
hold on
grid on
plot(tt,smean_ensvar,'b')
plot(tt,space_spread_skill_ratio,'k')
legend('SMSE','Spread','Spread/SMSE','location','best')
xlim([tt(1) tt(end)])
ttix = xticks;
if ttix(1) > da_spinup + 1
    ttix = [(da_spinup + 1) ttix];
end
if ttix(end) < da_spinup + Nt
    ttix = [ttix da_spinup + Nt];
end
set(gca,'xtick',ttix,'yscale','log')
xlabel('Cycle')
title('Spread-Skill at Each Timestep')


%% Synchronize y-axes of the plots
miny = min(ax1.YLim, ax2.YLim);
maxy = max(ax1.YLim, ax2.YLim);
ax1.YLim = [miny,maxy];
linkaxes([ax1,ax2],'y')
