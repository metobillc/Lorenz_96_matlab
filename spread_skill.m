%% This script calculates spread-skill ratio over time

%% Load prior, posterior, truth, and obs. Choose cycle to plot
outfolder='C:\Users\knisely\Desktop\Lorenz_96\';
% prior_alpha_1.10_tf_0.05_R_1.00_loc_n40_nc_50_1st_2_skip_4_seed_29418_GKSQ.mat
[fprior,pprior]=uigetfile([outfolder,'K500\prior*'],'Prior');
[fpost,ppost]=uigetfile([outfolder,'K500\post*'],'Posterior');
[ftruth,ptruth]=uigetfile([outfolder,'*L05*'],'Truth');
[fobs,pobs]=uigetfile([outfolder,'obs*'],'Obs');

prior=load_ensemble(pprior,fprior);  
[Nx,K,Nt]=size(prior);  % Space x Nmem x Time
post=load_ensemble(ppost,fpost);
truth=load_truth(ptruth,ftruth);
truth=truth(:,1:Nt); % Space x Time
obs=load_obs(pobs,fobs);
obs=obs(:,1:Nt); % Space x Time
% Skip unobserved gridpoints
string = strsplit(fprior,'_');
if length(string)>11
    first = str2num(string{13});
    skip = str2num(string{15});
else
    first = 1;
    skip = 1;
end

%% Ensemble means, min, and max (Space x Time)
prior_mean=squeeze(mean(prior,2)); 
prior_min = squeeze(min(prior,[],2));
prior_max = squeeze(max(prior,[],2));

post_mean=squeeze(mean(post,2));
post_min = squeeze(min(post,[],2));
post_max = squeeze(max(post,[],2));

%% RMSE, spread, and ratio calculations

RMSE = sqrt(mean(((obs-post_mean).^2),1));

% RMS of post, post_mean
%RMS_p = squeeze(sqrt(mean((bsxfun(@minus,post,mean(post,1))).^2,1)));
%RMS_pm = sqrt(mean((bsxfun(@minus,post_mean,mean(post_mean,1))).^2,1));
% Spread calc
%spread = sqrt((K/(K-1))*mean((bsxfun(@plus,RMS_p,RMS_pm)).^2,1));

MSE_for = zeros(Nx,K,Nt);
for t = 1:Nt
    MSE_for(:,:,t) = (bsxfun(@minus,post(:,:,t),post_mean(:,t))).^2;
end
spread = sqrt(mean(squeeze((K/(K-1)).*(mean(MSE_for,2))),1));

spr_sk = spread./RMSE;
plot(RMSE,'r')
hold on;
plot(spread,'b')
hold on;
plot(spr_sk,'k')
legend('RMSE','Spread','Spread/RMSE')
ylim([0 2])
xlabel('Cycle')