%% This script calculates spread-skill ratio over time

%% Load prior, posterior, truth, and obs. Choose cycle to plot
outfolder='C:\Users\knisely\Desktop\Lorenz_05\';
% prior_alpha_1.10_tf_0.05_R_1.00_loc_n40_nc_50_1st_2_skip_4_seed_29418_GKSQ.mat
[fprior,pprior]=uigetfile([outfolder,'K500\prior*'],'Prior');
[fpost,ppost]=uigetfile([outfolder,'K500\post*'],'Posterior');
[ftruth,ptruth]=uigetfile([outfolder,'*L05*'],'Truth');
[fobs,pobs]=uigetfile([outfolder,'obs*'],'Obs');

% Choose seed
%prompt = "Enter file name ";
%name='Out name';
%numlines=1;
%default={'Model#_#infl_F#_#loc'};
%answer=inputdlg(prompt,name,numlines,default);
%outname = answer{1};

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

%% MSE, variance, and ratio calculations
climnorm = 4.7820;

MSE = mean(((post_mean-truth).^2),1);
%display(mean(MSE./(climnorm.^2)))

MSE_for = zeros(Nx,K,Nt);
for t = 1:Nt
    MSE_for(:,:,t) = (bsxfun(@minus,post(:,:,t),post_mean(:,t))).^2;
end
meanvar = mean(squeeze((K/(K-1)).*(mean(MSE_for,2))),1);

spr_sk = meanvar./MSE;
plot(MSE,'r')
%hold on;
%plot(MSE./(climnorm.^2),'r--')
hold on;
plot(meanvar,'b')
hold on;
plot(spr_sk,'k')
hold on;
plot(smoothdata(spr_sk,'gaussian',20),'k--')

legend('MSE','Variance','Variance/MSE')
ylim([0 inf])
xlabel('Cycle')

%saveas(gcf,[outfolder,'\figs\spread_skill_',outname,'.fig']);
%saveas(gcf,[outfolder,'\figs\png\spread_skill_',outname,'.png']);