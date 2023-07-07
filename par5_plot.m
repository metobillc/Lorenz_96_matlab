% This script is for data exploration
% Bill Campbell
% Last modified 7/21/2022

%% Load prior, posterior, truth, and obs
outfolder='C:\Users\campbell\Documents\MATLAB\Lorenz_96_model\';
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

%% Forecast and analysis error covariance matrices
% at the final time
num=20;
skm1i=1./sqrt(num)./sqrt(K-1);
Pfsamp = zeros(Nx);
Pasamp = Pfsamp;
for t=Nt-num+1:Nt
    Zt = squeeze(prior(:,:,t));  % Space x Nmem
    Ztbar = repmat(mean(Zt,2),1,K);
    Zt = skm1i*(Zt - Ztbar);
    Pfsamp = Pfsamp + Zt * Zt';
    
    Xt = squeeze(post(:,:,t));
    Xtbar = repmat(mean(Xt,2),1,K);
    Xt = skm1i*(Xt - Xtbar);
    Pasamp = Pasamp + Xt * Xt';
end
CL = covloc(Nx,1,Nx/20);
Pfloc = CL.*Pfsamp;
Paloc = CL.*Pasamp;

figure (1);
imagesc(Pfsamp);
axis square;
colormap linspecer;
colorbar('v');
ttl=sprintf('K=%d Sample Forecast Error Covariance\n',K);
title(ttl);

figure (2);
imagesc(Pasamp);
axis square;
colormap linspecer;
colorbar('v');
ttl=sprintf('K=%d Sample Analysis Error Covariance\n',K);
title(ttl);

figure (3);
imagesc(Pfloc);
axis square;
colormap linspecer;
colorbar('v');
ttl=sprintf('K=%d Localized Sample Forecast Error Covariance\n',K);
title(ttl);

figure (4);
imagesc(Paloc);
axis square;
colormap linspecer;
colorbar('v');
ttl=sprintf('K=%d Localized Sample Analysis Error Covariance\n',K);
title(ttl);

%% Time series plots
% figure (5);
% for t=1:ncycles,
%     clf
%     hold on
%     mp1=plot(1:ngrid,truth(:,t),'k-',first:skip:ngrid,obs(first:skip:ngrid,t),'bo',...
%         1:ngrid,prior_mean(:,t),'r--',1:ngrid,post_mean(:,t),'g--');
%     set(mp1,'linewidth',2);
%     mp2=patch([1:ngrid ngrid:-1:1],[prior_min(:,t)' prior_max(ngrid:-1:1,t)'],...  
%         'y','facealpha',0.5,'edgecolor','y');
%     grid on
%     ylabel('Amplitude');
%     xlabel('Gridpoints');
%     set(gca,'xlim',[-10 ngrid+10],'ylim',[-15 20],'linewidth',1.5);
%     ttlstr=sprintf('State vector at t=%d\n',t);
%     title(ttlstr);
%     legend('Truth','Obs','Prior Mean','Post Mean','Prior Envelope');
%     hold off;
%     pause;
% end

