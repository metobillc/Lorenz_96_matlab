function plot_mse_ring(tmse_prior,tmse_post,tstd_prior,tstd_post, ...
    tmean_ensvar_prior,tmean_ensvar_post,tstd_ensvar_prior, ...
    tstd_ensvar_post,ci,Ncycles)

% Plot ring plot of tmse_prior&post with error bars based on tstd_prior&post, 
% and spatial_mean_err,std in the title string.
% 
% Also ring plot of tmean_ensvar_prior&post with error bars based on 
% tstd_ensvar_prior&post and spatial_mean_ensvar,ensvar_std in the 
% title string

Nx = size(tmse_post,2);

% Find std. err., establish error bars, and take some means
tstd_prior = tstd_prior/sqrt(Ncycles);
tstd_post = tstd_post/sqrt(Ncycles);
tstd_ensvar_prior = tstd_ensvar_prior/sqrt(Ncycles);
tstd_ensvar_post = tstd_ensvar_post/sqrt(Ncycles);

cifac = abs(norminv(1-ci));
prior_upper = tmse_prior + cifac*tstd_prior;
post_upper = tmse_post + cifac*tstd_post;
ensvar_prior_upper = tmean_ensvar_prior + cifac*tstd_ensvar_prior;
ensvar_post_upper = tmean_ensvar_post + cifac*tstd_ensvar_post;
prior_lower = tmse_prior - cifac*tstd_prior;
post_lower = tmse_post - cifac*tstd_post;
ensvar_prior_lower = tmean_ensvar_prior - cifac*tstd_ensvar_prior;
ensvar_post_lower = tmean_ensvar_post - cifac*tstd_ensvar_post;

spatial_mean_err = mean(tmse_post);
spatial_mean_stderr = mean(tstd_post);
spatial_mean_ensvar = mean(tmean_ensvar_post);
spatial_mean_ensvar_stderr = mean(tstd_ensvar_post);

%% Plot

figure('position',[100 200 1100 450]);
plot(1:Nx,tmse_prior,'b-',1:Nx,prior_upper,'b:',1:Nx,prior_lower,'b:', ...
     1:Nx,tmse_post,'r-',1:Nx,post_upper,'r:',1:Nx,post_lower,'r:')
set(gca,'fontsize',12,'fontweight','bold')
xlim([1,Nx]);
ylim([0,0.5]); %inf or 0.5
nticks = 15;
gap = floor(Nx/nticks);
xticks([1 gap:gap:Nx]);
grid on
xlabel('Ring')
ylabel('Time MSE')
legend({'Prior','','','Post','',''},'location','northwest')
ttl1 = sprintf('Posterior Mean MSE: %8.5f',spatial_mean_err);
ttl2 = sprintf('Std. Err.: %8.5f',spatial_mean_stderr);
title({ttl1,ttl2})

figure('position',[100 200 1100 450]);
plot(1:Nx,tmean_ensvar_prior,'b-',1:Nx,ensvar_prior_upper,'b:',1:Nx,ensvar_prior_lower,'b:', ...
     1:Nx,tmean_ensvar_post,'r-',1:Nx,ensvar_post_upper,'r:',1:Nx,ensvar_post_lower,'r:')
set(gca,'fontsize',12,'fontweight','bold')
xlim([1,Nx]);
ylim([0,0.5]); %inf or 0.5
nticks = 15;
gap = floor(Nx/nticks);
xticks([1 gap:gap:Nx]);
grid on
xlabel('Ring')
ylabel('Time Mean Ens. Var.')
legend({'Prior','','','Post','',''},'location','northwest')
ttl1 = sprintf('Post. Mean Ensemble Variance: %8.5f',spatial_mean_ensvar);
ttl2 = sprintf('Std. Err.: %8.5f',spatial_mean_ensvar_stderr);
title({ttl1,ttl2})

end