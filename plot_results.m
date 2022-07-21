function plot_results(first,current_cycle,prior_mse_norm,prior_varse_norm,...
    post_mse_norm,post_varse_norm,Ncycles,printcycle,ci,K)
persistent h old_prior_mse_norm old_prior_stderr old_post_mse_norm old_post_stderr cifac
if first==true
    h = figure;
    xlabel('Cycles');
    ylabel('Normalized Mean Squared Error');
    ttl1=sprintf('K=%d: Normalized Mean Error Estimate = (%8.6f,%8.6f)',...
        K,prior_mse_norm,post_mse_norm);
    ttl2=sprintf('(%4.1f%% Confidence Intervals)',100*ci);
    title({ttl1,ttl2});
    cifac = abs(norminv((1-ci)/2));
    set(gca,'xlim',[1 Ncycles]);
    grid on
end
figure(h);
hold on
prior_stderr = sqrt(prior_varse_norm/current_cycle);
prior_upper = prior_mse_norm + cifac*prior_stderr;
prior_lower = prior_mse_norm - cifac*prior_stderr;
post_stderr = sqrt(post_varse_norm/current_cycle);
post_upper = post_mse_norm + cifac*post_stderr;
post_lower = post_mse_norm - cifac*post_stderr;
da_spinup = 300;
if current_cycle >= da_spinup
    plot(current_cycle,prior_mse_norm,'m.',...
        current_cycle,prior_upper,'g',current_cycle,prior_lower,'g');
end
plot(current_cycle,post_mse_norm,'b*',...
    current_cycle,post_upper,'r^',current_cycle,post_lower,'rv');
%line([current_cycle,current_cycle],[mse_norm+cifac*stderr,mse_norm-cifac*stderr]);
if first==false
    times = [current_cycle-printcycle, current_cycle];
    old_prior_upper = old_prior_mse_norm + cifac*old_prior_stderr;
    old_prior_lower = old_prior_mse_norm - cifac*old_prior_stderr;
    if current_cycle >= da_spinup + printcycle
        plot(times,[old_prior_mse_norm,prior_mse_norm],'m-.',...
            times,[old_prior_upper,prior_upper],'m:',...
            times,[old_prior_lower,prior_lower],'m:');
    end
    old_post_upper = old_post_mse_norm + cifac*old_post_stderr;
    old_post_lower = old_post_mse_norm - cifac*old_post_stderr;
    plot(times,[old_post_mse_norm,post_mse_norm],'b-',...
        times,[old_post_upper,post_upper],'r--',...
        times,[old_post_lower,post_lower],'r--');
end
ttl1=sprintf('K=%d: Normalized Posterior Mean Squared Error Estimate',K);
ttl2=sprintf('%8.6f +/- %8.6f (%4.1f%% Confidence Intervals)',...
    post_mse_norm,cifac*post_stderr,100*ci);
title({ttl1,ttl2});
old_prior_mse_norm = prior_mse_norm;
old_prior_stderr = prior_stderr;
old_post_mse_norm = post_mse_norm;
old_post_stderr = post_stderr;
drawnow
end 
