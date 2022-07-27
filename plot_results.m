function plot_results(first,current_cycle,prior_mse_norm,prior_varse_norm,...
    post_mse_norm,post_varse_norm,Ncycles,printcycle,ci,K)
persistent h old_prior_mse_norm old_prior_stderr old_post_mse_norm old_post_stderr cifac
if first==true
    h = figure;
    xlabel('Cycles');
    ylabel('Normalized Mean Squared Error');
    ttl1=sprintf('K=%d: Normalized Mean Error Estimate = (%8.6f,%8.6f)',...
        K,prior_mse_norm,post_mse_norm);
    ttl2=sprintf('(%4.1f%% 1-Sided CI)',100*ci);
    title({ttl1,ttl2});
    % One-sided is more appropriate for a squared quantity
    cifac = abs(norminv(1-ci));
    xlim([1 Ncycles]);
    ylim([0,Inf]);
    grid on
end

% Need to compute effective N, which is less than current_cycle
% because of temporal autocovariance
% Experimentally, the lag-1 autocovariance os the mean squared error
% (analysis minus truth).^2 is ~ 0.40
rho = 0.40;
factor = (1 - rho) / (1 + rho);
Neff = current_cycle * factor;
% Now compute the standard error as the stdev / sqrt(Neff)
% This is not quite correct either, as we are dealing with positive
% definite quantities, which are non-Gaussian close to zero
prior_stderr = sqrt(prior_varse_norm/Neff);
prior_upper = prior_mse_norm + cifac*prior_stderr;
post_stderr = sqrt(post_varse_norm/Neff);
post_upper = post_mse_norm + cifac*post_stderr;
da_spinup = 400;
figure(h);
hold on
if current_cycle >= da_spinup
    plot(current_cycle,prior_mse_norm,'m.',current_cycle,prior_upper,'g');
end
plot(current_cycle,post_mse_norm,'b*',current_cycle,post_upper,'r^');
ylim([0,Inf]);
%line([current_cycle,current_cycle],[mse_norm+cifac*stderr,mse_norm-cifac*stderr]);
if first==false
    times = [current_cycle-printcycle, current_cycle];
    old_prior_upper = old_prior_mse_norm + cifac*old_prior_stderr;
    if current_cycle >= da_spinup + printcycle
        plot(times,[old_prior_mse_norm,prior_mse_norm],'m-.',...
            times,[old_prior_upper,prior_upper],'m:');
    end
    old_post_upper = old_post_mse_norm + cifac*old_post_stderr;
    plot(times,[old_post_mse_norm,post_mse_norm],'b-',...
        times,[old_post_upper,post_upper],'r--');
end
ttl1=sprintf('K=%d: Normalized Posterior Mean Squared Error Estimate',K);
ttl2=sprintf('%8.6f + %8.6f (%4.1f%% 1-Sided CI)',...
    post_mse_norm,cifac*post_stderr,100*ci);
title({ttl1,ttl2});
old_prior_mse_norm = prior_mse_norm;
old_prior_stderr = prior_stderr;
old_post_mse_norm = post_mse_norm;
old_post_stderr = post_stderr;
drawnow
end 
