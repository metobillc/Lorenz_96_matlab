function plot_results(first,ncycle,aerr,variance,Ncycles,printcycle,ci,K)
persistent h old_aerr old_stderr cifac
if first==true,
    h = figure;
    xlabel('Cycles');
    ylabel('Mean Error');
    ttlstr=sprintf('K=%d: Mean Error Estimate = %6.4f with %4.1f%% Confidence Intervals',K,aerr,100*ci);
    title({ttlstr});
    cifac = abs(norminv((1-ci)/2));
    set(gca,'xlim',[1 Ncycles]);
    grid on
end
figure(h);
hold on
stderr = sqrt(variance/ncycle);
plot(ncycle,aerr,'b.',ncycle,aerr+cifac*stderr,'r^',ncycle,aerr-cifac*stderr,'rv');
%line([ncycle,ncycle],[aerr+cifac*stderr,aerr-cifac*stderr]);
if first==false,
    plot([ncycle-printcycle,ncycle],[old_aerr,aerr],'b-',...
        [ncycle-printcycle,ncycle],[old_aerr+cifac*old_stderr,aerr+cifac*stderr],'r--',...
        [ncycle-printcycle,ncycle],[old_aerr-cifac*old_stderr,aerr-cifac*stderr],'r--');
end
ttl1=sprintf('K=%d: Mean Error Estimate = %8.6f +/- %8.6f',K,aerr,cifac*stderr);
ttl2=sprintf('(%4.1f%% Confidence Intervals)',100*ci);
title({ttl1,ttl2});
old_aerr = aerr;
old_stderr = stderr;
drawnow
end 
