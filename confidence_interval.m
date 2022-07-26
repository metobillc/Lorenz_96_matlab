%% Confidence interval accounting for autocorrelation
function [lower, upper] =...
    confidence_interval(Nt, tmean, tstd, trho, ci_pct)
   if (~exist('ci_pct','var'))
       ci_pct = 95; % 95 percent confidence interval
   else
       ci_pct = abs(ci_pct);
       if ci_pct < 1.0
           ci_pct = 100* ci_pct;
       end
       if ci_pct >= 100.0
           fprintf('Confidence interval percent must be < 100%...'+...
               'defaulting to 95%.\n');
           ci_pct = 95;
       end
   end
   Nx = size(tmean,2);
   if (size(tstd,2) ~= Nx || size(trho,2) ~= Nx)
       error('Mean, stddev, and rho all must have Nx columns, aborting...')
   end
   cifac = norminv(ci_pct/100.);
   factor = (1 - trho) ./ (1 + trho); % 1 x Nx
   Neff = Nt * factor; % 1 x Nx
   stderr = tstd ./ sqrt(Neff); % 1 x Nx
   lower = tmean - cifac*stderr; % 1 x Nx
   upper = tmean + cifac*stderr; % 1 x Nx
end