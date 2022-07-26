%% Confidence interval accounting for autocorrelation
function [lower, upper, varargout] =...
    confidence_interval(Nt, tmean, tstd, trho, ci_pct)
   % tmean, tstd, trho expected to be 1 x Nx
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
   ci = ci_pct/100;
   % Two-sided CI
   cifac = abs(norminv((1 - ci)/2));
   factor = (1 - trho) ./ (1 + trho); % 1 x Nx
   Neff = Nt * factor; % 1 x Nx
   stderr = tstd ./ sqrt(Neff); % 1 x Nx
   lower = tmean - cifac*stderr; % 1 x Nx
   upper = tmean + cifac*stderr; % 1 x Nx
   if nargout >= 3
       varargout{1} = Neff;
   end
   if nargout >= 4
       varargout{2} = cifac;
   end
   if nargout >= 5
       varargout{3} = stderr;
   end
end