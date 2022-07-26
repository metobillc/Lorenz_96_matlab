%% Lagged autocovariance
function rho = lagged_autocovariance(tseries, lag)
    % Expect tseries to have Nt rows and 1 or more columns
    if ~exist('lag','var')
        lag = 1;
    end
    % Each column is a time series
    Nt = size(tseries,1);
    if Nt == 1
        tseries = tseries.';
        Nt = size(tseries,1);
    end
    if Nt < 2
        fprintf('Time series must have at least 2 times, returning 0');
        rho = zeros(size(tseries));
        return
    end
    tplus = tseries(1+lag:end,:); % (Nt-lag) x Nx
    [tplus_mean, tplus_std] = mynanstats(tplus); % 1 x Nx
    tplus_prime = tplus - tplus_mean; % (Nt-lag) x Nx
    tplus_pn = tplus_prime ./ tplus_std; % (Nt-lag) x Nx

    tminus = tseries(1:end-lag,:); % (Nt-lag) x Nx
    [tminus_mean, tminus_std] = mynanstats(tminus); % 1 x Nx
    tminus_prime = tminus - tminus_mean; % (Nt-lag) x Nx
    tminus_pn = tminus_prime ./ tminus_std; % (Nt-lag) x Nx
    
    rho = mean(tplus_pn.*tminus_pn); % 1 x Nx
end