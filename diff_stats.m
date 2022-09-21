function [AmB,OmHB,OmHA,AmT,BmT,OmHT] = diff_stats(A,B,O,T,H,spinup)
%DIFF_STATS Diffs between Analysis, Background, Obs, and Truth
% If we are interested in the variances as well, must do paired
% differences, not just mean differences, to account for correlation
% For now, just calculate the means.

%% Check dimensions of input
[Nx, ~, Ncycles] = size(A);
if ~isequal(size(B), size(A))
    error("A and B must have same size, aborting...");
end
if ~isequal(size(T), [Nx, Ncycles])
    error("T and ensmean of A must have same size, aborting...");
end
if ~isequal(size(O), size(H*T))
    error("O and H*T must have same size, aborting...");
end
if Ncycles - spinup < 2
    error("Ncycles must be at least 2 greater than spinup, aborting...");
end

%% Compute means discarding first spinup cycles
ensmA = squeeze(mean(A(:,:,spinup+1:end),2)); % Nx x Ncycles-spinup
ensmB = squeeze(mean(B(:,:,spinup+1:end),2)); % Nx x Ncycles-spinup
O = O(:,spinup+1:end); % Nobs x Ncycles-spinup
T = T(:,spinup+1:end); % Nx x Ncycles-spinup

%% Time mean differences as a function of position 
AmB = mean(ensmA - ensmB, 2); % Nx x 1
OmHB = mean(O - H*ensmB,2); % Nobs x 1
OmHA = mean(O - H*ensmA,2); % Nobs x 1
AmT = mean(ensmA - T, 2); % Nx x 1
BmT = mean(ensmB - T,2); % Nx x 1
OmHT = mean(O - H*T,2); % Nobs x 1

end
