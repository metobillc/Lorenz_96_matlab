function [AmB,AmBsmooth,OmHB,OmHBsmooth,OmHA,OmHAsmooth,AmT,AmTsmooth,BmT,BmTsmooth,OmHT,OmHTsmooth] = diff_stats(A,B,O,T,H,spinup)
%DIFF_STATS Diffs between Analysis, Background, Obs, and Truth
% If we are interested in the variances as well, must do paired
% differences, not just mean differences, to account for correlation
% For now, just calculate the means.

%% Check dimensions of input
[Nx, ~, Ncycles] = size(A);
[Nobs] = size(O,1);
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

%% Caluculate cyclical smoothed mean differences for BC purposes
AmBhat = [AmB(Nx-10:Nx,:);AmB;AmB(1:10,:)];
OmHBhat = [OmHB(Nobs-10:Nobs,:);OmHB;OmHB(1:10,:)];
OmHAhat = [OmHA(Nobs-10:Nobs,:);OmHA;OmHA(1:10,:)];
AmThat = [AmT(Nx-10:Nx,:);AmT;AmT(1:10,:)];
BmThat = [BmT(Nx-10:Nx,:);BmT;BmT(1:10,:)];
OmHThat = [OmHT(Nobs-10:Nobs,:);OmHT;OmHT(1:10,:)];

AmBsmooth = sgolayfilt(AmBhat,9,21);
AmBsmooth = AmBsmooth(11:Nx+10,:);
OmHBsmooth = sgolayfilt(OmHBhat,9,21);
OmHBsmooth = OmHBsmooth(11:Nobs+10,:);
OmHAsmooth = sgolayfilt(OmHAhat,9,21);
OmHAsmooth = OmHAsmooth(11:Nobs+10,:);
AmTsmooth = sgolayfilt(AmThat,9,21);
AmTsmooth = AmTsmooth(11:Nx+10,:);
BmTsmooth = sgolayfilt(BmThat,9,21);
BmTsmooth = BmTsmooth(11:Nx+10,:);
OmHTsmooth = sgolayfilt(OmHThat,9,21);
OmHTsmooth = OmHTsmooth(11:Nobs+10,:);

end
