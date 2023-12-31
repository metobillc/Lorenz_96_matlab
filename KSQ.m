function [Xa,score,ensvar,varargout] = KSQ(Zfull,alpha,K,H,R,y,Xt,varargin)
% Liz Satterfield 10/24/2019
% Add localization (new parameter)
% Return prior and posterior error variances along with squared errors
% Model and obs bias correction by Bill Campbell 6/26/2023

if nargin > 7
    CL = varargin{1};
end
if nargin > 8
    biascor = varargin{2};
else
    biascor.apply_to_standard_obs = false;
    biascor.apply_to_simobs_only = false;
    biascor.apply_to_model = false;
end
if nargin > 9
    obs = varargin{3};
    std_idx = obs.standard_idx;
else
    std_idx = 1:size(y, 1); % y is Nobs x Ncycles
end

%% Model bias correction
if biascor.apply_to_model
    % Add smoothed mean of A - B from another run to model state
    Zfull = Zfull + repmat(biascor.AmB_smooth_post,1,K);
end

%% Ensemble Transform Kalman Filter
xb_bar = mean(Zfull,2); % Nx x 1
Xb = sqrt(alpha).*(Zfull - repmat(xb_bar,1,K)); % Nx x K

%% Bias correction for obs and simulated obs
if biascor.apply_to_simobs_only && ~biascor.apply_to_model  % Do not correct twice if bc_model true
    yb_bar  = H * (xb_bar + biascor.AmB_smooth_post); % Nobs x 1
else
    yb_bar = H * xb_bar; % Nobs x 1
end
% Apply obs bias correction directly to obs
if biascor.apply_to_standard_obs  % Subtract mean of O - H*B from another run
    % Multiplicative and additive obs biases
    y(std_idx,:) = y(std_idx,:) - biascor.OmHB_smooth_post(std_idx);
end

%% Create and update ensemble
skm1 = sqrt(K-1);
Zf = Xb ./ skm1; % Nx x K
Pf = Zf * Zf.'; % Nx x Nx

% Only handles diagonal R for now
Rsqinv = diag(1./sqrt(diag(R))); % Nobs x Nobs
Htilde = Rsqinv * H; % Nobs x Nobs
A = Htilde * Zf;
% K x K
[C, Gamma] = eig(A.' * A);
T = C * sqrt(pinv(Gamma + eye(K)));
zft = Zf * T;
Pa = zft * zft.';
Za = zft * C.';
Xa = Za .* skm1;

%update mean
%traditional KF update
%xa_bar=xb_bar+Pf*H'*inv(H*(Pf)*H'+R)*(y-yb_bar);
%traditional KF update with localization
Pfloc = CL .* Pf;
PflocHt = Pfloc * H.';
repr = H * PflocHt + R;
W = repr \ (y - yb_bar); % More stable than inv(repr)*(y - yb_bar)
xa_bar = xb_bar + PflocHt * W;
%ETKF update xa_bar=xb_bar+Pa*H'*Rinv*(y - yb_bar)
% xa_bar=xb_bar+Zf*C*pinv(Gamma+eye(K))*C'*Zf'*Htilde'*Rsqinv*(y - yb_bar);

% Recenter Analysis ensemble
Xa = Xa + repmat(xa_bar,1,K); % Nx x K

%% MSE (score) of prior and posterior ensemble means
prierr = xb_bar - Xt;
priscore = prierr.^2; % Nx x 1
posterr = xa_bar - Xt;
postscore = posterr.^2; % Nx x 1
score = [priscore; postscore]; % 2*Nx x 1
err = [prierr; posterr];
if nargout==4
    varargout{1} = err;
end

%% Prior and posterior ensemble variances
privar = diag(Pf);
postvar = diag(Pa);
ensvar = [privar; postvar]; % 2*Nx x 1
end
