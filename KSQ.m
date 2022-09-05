function [Xa,score,ensvar,varargout] = KSQ(Zfull,alpha,K,H,R,y,Xt,varargin)
% Liz Satterfield 10/24/2019
% Add localization (new parameter)
% Return prior and posterior error variances along with squared errors
% Model and obs bias correction by Bill Campbell 9/4/2022

if (nargin > 7)
    CL = varargin{1};
end
if (nargin == 8)
    obs_parms = varargin{2};
    bc_obs = obs_parms.bc_obs;
    bc_simobs = obs_parms.bc_simobs;
else
    bc_obs = 0;
    bc_simobs = 0;
end
%% Ensemble Transform Kalman Filter

xb_bar = mean(Zfull,2); % Nx x 1
Xb = sqrt(alpha).*(Zfull - repmat(xb_bar,1,K)); % Nx x K

%% Model and obs bias correction
% Apply "model bias" correction only to the simulated obs,
% NOT to the background itself, i.e. leave xb_bar unchanged
if (bc_simobs == 1)
   xb_bar_bc = apply_model_bias_correction_to_obs(xb_bar);
   yb_bar  = H * xb_bar_bc; % Nobs x 1
else
    yb_bar = H * xb_bar; % Nobs x 1
end
% Apply obs bias correction directly to obs
if (bc_obs == 1) 
    y = apply_obs_bias_correction(y);
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
repr = H*PflocHt + R;
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

function xb_bar_bc = apply_model_bias_correction_to_obs(xb_bar)
    % Can e.g. add long-term mean of analysis minus background here
    anal_minus_bk = 0.0; % Instead read from a file
    % Add to background the mean of A-B from previous model bias run with
    % perfect(i.e. unbiased) obs and no bias correction
    xb_bar_bc = xb_bar + anal_minus_bk;
end

function y_bc = apply_obs_bias_correction(y)
    obs_minus_bk = 0.0; % Instead read from file
    % Subtract from ob the mean of O-B from previous obs bias run with
    % perfect(i.e. identical model params e.g. forcing) model and no bias correction 
    y_bc = y - obs_minus_bk;
end
