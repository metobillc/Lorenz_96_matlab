function [Xa,score,ensvar] = KSQ(Zfull,alpha,K,H,R,y,Xt,varargin)
% Liz Satterfield 10/24/2019
% Add localization (new parameter)
% Return prior and posterior error variances along with squared errors
%% Ensemble Transform Kalman Filter


%% Localization (new section)
xb_bar = mean(Zfull,2); % Nx x 1
Xb=sqrt(alpha)*(Zfull - repmat(xb_bar,1,K)); % Nx x K
Nx=length(xb_bar);
Rinv=R^(-1); %Nobs x Nobs
Rsqinv = diag(1./sqrt(diag(R))); % Nobs x Nobs
Htilde=Rsqinv*H; %Nobs x Nobs
CL = varargin{1};
Yb=H*Xb; %Nobs x K
yb_bar=H*xb_bar; %Nobs x 1
Zf=Xb/sqrt(K-1); %Nx x K
Pf=Zf*Zf'; %Nx x Nx

%update ensemble
A=Htilde*Zf; 
% K x K
[C Gamma]=eig(A'*A);
T = C*sqrt(pinv(Gamma + eye(K)));
Pa=(Zf*T)*(Zf*T)';
Za=Zf*T*C';
Xa=Za*sqrt(K-1);

%update mean
%traditional KF update
%xa_bar=xb_bar+Pf*H'*inv(H*(Pf)*H'+R)*(y-yb_bar);
%traditional KF update with localization
xa_bar=xb_bar+(CL.*Pf)*H'*inv(H*(CL.*Pf)*H'+R)*(y-yb_bar);
%ETKF update xa_bar=xb_bar+Pa*H'*Rinv*(y-yb_bar)
% xa_bar=xb_bar+Zf*C*pinv(Gamma+eye(K))*C'*Zf'*Htilde'*Rsqinv*(y-yb_bar);

% Recenter Analysis ensemble
Xa = Xa + repmat(xa_bar,1,K); % Nx x K

%% MSE (score) of prior and posterior means
priscore = (xb_bar - Xt).^2; % Nx x 1
postscore = (xa_bar - Xt).^2; % Nx x 1
score = [priscore; postscore]; % 2*Nx x 1

privar = diag(Pf);
postvar = diag(Pa);
ensvar = [privar; postvar]; % 2*Nx x 1
% 
% mean(var(Xa'))
% mean(privar)
% mean(priscore)
% mean(postvar)
% mean(postscore)

%%%%%%%%%%%%%plotting
% figure(1);
% obs_plot=sum(H);
% ix=find(obs_plot==1);
% ih=find(obs_plot==0);
% if ~isempty(ih)
%     obs_plot(ih)=nan
% end
% obs_plot(ix)=y;
% h=plot(Xt(:),'k-');
% set(h,'LineWidth',2)
% hold on
% axis([1 Nx -20 20]);
% for j=1:K
%     plot(Xa(:,j),'b--');
% end
% plot(obs_plot,'r*');
% hold off
% pause(0.01)
end
