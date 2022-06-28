% Nature_plot
% Plot time series of Lorenz 96 Model II nature runs and observations
% Bill Campbell
% Last modified 9/1/16
%% Get nature run file
outfolder = 'F:\Matlab\work\MultimodelEnsembles\*.mat'; % Local hard drive
[fname,pname]=uigetfile(outfolder,'Load Lorenz 96 Model II full nature run:'); % truth trajectory
load([pname fname]);
[nT,N]=size(Xt);
string=strsplit(fname,'_');
K=str2num(string{3}(2:end));
[xtbar,xtstd]=mynanstats(Xt(:));
[fname,pname]=uigetfile(outfolder,'Load Lorenz 96 Model II nature run observations:'); % truth trajectory
load([pname fname]);
% Should have observed every variable at every timestep -- can subsample
% (thin) later if desired
% yobs is N by Nt, want it the same dimensions as Xt
yobs = yobs';

%% Animate gridpoints in time
figure;
plot(1:N,Xt(1,:),'b-',1:N,yobs(1,:),'ro');
set(gca,'ylim',[xtbar-4*xtstd,xtbar+4*xtstd],'xlim',[1 N],'linewidth',2);
grid on
ttlstr=sprintf('Nature Run Animation N=%d, K=%d\n',N,K);
title(ttlstr);
ylabel('Amplitude');
xlabel('Gridpoint');
disp('Hit any key to begin animation');
pause;
for t=2:nT,
    plot(1:N,Xt(t,:),'b-',1:N,yobs(t,:),'ro');
    set(gca,'ylim',[xtbar-4*xtstd,xtbar+4*xtstd],'xlim',[1 N]);
    grid on
    title(ttlstr);
    ylabel('Amplitude');
    xlabel('Gridpoint');
    drawnow
%    pause(0.05);
end
