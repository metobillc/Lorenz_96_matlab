%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute background error covariance from
% true trajectory
% Bill Campbell
% Last modified 6/28/2022

%% Get nature run
needed = {'Xt','spinup','tF'};
reload = false;
for i = 1:length(needed)
    if ~exist(needed{i},'var')
        reload = true;
        fprintf('%s variable not present\n',needed{i})
        break
    end
end
if reload
    outfolder = 'C:\Users\knisely\Desktop\Lorenz_05\'; % Local hard drive
    [ftruth,pathname]=uigetfile([outfolder,'*L05*'],...
                      'Load Lorenz 96 Model I, II or III full nature run:'); % truth trajectory
    % Extract spinup from truth file
    string=strsplit(ftruth,'_');
    idx=8; tF = str2double(string{idx}(3:end));
    idx=9; spinup = str2double(string{idx}(7:end));
    tic;
    [Xt,abstol,reltol] = load_truth(pathname,ftruth);
    toc;
end
tic;
[N,Ncycles] = size(Xt);

%% Compute time mean state and perturbation
sqncm1 = sqrt(Ncycles - 1);
Zt = Xt./sqncm1;
% tF==1 corresponds to 5 days, so 1 day is tF==0.20
oneday = 0.20;
shift1dy = oneday/tF;
shift6hr = shift1dy/4;
diff6hr = Zt(:,(spinup+shift6hr):Ncycles) - Zt(:,spinup:(Ncycles-shift6hr));
diff1dy = Zt(:,(spinup+shift1dy):Ncycles) - Zt(:,spinup:(Ncycles-shift1dy));
Zts = Zt(:,(spinup+1):Ncycles);
Ztsbar = mean(Zts,2);
Zt6bar = mean(diff6hr,2);
Zt24bar = mean(diff1dy,2);
Ztsprime = Zts - repmat(Ztsbar,1,Ncycles-spinup);
Zt6prime = diff6hr - repmat(Zt6bar,1,(Ncycles-spinup+1-shift6hr));
Zt24prime = diff1dy - repmat(Zt24bar,1,(Ncycles-spinup+1-shift1dy));

%% Compute climatological background error covariance
B6 = Zt6prime * Zt6prime.';
B24 = Zt24prime * Zt24prime.';
Bts = Ztsprime * Ztsprime.';
toc

%% Plot
plotit = true;
if plotit
    figure('position',[900,400,1500,500]);
    subplot(131);
    im6 = imagesc(flipud(B6));
    axis square;
    colorbar('v');
    colormap linspecer;
    %set(gca,'clim', [-0.10, 1]);
    wavenum_est = 6;
    dist = floor(N/wavenum_est);
    tix = [1 dist:dist:(N-dist) N];
    set(gca,'xtick',tix);
    set(gca,'ytick',tix);
    set(gca, 'YTickLabel',flip(tix))
    grid on
    set(gca,'fontweight', 'bold');
    set(gca,'fontsize', 14);
    title('Nature Run 6hr B');
    subplot(132);
    im24 = imagesc(flipud(B24));
    axis square;
    colorbar('v');
    colormap linspecer;
    %set(gca,'clim', [-0.10, 1]);
    set(gca,'xtick',tix);
    set(gca,'ytick',tix);
    set(gca, 'YTickLabel',flip(tix))
    grid on
    set(gca,'fontweight', 'bold');
    set(gca,'fontsize', 14);
    title('Nature Run 1dy B');
    subplot(133);
    ims = imagesc(flipud(Bts));
    axis square;
    colorbar('v');
    colormap linspecer;
    %set(gca,'clim', [-0.10, 1]);
    set(gca,'xtick',tix);
    set(gca,'ytick',tix);
    set(gca, 'YTickLabel',flip(tix))
    grid on
    set(gca,'fontweight', 'bold');
    set(gca,'fontsize', 14);
    title('Nature Run Climo B');
    % Set 1dy to same colorscale as climatological
    set(get(im24,'Parent'),'clim',get(get(ims,'Parent'),'clim'));
    sgtitle('Background Error Covariance Matrices',...
        fontsize=20,fontweight='Bold')
end
