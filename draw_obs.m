%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw_obs from true trajectory, with option
% for a state-dependent bias and error of 
% representation via spatial convolution
% Bill Campbell
% Last modified 7/29/2022

%% Get nature run
outfolder = 'C:\Users\knisely\Desktop\Lorenz_05\'; % Local hard drive
[ftruth,pathname]=uigetfile([outfolder,'*L05*'],...
                  'Load Lorenz 96 Model I, II or III full nature run:'); % truth trajectory
% Extract time step from truth file, e.g.
% L05M2_N960_K32_F8.00_I1_b1.00_c1.00_tf0.05_spinup100_tsteps10000_seed51422
% Extract tf, K, F, and seed from truth filename
tic;
string=strsplit(ftruth,'_');
idx=1;     model = str2double(string{idx}(end));
idx=idx+1; Nx = str2double(string{idx}(2:end));
idx=idx+1; Kparm = str2double(string{idx}(2:end));
idx=idx+1; F = str2double(string{idx}(2:end));
idx=idx+1; Iparm = str2double(string{idx}(2:end));
idx=idx+1; b = str2double(string{idx}(2:end));
idx=idx+1; c = str2double(string{idx}(2:end));
idx=idx+1; tF = str2double(string{idx}(3:end));
idx=idx+1; spinup =str2double(string{idx}(7:end));
idx=idx+1; Ncycles = str2double(string{idx}(7:end));
idx=idx+1; tseed = str2double(string{idx}(5:end-4));
[Xt,abstol,reltol] = load_truth(pathname,ftruth);
[N,Ncycles] = size(Xt);

%% Set up R and H for DA
name='Lorenz 96 Model I,II or III Nature Run Observations';
numlines=1;
prompt  = {'First var observed','Skip','Ob errors','Bias Factors',...
           'Reperr Kernel','Seeds','Plotit'};
default = {'1','1','[1.0]','[1.0]',...
           '[0 0 1 0 0]','[4207416]','1'};
answer=inputdlg(prompt,name,numlines,default);i=1;
% First observed variable
first=str2double(answer{i});i=i+1;
% Variables skipped (unobserved)
skip=str2double(answer{i});i=i+1;
mask=first:skip:N;
oblist=zeros(1,N);
oblist(mask)=1; % 1=observed, 0=unobserved
% Observation error variance
errlist=str2num(answer{i});i=i+1;
rlist=(errlist.').*ones(size(oblist));
% Obs bias factors
biasfaclist=str2num(answer{i});i=i+1;
% Error of representation spatial convolution
repr = str2num(answer{i});i=i+1;
repr_norm = repr./sum(repr);
% seed for random number generator
seedlist = str2num(answer{i});i=i+1;
% Plot option
plotit = str2double(answer{i});i=i+1;

%% Loop over seeds
for is=1:length(seedlist)
    seed = seedlist(is);
    % initialize random number generators
    rng(seed,"twister");
    % Loop over ob error sizes
    for ir=1:size(rlist,1)
        % Set up observation operator and ob error covariance matrix
        oberr = rlist(ir,:);
        [Nobs,H,R] = forward(oblist,oberr);
        % Create simulated obs with error covariance R and forward model H
        Rsqrt = chol(R, 'lower'); % R need not be diagonal
        yobs0 = H * Xt + Rsqrt * randn(Nobs,Ncycles);
        % Error of representation via spatial convolution
        yobs0 = conv2(repr_norm, 1, yobs0, 'same');
        for ibf=1:length(biasfaclist)
            bias_factor = biasfaclist(ibf);
            yobs = bias_factor.*yobs0;
            % Save obs
            fobs = [outfolder,...
                '\obs_tf',num2str(tF,'%4.2f\n'),...
                '_nc',num2str(Ncycles,'%d\n'),...
                '_R',num2str(R(1,1),'%4.2f\n'),...
                '_N',num2str(N,'%d\n'),...
                '_K',num2str(Kparm,'%d\n'),...
                '_F',num2str(F,'%5.2f\n'),...
                '_I',num2str(Iparm,'%d\n'),...
                '_b',num2str(b,'%5.2f\n'),...
                '_c',num2str(c,'%5.2f\n'),...
                '_tseed',num2str(tseed,'%d\n'),...
                '_oseed',num2str(seed,'%d\n'),...
                '_biasfac',num2str(bias_factor,'%5.3f\n'),...
                '_rep',num2str(repr,'%1d'),...
                '.mat'];
            save(fobs,'yobs','oblist','abstol','reltol');
            if plotit
                plot_obs(mask, Ncycles, Xt, yobs, N, bias_factor, repr);
            end
        end
    end
end
toc

function plot_obs(mask, Ncycles, Xt, yobs, N, bias_factor, repr)
    % plotting the states versus time should be your first check to see if the
    % result seems reasonable
    figure;
    subplot(211)
    node = 1;
    obloc = find(mask==node, 1);
    plot(1:Ncycles, Xt(node, :))
    grid on
    if ~isempty(obloc)
        hold on
        plot(1:Ncycles, yobs(node, :), 'r.')
    else
        text(18,18,sprintf('No obs at node %d',node));
    end
    ylabel('Time series, Location 1')
    xlabel('Time [s]')
    ylim([-10 20])
    ttl1=sprintf('Time series at node %d',node);
    repstr=strrep(num2str(repr,'%1d'),' ','');
    ttl2=sprintf('Bias factor %5.2f, Rep5 [%s]',bias_factor,repstr);
    title({ttl1,ttl2})
    subplot(212)
    ts = Ncycles;
    plot(1:N, Xt(:, ts))
    grid on
    hold on
    plot(mask, yobs(:, ts), 'ro')
    ylabel('Ring Map')
    xlabel('Location')
    ylim([-10 20])
    ttl1=sprintf('Ring map at time %d',ts);
    title({ttl1,ttl2})
end