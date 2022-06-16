%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw_obs from true trajectory
% Bill Campbell
% Last modified 6/15/2022

%% Get nature run
outfolder = 'C:\Users\campbell\Documents\MATLAB\Lorenz_96_model'; % Local hard drive
[ftruth,pathname]=uigetfile(outfolder,'Load Lorenz 96 Model I, II or III full nature run:'); % truth trajectory
% Extract time step from truth file, e.g.
% L05M2_N960_K32_F8.00_I1_b1.00_c1.00_tf0.05_spinup100_tsteps10000_seed51422
% Extract tf, K, F, and seed from truth filename
string=strsplit(ftruth,'_');
model=str2double(string{1}(end))
Kparm = str2double(string{4}(2:end))
F = str2double(string{5}(2:end))
Iparm = str2double(string{6}(2:end))
b = str2double(string{7}(2:end))
c = str2double(string{8}(2:end))
tF = str2double(string{9}(3:end))
tseed = str2double(string{12}(5:end-4))
tic;
[Xt,abstol,reltol] = load_truth(pathname,ftruth);
[N,Ncycles] = size(Xt);

%% Set up R and H for DA
prompt={'First var observed','Skip','Ob error','Seedlist','Plotit'};
name='Lorenz 96 Model I,II or III Nature Run Observations';
numlines=1;
default={'1','1','1.0','4207416','1'};
answer=inputdlg(prompt,name,numlines,default);i=1;
% First observed variable
first=str2double(answer{i});i=i+1;
% Variables skipped (unobserved)
skip=str2double(answer{i});i=i+1;
mask=first:skip:N;
oblist=zeros(1,N);
oblist(mask)=1; % 1=observed, 0=unobserved
% Observation error variance
err=str2double(answer{i});i=i+1;
rlist=err.*ones(size(oblist));
% seed for random number generator
seedlist = str2double(answer{i});i=i+1;
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
        [Nobs,H,R] = forward(oblist,rlist(ir,:));
        % Create simulated obs with error covariance R and forward model H
        yobs = H*Xt + sqrt( R )*randn(Nobs,Ncycles);
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
            '_oseed',num2str(seed,'%d\n'),'.mat'];
        save(fobs,'yobs','oblist','abstol','reltol');
    end
end
toc
if plotit
    % plotting the states versus time should be your first check to see if the
    % result seems reasonable
    figure;
    subplot(211)
    node = 1;
    obloc = find(mask==node);
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
    title(sprintf('Time series at node %d',node))
    subplot(212)
    ts = 1000;
    plot(1:N, Xt(:, ts))
    grid on
    hold on
    plot(mask, yobs(:, ts), 'b.')
    ylabel('Ring Map')
    xlabel('Location')
    ylim([-10 20])
    title(sprintf('Ring map at time %d',ts))
end
