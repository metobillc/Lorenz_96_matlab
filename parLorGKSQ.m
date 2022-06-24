%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parLorGKSQ with prior inflation
% Using Lorenz 96 Model I, II or III
% Global Kalman Square Root filter
% Bill Campbell
% Last modified 6/15/2022

start = datestr(clock);
fprintf('Started: %s\n',start);
plotit = 0;
%% Parameter input
prompt={'Solver','Cycles','First var observed','Grid Skip','Time Skip',...
    'Ob error','Ensemble size(s)','Confidence level','Prior inflation(s)',...
    'Cycling interval','Seedlist','Cycles per print','Localization Type',...
    'Localization radius','Save state'};
name='Filter Input';
numlines=1;
altdefault={'@ode45','200','1','1','1',...
   '1.0','[500]','[0.95]','[1.01]',...
   '0.05','29418','20','none',...
   '[10]','1'};
% altdefault={'@ode45','500','1','1','1',...
%     '1','100','[0.95]','[1.02]',...
%     '0.05','8142','10','none',...
%     '[10]','1'};
answer=inputdlg(prompt,name,numlines,altdefault);i=1;
%answer=inputdlg(prompt,name,numlines,default);i=1;
% Built in ODE solver to use
solver = eval(answer{i});i=i+1;
% Number of cycles to run
Ncycles = str2num(answer{i});i=i+1;
% First observed variable
first = str2num(answer{i});i=i+1;
% Variables skipped (unobserved)
skip = str2num(answer{i});i=i+1;
% Time steps skipped (completely unobserved)
tskip = str2num(answer{i});i=i+1;
% observation error variance
oberr = str2num(answer{i});i=i+1;
% number of ensemble members --- 32, 64, 128, 256, 512, 1024
Kvec = str2num(answer{i});i=i+1;
% Confidence level (0.95 is typical)
ci = str2num(answer{i});i=i+1;
% inflation parameters (tune these for each experiment)
alphavec = str2num(answer{i});i=i+1;
% cycling interval -- should get this from nature run
tFvec = str2num(answer{i});i=i+1;
fstr = 'Global Kalman Square Root';
seedlist = str2num(answer{i});i=i+1;
printcycle = str2num(answer{i});i=i+1;
% Localization type
switch lower(answer{i})
    case {'gaspari','gasparicohn','g','gc','gaspari-cohn','g-c'}
        loctype = 1;
        locstr = 'g';
    case {'boxcar','box','flat','b'}
        loctype = 2;
        locstr = 'b';
    case {'triangular','tri','t','tent'}
        loctype = 3;
        locstr = 't';
    otherwise % No localization
        loctype = 0;
        locstr = 'n';
end
i=i+1;
% Localization radius (in gridpoints)
locradlist = abs(round(str2num(answer{i})));i=i+1;
savestate = logical(str2num(answer{i}));i=i+1;
disp([prompt',answer]);

%% Nature run and observations of the truth
% Truth file e.g. L05M3_N960_K32_F15.00_I12_b10.00_c2.50_tf0.05_spinup100_tsteps10000_seed51422
% Obs file e.g. obs_tf0.05_nc10000_R1.00_N960_K32_F15.00_I12_b10.00_c2.50_tsteps960_tseed51422_oseed73033
% Nature run location
outfolder = 'C:\Users\campbell\Documents\MATLAB\Lorenz_96_model\'; % Local hard drive
% Load truth and compatible observations
compatible = false;
while ~compatible
    [ftruth,truepath]=uigetfile([outfolder,'*L05*'],'Choose truth file:'); % truth trajectory
    [fobs,obspath]=uigetfile([outfolder,'obs*'],'Choose obs file:'); % obs trajectory
    % Check compatibility of truth and obs
    compatible = new_check_obs_compatibility(ftruth,fobs);
end
% Extract Lorenz 96 Model I, II or III parameters from nature run filename
string = strsplit(ftruth,'_');
parms.K = str2double(string{3}(2:end));
parms.F = str2double(string{4}(2:end));
parms.I = str2double(string{5}(2:end));
parms.b = str2double(string{6}(2:end));
parms.c = str2double(string{7}(2:end));

fname = cell(length(alphavec),1);
fid = zeros(length(alphavec),1);
tot_err = zeros(length(alphavec),1);
tot_var = tot_err;
NPROC = feature('numcores'); % This machine has 4 processors
reserve = 1;
allprocs = NPROC - reserve; % Use all processes available
% If have enough processors, parallelize over both the alpha parameter and
% number of ensemble members
parprocs = min([allprocs,length(alphavec)]); % Devote max resources to parameters, devote remainder to ensemble members
% To run with only one processor, uncomment the two lines below
% parprocs = 0;
% myprocs = 1;
myprocs = parprocs; % Can be less then allprocs for small ensembles
if myprocs > 1 && isempty(gcp)
    parpool(myprocs);
end
%profile on
for itf = 1:length(tFvec) % Forward model time steps
    tF = tFvec(itf);
    % Load nature run and observations for this time step
    ftruth = regexprep(ftruth,'tf_....',['tf_',num2str(tF,'%4.2f\n')]);
    fobs = regexprep(fobs,'tf_....',['tf_',num2str(tF,'%4.2f\n')]);
    [Xt,abstol,reltol] = load_truth(truepath,ftruth); % full nature run
    parms.abstol = abstol;
    parms.reltol = reltol;
    y = load_obs(obspath,fobs);  % observations consistent with nature run
    % Generate initial ensemble from nature run climatology
    % Equally spaced through time
    ics = floor(linspace(1,size(Xt,2),max(Kvec)+1)); % Kmax+1 x 1
    XIC = Xt(:,ics(2:end)); % Nx x Kmax
    %keyboard;
    if (Ncycles > size(y,2))
        fprintf('Error -- obs file not long enough: Ncycles=%d\n > len(obs)',...
            Ncycles,size(y,2));
        break;
    else % Only use needed portion of nature run and obs
        Xt = Xt(:,1:Ncycles);  % Nx x Ncycles
        y = y(:,1:Ncycles);  % Nx x Ncycles
    end
    % Set oblist and rlist here
    Nx = size(y,1);
    oblist = zeros(1,Nx);
    rlist = oblist;
    % Observed gridpoints
    gmask = first:skip:Nx;
    oblist(gmask)=1;
    rlist(gmask)=oberr;
    % Set up observation operator and ob error covariance matrix
    % H will be Nobs x Nx, R will be Nobs x Nobs, Nobs <= Nx
    [~,H,R] = forward(oblist,rlist);
    % Observed time steps -- implement later
    %tmask = first:tskip:Ncycles; this does not work
    % Compress the observation vector to match compressed R
    y = y(gmask,:);  % Nobs x Ncycles
    for iK = 1:length(Kvec) % Ensemble sizes
        K = Kvec(iK);
        if ~exist([outfolder,'K',num2str(K,'%d\n')], 'dir')
            mkdir([outfolder,'K',num2str(K,'%d\n')]);
        end
        for iseed = 1:length(seedlist)
            myseed = seedlist(iseed);
            for irad = 1:length(locradlist) % Localization radii
                locrad = locradlist(irad);
                for ialpha = 1:length(alphavec) % Must set filenames outside parfor loop if myprocs < length(alphavec)
                    alpha = alphavec(ialpha);
                    fname{ialpha,:} = [outfolder,'K',num2str(K,'%d\n'),...
                        '\alpha_',num2str(alpha,'%4.2f\n'),...
                        '_tF_',num2str(tF,'%4.2f\n'),...
                        '_loc_',locstr,num2str(locrad,'%d\n'),...
                        '_nc_',num2str(Ncycles,'%d\n'),...
                        '_1st_',num2str(first,'%d\n'),...
                        '_skip_',num2str(skip,'%d\n'),...
                        '_',num2str(myseed,'%d\n'),'.txt'];
                end
                if parprocs > 1
                    parfor (ialpha = 1:length(alphavec),parprocs) % Prior inflation factors
                        alpha = alphavec(ialpha);
                        fid(ialpha) = fopen(fname{ialpha,:},'at');
                        fprintf(fid(ialpha),'K=%d, tF=%4.2f, loc %s, rad %d, %d cycles\t%s\t%s\t%s\t%s\t%s\n',...
                            K,tF,locstr,locrad,Ncycles,fstr); %#ok<*PFBNS>
                        tstart=tic;
                        % Execute a single DA cycle by callint parDA_GKSQ
                        [tot_err(ialpha,:),tot_avar(ialpha,:)] = parDA_GKSQ(solver,XIC,Xt,y,outfolder,H,R,...
                            first,skip,tF,ci,alpha,K,myseed,fstr,savestate,...
                            printcycle,loctype,locstr,locrad,parms);
                        tend=toc(tstart);
                        thour=floor(tend/3600);
                        tminute=floor((tend-thour*3600)/60);
                        tsecond=round(mod(tend,60));
                        fprintf('Elapsed time %0.2d:%0.2d:%0.2d\n',thour,tminute,tsecond);
                        fprintf('%4.2f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                            alpha,tot_err(ialpha,:));
                        fprintf(fid(ialpha),'%4.2f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                            alpha,tot_err(ialpha,:));
                        fclose(fid(ialpha));
                    end % alpha
                else
                    for ialpha = 1:length(alphavec) % Prior inflation factors
                        alpha = alphavec(ialpha);
                        fid(ialpha) = fopen(fname{ialpha,:},'at');
                        fprintf(fid(ialpha),'K=%d, tF=%4.2f, loc %s, rad %d, %d cycles\t%s\t%s\t%s\t%s\t%s\n',...
                            K,tF,locstr,locrad,Ncycles,fstr);
                        tstart=tic;
                        % Execute a single DA cycle by callint parDA_GKSQ
                        [tot_err(ialpha,:),tot_avar(ialpha,:)] = parDA_GKSQ(solver,XIC,Xt,y,outfolder,H,R,...
                            first,skip,tF,ci,alpha,K,myseed,fstr,savestate,...
                            printcycle,loctype,locstr,locrad,parms);
                        tend=toc(tstart);
                        thour=floor(tend/3600);
                        tminute=floor((tend-thour*3600)/60);
                        tsecond=round(mod(tend,60));
                        fprintf('Elapsed time %0.2d:%0.2d:%0.2d\n',thour,tminute,tsecond);
                        fprintf('%4.2f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                            alpha,tot_err(ialpha,:));
                        fprintf(fid(ialpha),'%4.2f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                            alpha,tot_err(ialpha,:));
                        fclose(fid(ialpha));
                    end % alpha
                end % myprocs > 0
                % Plot time series of truth and ensemble mean
                if plotit
                    myfilt = 'GKSQ';
                    posterior = [outfolder,'K',num2str(K,'%d\n'),'\',...
                        'posterior_alpha_',num2str(alpha,'%4.2f\n'),...
                        '_tf_',num2str(tF,'%4.2f\n'),...
                        '_R_',num2str(R(2,2),'%4.2f\n'),...
                        '_loc_',locstr,num2str(locrad,'%d\n'),...
                        '_nc_',num2str(Ncycles,'%d\n'),...
                        '_1st_',num2str(first,'%d\n'),...
                        '_skip_',num2str(skip,'%d\n'),...
                        '_seed_',num2str(myseed,'%d\n'),...
                        '_',myfilt,'.mat'];
                    Xa = load(posterior);
                    fn = fieldnames(Xa);
                    Xa = Xa.(fn{1});
                    figure;
                    plot(Xt(:,1),'k-');
                    hold on; grid on;
                    plot(mean(Xa(:,:,1),2),'r-');
                    legend('Truth','Ensmean');
                    title(sprintf('%s time=1',myfilt));
                    for t=2:Ncycles
                        pause(0.5);
                        clf;
                        plot(Xt(:,t),'k-');
                        hold on; grid on;
                        plot(mean(Xa(:,:,t),2),'r-');
                        legend('Truth','Ensmean');
                        title(sprintf('%s time=%d',myfilt,t));
                    end % for t=2:Ncycles
                    pause;
                  
                end % if plotit
            end % localization radii
        end % seeds
    end % ensemble sizes
end % time steps
% if ~isempty(gcp('nocreate')),
%     delete(gcp);
% end
finish = datestr(clock);
fprintf('Started: %s, Finished %s\n',start,finish);
disp([prompt',answer]);
fclose('all');
%profile viewer
