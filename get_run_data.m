%% Nature run, prior, posterior, parameters
function [Xt, Xb, Xa, trueparms, H, R, yobs] = get_run_data()
    % Nature run location
    outfolder = 'C:\Users\knisely\Desktop\Lorenz_05\'; % Local hard drive
    % Load truth, prior, derive posterior and parms names, load those
    [ftruth,tpath]=uigetfile([outfolder,'*L05*'],'Choose truth file:'); % truth trajectory
    [fobs,opath]=uigetfile([outfolder,'obs*'],'Choose truth file:'); % obs
    [fprior,ppath]=uigetfile([outfolder,'K*/prior*'],'Choose prior file:'); % prior

    truth = fullfile(tpath,ftruth);
    fprintf('Loading truth...\n');
    load(truth,'Xt');  % Ntfull x Nx
    Xt = Xt.'; % Nx x Ntfull

    obs = fullfile(opath,fobs);
    fprintf('Loading obs...\n');
    load(obs,'yobs');  % Nx x Ntfull

    prior = fullfile(ppath,fprior);
    fprintf('Loading prior...\n');
    load(prior,'ZKSQ');
    Xb = ZKSQ; % Nx x Nens x Nt
    [pri_loc, pri_ens, pri_time] = size(Xb);

    posterior = strrep(prior,'prior','posterior');
    fprintf('Loading posterior...\n');
    load(posterior,'XKSQ');
    Xa = XKSQ; % Nx x Nens x Nt
    [post_loc, post_ens, post_time] = size(Xa);

    if any([pri_loc, pri_ens, pri_time]) ~= any([post_loc, post_ens, post_time])
        error('Prior and posterior have incompatible sizes');
    end
    trueparmsname = strrep(prior,'prior','trueparms');
    load(trueparmsname,'trueparms','H','R');

    % Also need nature run spinup
    string = strsplit(ftruth,'_');
    spinup = str2double(string{9}(7:end));
    Xt = Xt(:,spinup+1:spinup+post_time); % Nx x Nt
    yobs = yobs(:,spinup+1:spinup+post_time); % Nx x Nt

    fprintf('Data read finished.\n')
end