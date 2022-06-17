%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nature_run -- generate true trajectory
% for Lorenz 2005 models I, II, and III
% See also: draw_obs
% Bill Campbell
% Last modified 6/15/2022

start = datestr(clock);
fprintf('Started: %s\n',start);
outfolder = 'C:\Users\campbell\Documents\Lorenz_96_model'; % Local hard drive
prompt={'State variables','Kparm','Forcing','Iparm','Small damp','Coupling',...
    'Cycles','Spinup','Timestep','ODE AbsTol','ODE RelTol','Seedlist','Plotit'};
name='Lorenz 96 Model I, II and III Nature Run Parameters';
numlines=1;
default={'[960]','[32]','[15]','[12]','[10]','[2.5]',...
    '1000','100','[0.05]','[1e-6]','[1e-3]','[4921760]','[1]'};
answer=inputdlg(prompt,name,numlines,default);i=1;
% Number of state variables
Nxlist = str2num(answer{i});i=i+1;
% Number of K parameters (correlation length)
Kparmlist = str2num(answer{i});i=i+1;
% Forcing
Flist = str2num(answer{i});i=i+1;
% Number of I parameters (small scale correlation length)
Iparmlist = str2num(answer{i});i=i+1;
% Number of small scale damping parameters
blist = str2num(answer{i});i=i+1;
% Number of scale coupling parameters
clist = str2num(answer{i});i=i+1;
% Number of cycles to run
Ncycles = str2num(answer{i});i=i+1;
% Number of cycles to run
spinup = str2num(answer{i});i=i+1;
if isempty(spinup)
    spinup = 0;
end
% Time step (cycling interval)
tflist = str2num(answer{i});i=i+1;
% ODE Absolute Tolerance
abstol = str2num(answer{i});i=i+1;
if isempty(abstol)
    abstol = 1e-6;
end
% ODE Relative Tolerance
reltol = str2num(answer{i});i=i+1;
if isempty(reltol)
    reltol = 1e-3;
end
% seed for random number generator
seedlist = str2num(answer{i});i=i+1;
% Plot option
plotit = str2num(answer{i});i=i+1;
disp([prompt',answer]);
% Set ode solver options
% Save values in output file, but not in the filename
options=odeset('RelTol',reltol,'AbsTol',abstol); % defaults are RelTol=1e-3, AbsTol=1e-6
% Should abstol and reltol be saved in the filename, or in the file?
%% Loop over ring sizes
for in=1:length(Nxlist)
    Nx = Nxlist(in);
    % Loop over correlation lengths (K parameters)
    for ik=1:length(Kparmlist)
        Kparm = Kparmlist(ik);
        % create a structure to hold all the constants
        parms.K = Kparm;
        % Loop over forcing
        for jf=1:length(Flist)
            F = Flist(jf);
            parms.F = F;
            for ip=1:length(Iparmlist)
                Iparm = Iparmlist(ip);
                parms.I = Iparm;
                for ib=1:length(blist)
                    b = blist(ib);
                    for ic=1:length(clist)
                        c = clist(ic);
                        % Loop over time step sizes
                        for it=1:length(tflist)
                            tF = tflist(it);
                            tsteps = 0:tF:(Ncycles+spinup-1)*tF;
                            % Loop over seeds
                            for is=1:length(seedlist)
                                seed = seedlist(is);
                                % initialize random number generators
                                rng(seed,'twister');
                                % Perturb initial state
                                Xt0 = F*(ones(Nx,1)+0.01*randn(Nx,1));
                                % Integrate
                                if Kparm==1 % Model I run
                                    fprintf('Model I: Nx=%d, F=%5.2f, Kparm=%d, Iparm=%d\n',Nx,F,Kparm,Iparm);
                                    fprefix='\L05M1_N';
                                    model = 1; b=1; c=1;
                                elseif Iparm==1 % Model II run
                                    fprintf('Model II: Nx=%d, F=%5.2f, Kparm=%d, Iparm=%d\n',Nx,F,Kparm,Iparm);
                                    fprefix='\L05M2_N';
                                    model = 2; b=1; c=1;
                                else % Model III run
                                    fprintf('Model III: Nx=%d, F=%5.2f, Kparm=%d, Iparm=%d, ',Nx,F,Kparm,Iparm);
                                    fprintf('b=%5.2f, c=%5.2f\n',b,c);
                                    fprefix='\L05M3_N_';
                                    model = 3;
                                end
                                % create an anonymous function with the required inputs for ode45(), i.e.
                                % (t, x). Note that parms is set to the values above on creation of this
                                % function.
                                % The is needed for two reasons:
                                % 1. You can only store a m-file function in a variable by using anonymous
                                %    functions and ode45 requires that the function be anonymous
                                % 2. ode45 requires that the function only has t and x as arguments, but
                                %    we'd like to pass in the values of parms from the variables declared
                                %    above. This allows us to do that.
                                % Point 2 above is better than using global variables to share variables in
                                % all function scopes or declaring these parameters directly in the right
                                % hand side function where they have limited access.
                                parms.b = b;
                                parms.c = c;
                                loranon = @(t, x) circ_lorenz2005(t, x, parms);
                                % integrate the equations with one of the available integrators, in this
                                % case the Runga-Kutta 4,5 method (good for simple, non-stiff systems).
                                tic;
                                [~, Zt] = ode45(loranon, tsteps, Xt0, options);
                                toc
                                % Discard spinup
                                Zt = Zt((spinup+1):end,:);
                                % Decompose Zt(space,time) into large and small scales
                                if model==3
                                    Xt = zeros(size(Zt));
                                    alpha=(3*Iparm^2+3)/(2*Iparm^3+4*Iparm);
                                    beta=(2*Iparm^2+1)/(2*Iparm^2+Iparm^4);
                                    even=(mod(Kparm,2)==0);
                                    N = size(Zt,2);
                                    for j=-Iparm:Iparm
                                        fac = (alpha - beta*abs(j)) / (1.0 + even*(abs(j)==Iparm));
                                        for n=1:N
                                            Xt(:,n) = Xt(:,n) + fac * Zt(:,mod(n+j-1,N)+1);
                                        end
                                    end
                                    % Subtract Xt from Zt to get the small scales of Zt, stored in Yt
                                    % Lorenz 2005, eq. 13b
                                    Yt = Zt - Xt;
                                end
                                if plotit
                                    % plotting the states versus time should be your first check to see if the
                                    % result seems reasonable
                                    nature_plots(Zt,tsteps,tF,spinup)
                                end
                                % Save nature run
                                ftruth = [outfolder,...
                                    fprefix,num2str(Nx,'%d\n'),...
                                    '_K',num2str(Kparm,'%d\n'),...
                                    '_F',num2str(F,'%5.2f\n'),...
                                    '_I',num2str(Iparm,'%d\n'),...
                                    '_b',num2str(b,'%4.2f\n'),...
                                    '_c',num2str(c,'%4.2f\n'),...
                                    '_tf',num2str(tF,'%4.2f\n'),...
                                    '_spinup',num2str(spinup,'%d\n'),...
                                    '_tsteps',num2str(length(tsteps)-spinup,'%d\n'),...
                                    '_seed',num2str(seed,'%d\n'),'.mat'];
                                if model==3
                                    save(ftruth,'Zt','Xt','Yt','abstol','reltol');
                                else
                                    save(ftruth,'Zt','abstol','reltol');
                                end
                            end % Seedlist loop
                        end % Timestep loop
                    end
                end
            end
        end % Forcing loop
    end % Kparm loop
end % Nx loop
finish = datestr(clock);
fprintf('Started: %s, Finished %s\n',start,finish);
disp([prompt',answer]);
