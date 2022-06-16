function [Xt,abstol,reltol] = load_truth(pathname,ftruth)
allvars = load([pathname,ftruth]); % full nature run
if isfield(allvars, 'Zt') % Zt is correct here, NOT Xt
    Xt = allvars.Zt;  % See prior comment
    Xt = Xt';
else
    error('Zt not found in nature run file %s',ftruth);
end
if isfield(allvars, 'abstol')
    abstol = allvars.abstol;
else
    abstol = 1e-6;
end
if isfield(allvars, 'reltol')
    reltol = allvars.reltol;
else
    reltol = 1e-3;
end
