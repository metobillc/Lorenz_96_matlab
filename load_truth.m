function [Zt,abstol,reltol] = load_truth(pathname,ftruth)
DEFAULT_ABSTOL = 1e-6
DEFAULT_RELTOL = 1e-3

allvars = load([pathname,ftruth]); % full nature run
if isfield(allvars, 'Zt')
    Zt = allvars.Zt;
    Zt = Zt.';
else
    error('Zt not found in nature run file %s',ftruth);
end
if isfield(allvars, 'abstol')
    abstol = allvars.abstol;
else
    abstol = DEFAULT_ABSTOL;
end
if isfield(allvars, 'reltol')
    reltol = allvars.reltol;
else
    reltol = DEFAULT_RELTOL;
end
