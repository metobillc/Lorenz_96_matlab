function [Zt,abstol,reltol,varargout] = load_truth(pathname,ftruth)
% Load arrays and parameters from truth run .mat file
% Bill Campbell
% Last modified 6/23/2022

MAXNARG = 5;
if nargout > MAXNARG
    error('load_truth has at most %d output arguments, not %d', ...
        MAXNARG, nargout)
end

fprintf('Loading truth from %s\n',[pathname,ftruth]);
allvars = load([pathname,ftruth]); % full nature run
if isfield(allvars, 'Zt')
    Zt = allvars.Zt;
    Zt = Zt.';
else
    error('Zt not found in nature run file %s',ftruth);
end

DEFAULT_ABSTOL = 1e-6;
DEFAULT_RELTOL = 1e-3;
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
if nargout >= 4
    if isfield(allvars, 'Xt')
	    Xt = allvars.Xt;
		varargout{1} = Xt.';
	else
	    error('Xt not found in %s',ftruth)
	end
end
if nargout == 5
    if isfield(allvars, 'yt')
	    Yt = allvars.Yt;
		varargout{2} = Yt.';
	else
	    error('Yt not found in %s',ftruth)
	end
end
end % function load_truth
