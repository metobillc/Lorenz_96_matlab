function [yobs,varargout] = load_obs(pathname,fobs)
% Load arrays and optionally parameters from observations .mat file
% Bill Campbell
% Last modified 6/23/2022

MAXNARG = 3;
if nargout > MAXNARG
    error('load_obs has at most %d output arguments, not %d', ...
        MAXNARG, nargout)
end

fprintf('Loading obs from %s\n',[pathname,fobs]);
allvars = load([pathname,fobs]); % obs file
if isfield(allvars, 'yobs')
    yobs = allvars.yobs;
else
    error('yobs not found in obs file %s',fobs);
end

if nargout >= 2
    if isfield(allvars, 'abstol')
		varargout{1} = allvars.abstol;
	else
	    error('abstol not found in %s',fobs)
	end
end

if nargout == 3
    if isfield(allvars, 'reltol')
		varargout{2} = allvars.reltol;
	else
	    error('reltol not found in %s',fobs)
	end
end
end % function load_obs
