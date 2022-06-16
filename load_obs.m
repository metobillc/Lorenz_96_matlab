function [y,varargout] = load_obs(pathname,fobs)
% Load arrays and optionally parameters from observations .mat file
% Bill Campbell
% Last modified 6/16/2022

allvars = load([pathname,fobs]); % obs file
if isfield(allvars, 'y')
    y = allvars.y;
else
    error('y not found in obs file %s',fobs);
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
