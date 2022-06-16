function y = load_obs(pathname,fobs)
allvars = load([pathname,fobs]); % obs file
if isfield(allvars, 'y')
    y = allvars.y;
else
    error('y not found in obs file %s',fobs);
end
