function [XZ,varargout] = load_ensemble(pathname,ftruth)
% Load arrays and parameters from ensemble output .mat file
% Bill Campbell
% Last modified 6/21/2022

if nargout > 2
    error('Max output arguments for load_ensemble is 2, not %d',nargout);
end

allvars = load([pathname,ftruth]); % full ensemble
if isfield(allvars, 'ZKSQ') % Prior ensemble
    XZ = allvars.ZKSQ;
    XZ = XZ.';
elseif isfield(allvars, 'XKSQ')  % Posterior ensemble
    XZ = allvars.XKSQ;
    XZ = XZ.';
    if nargout >= 1
        if isfield(allvars, 'ensvarKSQ') % Ensemble variance
            varargout{1} = allvars.ensvarKSQ;
        else
            error('ensvarKSQ not found in ensemble file %s\n',ftruth);
        end
    elseif nargout == 2
        if isfield(allvars, 'scoreKSQ') % Assimilation score
            varargout{1} = allvars.scoreKSQ;
        else
            error('scoreKSQ not found in ensemble file %s\n',ftruth);
        end
    end
else
    error('[Z,X]KSQ not found in ensemble file %s\n',ftruth);
end
end % function load_ensemble
