function [XZ,varargout] = load_ensemble(pathname,ftruth)
% Load arrays and parameters from ensemble output .mat file
% Bill Campbell
% Last modified 6/23/2022

MAXNARG = 5;
if nargout > MAXNARG
    error('load_ensemble has at most %d output arguments, not %d', ...
        MAXNARG, nargout)
end

fprintf('Loading ensemble from %s\n',[pathname,ftruth]);
allvars = load([pathname,ftruth]); % full ensemble
if isfield(allvars, 'ZKSQ') % Prior ensemble
    XZ = allvars.ZKSQ;
elseif isfield(allvars, 'XKSQ')  % Posterior ensemble
    XZ = allvars.XKSQ;
    if nargout >= 2
        if isfield(allvars, 'ensvarKSQ') % Ensemble variance
            varargout{1} = allvars.ensvarKSQ;
        else
            error('ensvarKSQ not found in ensemble file %s\n',ftruth);
        end
    elseif nargout == 3
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
