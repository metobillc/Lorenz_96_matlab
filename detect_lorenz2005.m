function [varargout] = detect_lorenz2005(parms)
% detect_lorenz2005.m detects models I, II, and III from the Lorenz 2005 paper
% prints an informational message
% Bill Campbell
% Last modified 6/23/2022
if nargout > 2
    error('detect_lorenz2005 has at most 2 output arguments, not %d',nargout)
end
%% Check arguments
Kparm=parms.K;
Iparm=parms.I;
b=parms.b;
c=parms.c;

if Kparm==1 % Model I run
    model = 1; b=1; c=1;
elseif Iparm==1 % Model II run
    model = 2; b=1; c=1;
else % Model III run
    model = 3;
end
fprintf('Model %d: Kparm=%d, Iparm=%d, b=%f, c=%f\n',model,Kparm,Iparm,b,c);
if nargout >= 1
    varargout{1} = model;
end
if nargout ==2
    prefix = sprintf('\\L05M%d_N',model);
    varargout{2} = prefix;
end
