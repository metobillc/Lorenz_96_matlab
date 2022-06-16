function detect_lorenz2005(parms)
% detect_lorenz2005.m detects models I, II, and III from the Lorenz 2005 paper
% prints an informational message
% Bill Campbell
% Last modified 6/15/2022
%% Check arguments
Kparm=parms.K;
Iparm=parms.I;
b=parms.b;
c=parms.c;

if Kparm==1 % Model I run
    model = 1;
elseif Iparm==1 % Model II run
    model = 2;
else % Model III run
    model = 3;
end
fprintf('L05 model %d: Kparm=%d, Iparm=%d, b=%f, c=%f\n',model,Kparm,Iparm,b,c);
