function detect_lorenz2005(Kparm,varargin)
% detect_lorenz2005.m detects models I, II, and III from the Lorenz 2005 paper
% returning the model number
% Bill Campbell
% Last modified 6/13/2022
%% Check arguments
Kparm=round(Kparm);
Iparm=1;
b=1;
c=1;
if nargin==4, % Model III
    if Kparm > 1, % Kparm==1 implies model I
        Iparm=round(varargin{1});
        if (Iparm<1 || Iparm>Kparm),
            error('Iparm=%d must be an integer between 1 and Kparm=%d, inclusive\n',Iparm,Kparm);
        end
        b=abs(varargin{2});
        c=abs(varargin{3});
        if Iparm==1,
            model=2;
        else
            model=3;
        end
    else
        model=1; % Kparm==1 implies model==1
    end
elseif nargin==1, % Model II
    if Kparm > 1,
        model=2;
    else
        model=1; % Kparm==1 implies model==1
    end
else
    error('detect_lorenz2005 takes either 1 (model I or II) or 4 (model III) arguments, not %d\n',nargin);
end
fprintf('L05 model %d: Kparm=%d, Iparm=%d, b=%f, c=%f\n',model,Kparm,Iparm,b,c);
