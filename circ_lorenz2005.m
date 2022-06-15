function dZ = circ_lorenz2005(~,Z,parms)
% circ_lorenz2005.m implements models I, II, and III from the Lorenz 2005 paper
% Bill Campbell
% Last modified 6/15/2022
%% Check arguments
N=size(Z,1);
K=round(parms.K);
if (K >= N/2)
    error('K=%d should be less than N/2=%d',K,N);
end
F = abs(parms.F); % Forcing must be positive
if isscalar(F)
    F = F.*ones(N,1);
end
Iparm=round(parms.I);
if (Iparm<1 || Iparm>K)
    error('Iparm=%d must be an integer between 1 and K=%d, inclusive\n',Iparm,K);
end
if K ==1
  %% Simpler code for model I
  dZ = - Z + F; % Lorenz 2005, model I damping and forcing terms
  dZ = dZ + (circshift(Z,-1) - circshift(Z,2)) .* circshift(Z,1);
  return
elseif Iparm==1
  model=2;
else
  model=3;
end
%% Initialize variables and compute parameters
% K even vs. odd
even=(mod(K,2)==0);
if model==2
    X = Z;
    Y = zeros(N,1);
else % model III
%% Model III: Compute the large scales of Z, and store them in X
% Lorenz 2005, eq. 13a
    X = zeros(N,1);
    alpha=(3*Iparm^2+3)/(2*Iparm^3+4*Iparm);
    beta=(2*Iparm^2+1)/(2*Iparm^2+Iparm^4);
    for j=-Iparm:Iparm
        fac = (alpha - beta*abs(j)) / (1.0 + even*(abs(j)==Iparm));
        X = X + fac * circshift(Z,-j);
    end
    % Subtract X from Z to get the small scales of Z, stored in Y
    % Lorenz 2005, eq. 13b
    Y = Z - X;
end
%% Now compute the time derivative of Zn
% Lorenz 2005, eq. 15
% Must calculate all of the weights first
W=zeros(N,1);
J = (K-1+even)/2;
for i=-J:J % Lorenz 2005, eq 9
    fac = 1.0/K/(1.0 + even*(abs(i)==J));
    W = W + fac * circshift(X,i);
end
b = parms.b;
c = parms.c;
%% Damping and forcing terms (pieces of Lorenz 2005, eq 15)
dZ = - X - b*Y + F; % Y=0 and X=Z for model II, so no need for if block
% Nonlinear advection terms for large scale
dZ = dZ - circshift(W,2*K) .* circshift(W,K);
for j = -J:J % Lorenz 2005, eq 10 summation
    fac = 1.0/K/(1.0 + even*(abs(j)==J));
    dZ = dZ + fac * circshift(W,K-j) .* circshift(X,-K-j);
end
%% Additional terms for model III
if model==3
    % Nonlinear advection terms for small scale
    dZ = dZ + b^2 * circshift(Y,1) .* (circshift(Y,-1) - circshift(Y,2));
    % Scale coupling terms
    dZ = dZ + c * (circshift(Y,1) .* circshift(X,-1) - circshift(Y,2) .* circshift(X,1));
end
