function CL = covloc(N,loctype,locrad,varargin)
% NxN Covariance localization matrix
N = round(N);
CL = eye(N);
if locrad < 2,
    return
end
% Periodic boundary conditions
if nargin==3,
    periodic = false;
else
    periodic = logical(varargin{1});
end
% loctypes are GC, boxcar, tent, none
switch loctype
    case 1 % Gaspari-Cohn
        cutoff = floor(locrad/2);
        z = 0:locrad-1;
        C0 = GaspariCohnFifthOrder(z,cutoff);
        % Now put C0 values on the off diagonals
        CL = eye(N);
        if periodic, % must do this first if locrad > N/2
            parfor i=1:locrad-1,
                dvec = ones(i,1);
                CL = CL + C0(i+1).*(diag(dvec,N-i)+diag(dvec,-(N-i)));
            end
        end
        parfor i=1:locrad-1,
            dvec = ones(N-i,1);
            CL = CL + C0(i+1).*(diag(dvec,i)+diag(dvec,-i));
        end
    case 2 % Boxcar
        CL = eye(N);
        if periodic,
            parfor i=1:locrad-1,
                dvec = ones(i,1);
                CL = CL + (diag(dvec,N-i)+diag(dvec,-(N-i)));
            end
        end
        parfor i=1:locrad-1,
            dvec = ones(N-i,1);
            CL = CL + diag(dvec,i)+diag(dvec,-i);
        end
    case 3 % Tent
        CL = eye(N);
        C0 = linspace(1,0,locrad+1);
        if periodic, % must do this first if locrad > N/2
            parfor i=1:locrad-1,
                dvec = ones(i,1);
                CL = CL + C0(i+1).*(diag(dvec,N-i)+diag(dvec,-(N-i)));
            end
        end
        parfor i=1:locrad-1,
            dvec = ones(N-i,1);
            CL = CL + C0(i+1).*(diag(dvec,i)+diag(dvec,-i));
        end
    otherwise % None
        CL = ones(N);
end
% Debug
% figure;
% imagesc(CL);
% axis square;
% colormap linspecer;
% colorbar('v');
% title(sprintf('Covariance Localization: loctype=%d, locrad=%d\n',loctype,locrad));
end