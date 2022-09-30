function [Nobs,H,R,varargout] = forward(oblist,rlist)
% R is always diagonal, we can relax this assumption later.
% May not observe every gridpoint, so 
% H will be 1 where oblist is 1, 0 otherwise
Nobsall = length(oblist);
idx = find(oblist);
Nobs = length(idx);
H = zeros(Nobs,Nobsall);
for i=1:Nobs
    H(i,idx(i))=oblist(idx(i));
end
% Again since we may not observe every gridpoint,
% size of R will vary.
if length(rlist(idx))==1
    R = rlist*eye(Nobs,Nobs);
elseif length(rlist(idx))==Nobs
    R = diag(rlist(idx));
else  % default to 0.1
    R = 0.1*eye(Nobs,Nobs);
end
if nargout >= 4
    varargout{1} = idx;
end
