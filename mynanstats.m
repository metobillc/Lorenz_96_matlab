function [avg,stdev,varargout] = mynanstats(x,flag)
%MYNANSTATS Average and standard deviation ignoring NaNs.
%   MYNANSTATS(X) returns the same standard deviation treating NaNs 
%   as missing values.
%
%   For column vectors, MYNANSTATS(X) is the average and standard deviation
%   of the non-NaN elements in X.  For matrices (and row vectors),
%   MYNANSTATS(X) returns row vectors containing the average and standard
%   deviation of each column, ignoring NaNs.
%
%   MYNANSTATS(X) normalizes the standard deviation by (N-1) where N
%   is the sequence length. This makes stdev.^2 the best unbiased
%   estimate of the variance if X is a sample from a normal distribution.
%
%   MYNANSTATS(X,1) normalizes by N and produces the second moment of
%   the sample about its mean.  MYNANSTATS(X,0) is the same as 
%   MYNANSTATS(X).
%
%   Optional output argument computest the RMS error, which is the square
%   root of the sum of the squares of the avg and stdev.
%
%   See also MYNANSTD, MYNANMEAN, MYNANMEDIAN, MYNANMIN, MYNANMAX, MYNANSUM.
%
% Bill Campbell
% Last modified 5/14/2014

if nargin<2, flag = 0; end

if isempty(x) % Check for empty input.
    avg = NaN;
    stdev = NaN;
    if (nargout==3)
        varargout{1}=NaN;
    end
    return
end

[rows,cols] = size(x);

if rows==1	% Only take mean if a COLUMN vector
  avg = x;
  stdev = 0;
  return
end

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);

% Find mean
avg = nanmean(x);
%avg = mynanmean(x);

count = rows-sum(nans);
if cols==1
   x = x - avg;
else
   x = x - avg(ones(rows,1),:);
end

% Replace NaNs with zeros.
x(i) = zeros(size(i));

% Protect against a column of all NaNs
i = find(count==0);
count(i) = ones(size(i));
if flag
   stdev = sqrt(sum(x.*x)./max(count,1));
else
   stdev = sqrt(sum(x.*x)./max(count-1,1));
end
stdev(i) = i + NaN;
if nargout==3
    varargout{1} = sqrt(avg.^2+stdev.^2);
end
