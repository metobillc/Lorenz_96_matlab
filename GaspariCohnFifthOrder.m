function C0 = GaspariCohnFifthOrder(z,c)
% The value of the fifth order piecewise rational function C0(z,1/2,c)
% in equation 4.10 in Gaspari & Cohn's 1999 QJRMS 166, #125, p. 748
% Bill Campbell
% Last modified 9/9/2016
% Now works for vector arguments -- can extend to matrix args if needed
if c <= 0.0,
  fprintf('The second argument (c) must be positive');
  C0 = NaN;
  return
end
z = abs(z);
x = z./c;
x1 = x(z<=c);
C1 = 1+x1.*(0+x1.*(-5/3+x1.*(5/8+x1.*(1/2+x1.*(-1/4)))));
x2 = x(z>c & z<=2*c);
C2 = -2./(3*x2)+4+x2.*(-5+x2.*(5/3+x2.*(5/8+x2.*(-1/2+x2.*(1/12)))));
x3 = x(z>2*c);
C3 = zeros(size(x3));
if isscalar(z) || isrow(z),
    C0 = [C1 C2 C3];
elseif iscolumn(z),
    C0 = [C1;C2;C3];
end