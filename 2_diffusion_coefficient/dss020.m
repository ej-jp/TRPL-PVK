function ux = dss020 (xl,xu,n,u,v)
%
% Function dss020 computes the first derivative, "ux", of a
% variable "u" over the spatial domain xl .le. x .le. xu from
% the classical five point, fourth order finite difference approximations.
%
% Inputs:
% -------
% xl = lower boundary value of x
% xu = upper boundary value of x
% n = number of grid points in the x domain
% including the boundary points
% u = one dimensional array containing the value of u at
% the n grid points for which the derivative is to be computed
%
% Outputs:
% -------
% ux = one dimensional array containing the numerical values of
% the derivatives of u at the n grid points
%
% Matlab version written by: Antonio Flores T./ February 2008
%
%
% Spatial increment
%
dx = (xu-xl)/(n-1);
r4fdx = 1/(12*dx);
if v > 0
%
% Finite difference approximation for positive V
%
nm1 = n-1;
ux (1) = r4fdx*(-25*u(1)+48*u(2)-36*u(3)+16*u(4)-3*u(5));
ux (2) = r4fdx*(-3*u(1)-10*u(2)+18*u(3)-6*u(4)+u(5));
ux (3) = r4fdx*(u(1)-8*u(2)+8*u(4)-u(5));
for i = 4:nm1,
ux(i) = r4fdx*(-u(i-3)+6*u(i-2)-18*u(i-1)+10*u(i)+3*u(i+1));
end
ux (n) = r4fdx*(3*u(n-4)-16*u(n-3)+36*u(n-2)-48*u(n-1)+25*u(n));
%
% Finite difference approximation for negative V
%
else
nm3 = n-3;
ux (1) = r4fdx*(-25*u(1)+48*u(2)-36*u(3)+16*u(4)-3*u(5));
for i = 2:nm3,
ux(i) = r4fdx*(-3*u(i-1)-10*u(i)+18*u(i+1)-6*u(i+2)+u(i+3));
end
ux (n-2) = r4fdx*(u(n-4)-8*u(n-3)+8*u(n-1)-u(n));
ux (n-1) = r4fdx*(-u(n-4)+6*u(n-3)-18*u(n-2)+10*u(n-1)+3*u(n));
ux (n) = r4fdx*(3*u(n-4)-16*u(n-3)+36*u(n-2)-48*u(n-1)+25*u(n));
end
%-- End of the dss020.m file --  if true
% code
  end
