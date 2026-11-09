function [x,w] = gauss_legendre_rule(n,a,b)
%GAUSS_LEGENDRE_RULE n-point Gauss-Legendre quadrature on [a,b].
% Self-contained Newton construction; no toolbox is required.

if nargin < 2
    a = -1;
    b = 1;
end
if n < 1 || n ~= floor(n) || ~(b > a)
    error('Need integer n >= 1 and b > a.');
end

m = ceil(n/2);
k = (1:m).';

% Good asymptotic initial guesses for the nonnegative roots.
z = cos(pi*(k-0.25)/(n+0.5));

for iter = 1:20
    [Pn,Pnm1] = legendre_pair(n,z);
    dPn = n*(z.*Pn-Pnm1)./(z.^2-1);
    dz = Pn./dPn;
    z = z-dz;
    if max(abs(dz)) < 5e-15
        break;
    end
end

[Pn,Pnm1] = legendre_pair(n,z); %#ok<ASGLU>
dPn = n*(z.*Pn-Pnm1)./(z.^2-1);
wz = 2./((1-z.^2).*dPn.^2);

if mod(n,2)==0
    x0 = [-z; flipud(z)];
    w0 = [wz; flipud(wz)];
else
    x0 = [-z(1:end-1); flipud(z)];
    w0 = [wz(1:end-1); flipud(wz)];
end

x = (a+b)/2 + (b-a)/2*x0;
w = (b-a)/2*w0;

[x,order] = sort(x);
w = w(order);
end

function [Pn,Pnm1] = legendre_pair(n,x)
if n==0
    Pn = ones(size(x));
    Pnm1 = zeros(size(x));
    return;
end

P0 = ones(size(x));
P1 = x;
if n==1
    Pn = P1;
    Pnm1 = P0;
    return;
end

for ell = 2:n
    P2 = ((2*ell-1).*x.*P1 - (ell-1).*P0)/ell;
    P0 = P1;
    P1 = P2;
end
Pn = P1;
Pnm1 = P0;
end
