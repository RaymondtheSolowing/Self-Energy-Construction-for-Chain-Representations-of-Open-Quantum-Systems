function [g,omega,kappa,info] = structured_chain_coeffs(Nmax,bath,nQuad)
%STRUCTURED_CHAIN_COEFFS Jacobi-chain coefficients from the continuum measure.
%
% The continuum measure is dmu = J(omega)/pi d omega.
% A Gauss-Legendre discretization of the measure is followed by a fully
% reorthogonalized Lanczos tridiagonalization of the multiplication operator.

if nargin < 3 || isempty(nQuad)
    nQuad = max(1600,30*Nmax);
end
if Nmax < 2 || Nmax ~= floor(Nmax)
    error('Nmax must be an integer >= 2.');
end

[x,w] = gauss_legendre_rule(nQuad,bath.omegaL,bath.omegaR);
J = structured_spectral_density(x,bath);
mu = w.*J/pi;

g2 = sum(mu);
if g2 <= 0
    error('The spectral measure has non-positive total weight.');
end
g = sqrt(g2);

Q = zeros(nQuad,Nmax);
Q(:,1) = sqrt(mu)/g;
omega = zeros(Nmax,1);
kappa = zeros(Nmax-1,1);

for n = 1:Nmax
    q = Q(:,n);
    z = x.*q;
    if n > 1
        z = z-kappa(n-1)*Q(:,n-1);
    end

    omega(n) = real(q'*z);
    z = z-omega(n)*q;

    % Full reorthogonalization keeps the first several dozen Jacobi
    % coefficients stable even for a strongly structured measure.
    z = z-Q(:,1:n)*(Q(:,1:n)'*z);

    if n < Nmax
        kappa(n) = norm(z);
        if kappa(n) <= 100*eps
            error('Lanczos chain terminated unexpectedly at site %d.',n);
        end
        Q(:,n+1) = z/kappa(n);
    end
end

info.quadratureNodes = nQuad;
info.totalMeasure = g2;
info.expectedTotalMeasure = bath.alpha;
info.couplingNormalizationError = abs(g2-bath.alpha);
end
