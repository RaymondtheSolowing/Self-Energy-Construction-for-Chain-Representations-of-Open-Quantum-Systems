function [J,dJ,normConst] = structured_spectral_density(omega,bath)
%STRUCTURED_SPECTRAL_DENSITY Finite-band Ohmic background with Gaussian feature.
%
% J(omega) = (2*pi*alpha/N) * omega *
%            [1 + A exp(-(omega-omega0)^2/(2*sigma^2))]
% on [omegaL,omegaR], and zero outside.
%
% N is chosen so that integral J(omega)/pi d omega = alpha.

required = {'omegaL','omegaR','alpha','bumpAmplitude','bumpCenter','bumpWidth'};
for k = 1:numel(required)
    if ~isfield(bath,required{k})
        error('bath.%s is required.',required{k});
    end
end

L = bath.omegaL;
R = bath.omegaR;
alpha = bath.alpha;
A = bath.bumpAmplitude;
omega0 = bath.bumpCenter;
sigma = bath.bumpWidth;

if ~(R > L) || alpha <= 0 || sigma <= 0
    error('Invalid structured-bath parameters.');
end

% Integral of omega over [L,R].
Ibase = 0.5*(R^2-L^2);

% Integral of omega*Gaussian over [L,R].
F = @(x) -sigma^2 .* exp(-(x-omega0).^2/(2*sigma^2)) + ...
    omega0*sigma*sqrt(pi/2) .* erf((x-omega0)/(sqrt(2)*sigma));
Igauss = F(R)-F(L);

normConst = 2*(Ibase + A*Igauss);
prefactor = 2*pi*alpha/normConst;

omegaArray = omega;
e = exp(-(omegaArray-omega0).^2/(2*sigma^2));
J = prefactor .* omegaArray .* (1 + A*e);

dJ = prefactor .* ( ...
    1 + A*e - omegaArray .* A .* e .* (omegaArray-omega0)/sigma^2);

outside = (omegaArray < L) | (omegaArray > R);
J(outside) = 0;
dJ(outside) = 0;
end
