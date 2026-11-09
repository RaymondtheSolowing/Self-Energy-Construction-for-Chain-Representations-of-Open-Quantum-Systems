function Sigma = continuum_self_energy_regularized(p,bath,nQuad)
%CONTINUUM_SELF_ENERGY_REGULARIZED Accurate Cauchy transform for the
% finite-band structured Ohmic/Gaussian bath in structured_spectral_density.m.
%
% Sigma(p) = (1/pi) integral_L^R J(nu)/(p+i*nu) dnu.
% For p = eta-i*omega, eta>0, let z=i*p=omega+i*eta, and analytically
% continue the formula for J(nu)/pi from the open spectral interval to z:
%
% Sigma(p) = -i * { integral_L^R [f(nu)-f(z)]/(nu-z) dnu
%                     + f(z) [log(R-z)-log(L-z)] },  f=J/pi.
%
% The pole-subtracted integrand is smooth even at small eta. A direct
% nQuad-node quadrature of J/(p+i*nu) can introduce an apparent error floor
% in the frequency-domain BSE closure-depth test when eta is small.
%
% This specialized implementation assumes the analytic Ohmic/Gaussian
% interior form used in structured_spectral_density; it is NOT a generic
% routine for an arbitrary tabulated spectral density.

if nargin < 3 || isempty(nQuad)
    nQuad = 1200;
end
if ~isnumeric(p) || any(~isfinite(p(:))) || any(real(p(:)) <= 0)
    error('p must be finite with Re(p)>0.');
end

L = bath.omegaL;
R = bath.omegaR;
z = 1i*p(:).';
if any(real(z) < L-1e-12) || any(real(z) > R+1e-12)
    error('Regularized routine is intended for frequencies in [omegaL,omegaR].');
end
if any(imag(z) > bath.bumpWidth/2)
    error('eta too large for stable analytic-Gaussian subtraction; use ordinary quadrature.');
end

[x,w] = gauss_legendre_rule(nQuad,L,R);
[~,~,normConst] = structured_spectral_density((L+R)/2,bath);
f = @(v) (2*bath.alpha/normConst).*v.* ...
    (1+bath.bumpAmplitude*exp(-(v-bath.bumpCenter).^2/ ...
    (2*bath.bumpWidth^2)));

fx = f(x);
fz = f(z);
result = zeros(size(z));
chunk = 200;
for i1 = 1:chunk:numel(z)
    ii = i1:min(i1+chunk-1,numel(z));
    zi = z(ii);
    smooth = (fx - fz(ii))./(x - zi);
    integralRegular = sum(w.*smooth,1);
    logTerm = log(R-zi)-log(L-zi);
    result(ii) = -1i*(integralRegular + fz(ii).*logTerm);
end
Sigma = reshape(result,size(p));
end
