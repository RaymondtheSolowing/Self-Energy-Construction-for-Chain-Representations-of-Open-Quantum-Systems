function [A0,info] = finite_chain_amplitude_spectral( ...
    N,t,epsilon,g,omega,kappa,gamma)
%FINITE_CHAIN_AMPLITUDE_SPECTRAL Survival amplitude of a finite chain.
%
% Basis [q,a1,...,aN]. gamma_n enters as -i*gamma_n/2.

if nargin < 7 || isempty(gamma)
    gamma = zeros(N,1);
end

t = t(:).';
gamma = gamma(:);

if numel(omega) < N
    error('omega must contain at least N entries.');
end
if N > 1 && numel(kappa) < N-1
    error('kappa must contain at least N-1 entries.');
end
if numel(gamma) ~= N
    error('gamma must contain N entries.');
end

diagH = [epsilon; omega(1:N)-1i*gamma/2];
H = diag(diagH);
H(1,2) = g;
H(2,1) = g;

for n = 1:N-1
    H(n+1,n+2) = kappa(n);
    H(n+2,n+1) = kappa(n);
end

isHermitian = max(abs(gamma)) == 0;
if isHermitian
    [V,D] = eig(real(H));
    E = diag(D);
    coeff = abs(V(1,:)).^2;
    coeff = coeff(:);
else
    [V,D] = eig(H);
    E = diag(D);
    e0 = zeros(N+1,1);
    e0(1) = 1;
    c0 = V\e0;
    coeff = V(1,:).' .* c0;
end

A0 = spectral_sum(t,E,coeff);

info.eigenvalues = E;
info.isHermitian = isHermitian;
end

function A = spectral_sum(t,E,coeff)
A = zeros(size(t));
if numel(t)==1
    A(1) = coeff.'*exp(-1i*E*t(1));
    return;
end

if abs(t(1)) < 100*eps
    dt = t(2)-t(1);
    if max(abs(diff(t)-dt)) < 1e-10*max(1,abs(dt))
        phase = ones(size(E));
        step = exp(-1i*E*dt);
        for m = 1:numel(t)
            A(m) = coeff.'*phase;
            phase = phase.*step;
        end
        return;
    end
end

chunkSize = 500;
for i1 = 1:chunkSize:numel(t)
    i2 = min(numel(t),i1+chunkSize-1);
    A(i1:i2) = coeff.'*exp(-1i*(E*t(i1:i2)));
end
end
