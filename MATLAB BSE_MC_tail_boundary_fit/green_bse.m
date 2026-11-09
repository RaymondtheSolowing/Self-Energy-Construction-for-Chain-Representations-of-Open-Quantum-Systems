function [G,D,Dp] = green_bse( ...
    p,N,epsilon,g,omega,kappa,omegaBar,kappaBar)
%GREEN_BSE Qubit resolvent with a supplied homogeneous BSE continuation.
%
% Sites 1,...,N+1 retain their actual mapped-chain coefficients.  The
% continuation beyond N+1 is replaced by the homogeneous fixed-point BSE.

if N < 1 || N ~= floor(N)
    error('N must be a positive integer.');
end
if numel(omega) < N+1 || numel(kappa) < N
    error('Need at least N+1 frequencies and N chain couplings.');
end

omega = omega(:);
kappa = kappa(:);
links = [g; kappa(1:N)];

if nargout > 2
    [Sigma,dSigma] = sigma_uniform_tail(p,omegaBar,kappaBar);
else
    Sigma = sigma_uniform_tail(p,omegaBar,kappaBar);
end

for j = N+1:-1:1
    denominator = p + 1i*omega(j) + Sigma;
    if nargout > 2
        dSigma = -links(j)^2 .* (1+dSigma) ./ denominator.^2;
    end
    Sigma = links(j)^2 ./ denominator;
end

D = p + 1i*epsilon + Sigma;
G = 1./D;
if nargout > 2
    Dp = 1+dSigma;
end
end
