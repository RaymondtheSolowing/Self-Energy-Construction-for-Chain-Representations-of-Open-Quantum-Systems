function SigmaPort = boundary_self_energy_chain( ...
    p,N,g,omega,kappa,gamma,SigmaAfterN)
%BOUNDARY_SELF_ENERGY_CHAIN Inward propagation of a boundary response.
%
% SigmaAfterN is the self-energy presented to site N by everything beyond N.

if nargin < 6 || isempty(gamma)
    gamma = zeros(N,1);
end
if nargin < 7 || isempty(SigmaAfterN)
    SigmaAfterN = zeros(size(p));
end

gamma = gamma(:);
if numel(gamma) ~= N
    error('gamma must contain N entries.');
end

Sigma = SigmaAfterN;
if N==0, SigmaPort=Sigma;return;end
for n = N-1:-1:1
    Sigma = kappa(n)^2 ./ ...
        (p + 1i*omega(n+1) + gamma(n+1)/2 + Sigma);
end

SigmaPort = g^2 ./ ...
    (p + 1i*omega(1) + gamma(1)/2 + Sigma);
end
