function SigmaPort = bse_port_self_energy(p,N,g,omega,kappa,omegaBar,kappaBar)
%BSE_PORT_SELF_ENERGY Fitted bar{Sigma}_N propagated to the system port.
% N counts explicit sites 1,...,N. The first omitted site N+1 is eliminated
% analytically into the boundary response; the continuation beyond N+1 is
% the homogeneous infinite tail. This agrees with green_bse.
% Does NOT include the free qubit denominator p+i*epsilon.

if N < 1 || N ~= floor(N) || numel(omega)<N+1 || numel(kappa)<N
    error('Need integer N>=1, omega(1:N+1), and kappa(1:N).');
end
SigmaAfterN = bse_boundary_self_energy(p,N,omega,kappa,omegaBar,kappaBar);
SigmaPort = boundary_self_energy_chain( ...
    p,N,g,omega,kappa,zeros(N,1),SigmaAfterN);
end
