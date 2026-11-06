function Sigma=mc_port_self_energy(p,N,g,omega,kappa,mc)
%MC_PORT_SELF_ENERGY MC closure attached to site N, propagated to the system port.
% N exact chain sites plus the M damped closure modes of mc.
assert(N>=0 && N==floor(N));
Sigma=boundary_self_energy_chain(p,N,g,omega,kappa,zeros(N,1),mc_tail_self_energy(p,mc));
end
