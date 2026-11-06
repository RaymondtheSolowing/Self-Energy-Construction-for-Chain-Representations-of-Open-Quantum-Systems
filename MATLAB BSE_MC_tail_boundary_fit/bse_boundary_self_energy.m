function response=bse_boundary_self_energy(p,N,omega,kappa,omegaBar,kappaBar)
%BSE_BOUNDARY_SELF_ENERGY Complete fitted response attached to explicit site N.
% bar{Sigma}_N is the continuation beyond the eliminated site N+1.
response=kappa(N)^2./(p+1i*omega(N+1)+sigma_uniform_tail(p,omegaBar,kappaBar));
end
