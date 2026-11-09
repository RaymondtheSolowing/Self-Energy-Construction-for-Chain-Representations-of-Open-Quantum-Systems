function [H,gamma,v]=mc_environment_matrix(N,g,omega,kappa,mc)
%MC_ENVIRONMENT_MATRIX Environment Hamiltonian, losses and system couplings.
% Basis: chain sites 1,...,N followed by the M damped closure modes of mc,
% which are all coupled to site N (to the qubit at N=0).
% H_int(system,bath)=sum_j(v_j sigma_+ b_j + conj(v_j)b_j^dagger sigma_-).
% mc is a closure model (mc_published_model or chi_mc_model), used as given.
M=mc.M;H=zeros(N+M);gamma=[zeros(N,1);mc.gamma(:)];v=zeros(N+M,1);
aux=N+(1:M);H(aux,aux)=mc.h;
if N==0
    v=mc.c(:);
else
    H(1:N,1:N)=diag(omega(1:N));
    for j=1:N-1
        H(j,j+1)=kappa(j);H(j+1,j)=kappa(j);
    end
    H(N,aux)=mc.c(:).';H(aux,N)=conj(mc.c(:));v(1)=g;
end
assert(norm(H-H','fro')<1e-13 && all(gamma>=0));
end
