function [amplitude,info]=mc_population_amplitude(t,epsilon,H,gamma,v)
%MC_POPULATION_AMPLITUDE Single-excitation validation of a selected bath.
% This routine is never called by the parameter-selection code.
N=size(H,1); Hfull=[epsilon,v.';conj(v),H]-1i*diag([0;gamma])/2;
[V,D]=eig(Hfull); e0=zeros(N+1,1); e0(1)=1;
info.eigenvectorCondition=cond(V);
assert(info.eigenvectorCondition<1e10,'Ill-conditioned eigenbasis; use a matrix-exponential propagator.');
weights=V(1,:).' .* (V\e0); energies=diag(D);
t=t(:).'; amplitude=zeros(size(t));
for j=1:500:numel(t)
    ii=j:min(j+499,numel(t));
    amplitude(ii)=weights.'*exp(-1i*energies*t(ii));
end
info.Hfull=Hfull;
info.initialError=abs(sum(weights)-1);
end
