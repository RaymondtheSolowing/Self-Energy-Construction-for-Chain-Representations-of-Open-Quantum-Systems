function tail=mc_tail_self_energy(p,mc)
%MC_TAIL_SELF_ENERGY Closure self-energy c.'*(p*I+1i*h+Gamma/2)^(-1)*conj(c).
% mc is a closure model from mc_parameters or mc_model_from_network. This is
% the response presented to the attachment site N.
A=1i*mc.h+diag(mc.gamma/2);c=mc.c(:);sz=size(p);p=p(:);
[V,D]=eig(A);
if rcond(V)>1e-10
    tail=(1./(p+diag(D).'))*((c.'*V).'.*(V\conj(c)));
else
    tail=zeros(size(p));I=eye(mc.M);
    for j=1:numel(p),tail(j)=c.'*((p(j)*I+A)\conj(c));end
end
tail=reshape(tail,sz);
end
