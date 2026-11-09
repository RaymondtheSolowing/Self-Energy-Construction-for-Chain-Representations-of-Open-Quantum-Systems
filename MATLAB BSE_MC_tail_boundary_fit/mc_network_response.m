function [R,dR]=mc_network_response(x,net)
%MC_NETWORK_RESPONSE Self-energy of an MC-form network in reference units.
% R(x)=c.'*(-1i*x*I+A)^(-1)*conj(c), A=1i*H+diag(gam)/2, with H real
% symmetric tridiagonal (diagonal h, off-diagonal t). Frequencies x are
% measured from the reference frequency in units of the reference coupling.
% dR(:,k) is dR/dq_k for q=[h;t;sg;real(c);imag(c)], where gam=sg.^2.
M=numel(net.h);x=x(:);c=net.c(:);
A=diag(1i*net.h(:)+net.gam(:)/2)+diag(1i*net.t(:),1)+diag(1i*net.t(:),-1);
[V,D]=eig(A);lambda=diag(D);
if rcond(V)>1e-10
    den=1./(-1i*x+lambda.');
    Gc=(den.*(V\conj(c)).')*V.';   % rows: (G*conj(c)).'
    Gcc=(den.*(V\c).')*V.';        % rows: (G*c).'=c.'*G, since A is complex symmetric
else
    Gc=zeros(numel(x),M);Gcc=Gc;I=eye(M);
    for f=1:numel(x)
        Z=-1i*x(f)*I+A;Gc(f,:)=(Z\conj(c)).';Gcc(f,:)=(Z\c).';
    end
end
R=Gc*c;
if nargout>1
    sg=sqrt(net.gam(:)).';
    P=Gcc.*Gc;
    dR=[-1i*P,-1i*(Gcc(:,1:M-1).*Gc(:,2:M)+Gcc(:,2:M).*Gc(:,1:M-1)), ...
        -sg.*P,Gc+Gcc,1i*(Gc-Gcc)];
end
end
