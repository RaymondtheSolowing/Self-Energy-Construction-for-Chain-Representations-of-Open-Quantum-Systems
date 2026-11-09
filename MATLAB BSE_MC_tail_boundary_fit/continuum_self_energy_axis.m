function Sigma=continuum_self_energy_axis(w,bath,nPV)
%CONTINUUM_SELF_ENERGY_AXIS Exact system-port self-energy on the real axis.
% Sigma(0^+-i*w)=J(w)+1i*Delta(w), Delta(w)=(1/pi)*PV int J(v)/(w-v) dv.
% Valid inside and outside [omegaL,omegaR]; outside, J=0 and Delta is a
% regular integral, with the same pole subtraction as
% continuum_amplitude_composite. Exact band edges are excluded.
if nargin<3||isempty(nPV),nPV=256;end
L=bath.omegaL;R=bath.omegaR;w=w(:);
assert(all(isfinite(w)) && all(abs(w-L)>1e-14) && all(abs(w-R)>1e-14), ...
    'Exclude the band edges from the real-axis grid.');
[x,q]=gauss_legendre_rule(nPV,L,R);
[~,~,nc]=structured_spectral_density((L+R)/2,bath);
c=2*bath.alpha/nc;xc=bath.bumpCenter;s=bath.bumpWidth;
Delta=zeros(size(w));
for j=1:numel(w)
    v=w(j);d=x-v;gv=exp(-.5*((v-xc)/s)^2);
    df=exp(-.5*((x-xc)/s).^2).*d+v*gv.*expm1(-d.*(x+v-2*xc)/(2*s*s));
    smooth=-df./d;smooth(abs(d)<1e-13)=-gv*(1-v*(v-xc)/(s*s));
    lg=log(abs((v-L)/(R-v)));
    Delta(j)=c*(v*lg-(R-L)+bath.bumpAmplitude*(sum(q.*smooth)+v*gv*lg));
end
Sigma=structured_spectral_density(w,bath)+1i*Delta;
end
