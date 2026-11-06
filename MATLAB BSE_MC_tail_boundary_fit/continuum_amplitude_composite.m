function [A,info]=continuum_amplitude_composite(t,epsilon,bath,nQuad,nPV,nEdge)
%CONTINUUM_AMPLITUDE_COMPOSITE Edge-resolved reference for small closure errors.
% Same continuum model as exact_continuum_amplitude, with geometric endpoint
% panels resolving the hard-cutoff logarithm. No spectral renormalization.
% This is an optional high-accuracy reference, not a change to the main driver.
if nargin<4, nQuad=4800; end
if nargin<5, nPV=256; end
if nargin<6, nEdge=24; end
L=bath.omegaL; R=bath.omegaR;
u=unique([0,10.^(-14:-1),.25,.5,.75,1-10.^(-1:-1:-14),1]);
E=[]; weights=[];
for j=1:numel(u)-1
    n=max(nEdge,ceil(nQuad*(u(j+1)-u(j))));
    [x,w]=gauss_legendre_rule(n,0,1);
    a=L+(R-L)*u(j); b=L+(R-L)*u(j+1);
    v=a+(b-a)*x; w=(b-a)*w; keep=v>L & v<R;
    E=[E;v(keep)]; weights=[weights;w(keep)]; %#ok<AGROW>
end
[x,w]=gauss_legendre_rule(nPV,L,R);
[~,~,nc]=structured_spectral_density((L+R)/2,bath);
c=2*bath.alpha/nc; xc=bath.bumpCenter; s=bath.bumpWidth;
Delta=zeros(size(E));
for j=1:numel(E)
    e=E(j); d=x-e; ge=exp(-.5*((e-xc)/s)^2);
    % Stable subtraction of x*Gaussian(x)-e*Gaussian(e).
    df=exp(-.5*((x-xc)/s).^2).*d+e*ge.*expm1(-d.*(x+e-2*xc)/(2*s*s));
    v=-df./d; v(abs(d)<1e-13)=-ge*(1-e*(e-xc)/(s*s));
    lg=log((e-L)/(R-e));
    Delta(j)=c*(e*lg-(R-L)+bath.bumpAmplitude*(sum(w.*v)+e*ge*lg));
end
J=structured_spectral_density(E,bath);
rho=J./(pi*((E-epsilon-Delta).^2+J.^2));
weights=weights.*rho;
[Eb,Zb]=continuum_bound_states(epsilon,bath,10);
info=struct('continuumNodes',numel(E),'dispersionNodes',nPV,'edgeNodes',nEdge, ...
    'continuumWeight',sum(weights),'boundEnergies',Eb,'boundWeights',Zb, ...
    'totalWeight',sum(weights)+sum(Zb),'renormalized',false);
info.normalizationError=abs(info.totalWeight-1);
assert(info.normalizationError<1e-11,'Composite continuum normalization failed.');
E=[E;Eb(:)]; weights=[weights;Zb(:)];
t=t(:).'; A=zeros(size(t));
if numel(t)>1 && abs(t(1))<100*eps && max(abs(diff(t)-diff(t(1:2))))<1e-10
    phase=ones(size(E)); step=exp(-1i*E*(t(2)-t(1)));
    for j=1:numel(t), A(j)=weights.'*phase; phase=phase.*step; end
else
    for j=1:200:numel(t)
        ii=j:min(j+199,numel(t)); A(ii)=weights.'*exp(-1i*(E*t(ii)));
    end
end
end
