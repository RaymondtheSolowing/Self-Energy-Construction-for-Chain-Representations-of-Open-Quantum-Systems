function [net,info]=fit_mc_network_chi(xIn,targetIn,xOut,targetOut,net0,opts)
%FIT_MC_NETWORK_CHI Fit an MC-form network by minimizing the cut-invariant error.
% Objective: mean of chi(x)^2 over the in-band grid xIn, plus the weak
% absolute penalty |R-target|^2, weighted by opts.outOfBandWeight^2, on
% the out-of-band grid xOut. All quantities are in reference units: x is
% measured from the reference frequency in units of the reference coupling,
% and the targets are self-energies divided by that coupling.
% opts.fixOnsite=true keeps h=net0.h (the published MC has h=0).
if nargin<6,opts=struct();end
fixOnsite=get_opt(opts,'fixOnsite',false);lam=get_opt(opts,'outOfBandWeight',.3);
sgMin=get_opt(opts,'sgMin',1e-2);
M=numel(net0.h);xIn=xIn(:);xOut=xOut(:);targetIn=targetIn(:);targetOut=targetOut(:);
sg=sqrt(max(net0.gam(:),0));sg(sg<sgMin)=sgMin;
qFull=[net0.h(:);net0.t(:);sg;real(net0.c(:));imag(net0.c(:))];
free=true(5*M-1,1);if fixOnsite,free(1:M)=false;end
scale=1/sqrt(numel(xIn));
fun=@(q)residuals(q);
[qFree,lm]=lm_least_squares(fun,qFull(free),opts);
qFull(free)=qFree;net=unpack(qFull,M);
chi=cut_invariant_error(mc_network_response(xIn,net),targetIn);
info=struct('chiRms',sqrt(mean(chi.^2)),'chiMax',max(chi),'cost',lm.cost, ...
    'initialCost',lm.initialCost,'iterations',lm.iterations,'evaluations',lm.evaluations);
    function [r,Jr]=residuals(qf)
        qq=qFull;qq(free)=qf;n=unpack(qq,M);
        [R,dR]=mc_network_response([xIn;xOut],n);
        sgn=sign(qq(2*M:3*M-1));sgn(sgn==0)=1;
        dR(:,2*M:3*M-1)=dR(:,2*M:3*M-1).*sgn.';
        k=numel(xIn);Ri=R(1:k);Ro=R(k+1:end);
        s=Ri+conj(targetIn);D=abs(s);e=(Ri-targetIn)./D;o=lam*(Ro-targetOut);
        r=scale*[real(e);imag(e);real(o);imag(o)];
        dEr=1./D-(Ri-targetIn).*real(s)./D.^3;dEi=1i./D-(Ri-targetIn).*imag(s)./D.^3;
        dE=dEr.*real(dR(1:k,:))+dEi.*imag(dR(1:k,:));dO=lam*dR(k+1:end,:);
        Jr=scale*[real(dE);imag(dE);real(dO);imag(dO)];Jr=Jr(:,free);
    end
end

function net=unpack(q,M)
net=struct('h',q(1:M),'t',q(M+1:2*M-1),'gam',q(2*M:3*M-1).^2, ...
    'c',q(3*M:4*M-1)+1i*q(4*M:5*M-1));
end

function value=get_opt(opts,name,default)
if isfield(opts,name),value=opts.(name);else,value=default;end
end
