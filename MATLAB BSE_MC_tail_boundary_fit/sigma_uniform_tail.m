function [SigmaInf,dSigmaInf] = sigma_uniform_tail(p,omegaInf,kappaInf)
%SIGMA_UNIFORM_TAIL Physical surface response of a uniform half-chain.

if ~isscalar(kappaInf) || ~isreal(kappaInf) || kappaInf <= 0
    error('kappaInf must be a positive real scalar.');
end

z = p + 1i*omegaInf;
rootTerm = sqrt(z.^2 + 4*kappaInf^2);

left = real(z) < 0;
rootTerm(left) = -rootTerm(left);
onAxis = real(z) == 0;
wrongImag = onAxis & (imag(rootTerm).*imag(z) < 0);
rootTerm(wrongImag) = -rootTerm(wrongImag);

SigmaInf = 2*kappaInf^2 ./ (z + rootTerm);
if nargout > 1
    dSigmaInf = -SigmaInf ./ rootTerm;
end
end
