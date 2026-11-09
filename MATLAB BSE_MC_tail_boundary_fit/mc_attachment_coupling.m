function link=mc_attachment_coupling(N,e)
%MC_ATTACHMENT_COUPLING Scale of the couplings of the MC modes to site N.
% e.mcPortLink='kappaN' (default): kappa_N, the attachment coupling of the
%   retained chain, so the closure takes the place of the link N -> N+1.
% e.mcPortLink='kappaInf': kappaInf, the rescaling rule of the published MC.
% At N=0 the closure is attached to the qubit and the scale is g.
mode='kappaN';if isfield(e,'mcPortLink'),mode=e.mcPortLink;end
if N==0
    link=e.g;
elseif strcmp(mode,'kappaN')
    link=e.kappa(N);
elseif strcmp(mode,'kappaInf')
    link=e.kappaInf;
else
    error('Unknown mcPortLink "%s".',mode);
end
end
