function mc=chi_mc_model(kind,N,e,universal,targets,cfg,start)
%CHI_MC_MODEL MC-form closure at site N fitted by the cut-invariant error.
% The network keeps the structure of the published MC: M damped modes with
% nearest-neighbour hoppings, each coupled to site N, with couplings scaled
% by mc_attachment_coupling(N,e). All network parameters are fitted.
% kind='extended': fit to the chain-extended BSE closure at site N,
%   kappa_N^2/(p+1i*omega_(N+1)+bar{Sigma}_N(p)), on the band interior of
%   bar{Sigma}_N. It uses only omega(1:N+1) and kappa(1:N), so it does not need J(omega).
% kind='exact': fit to the exact tail Sigma_exact_(N+1 -> N) in the band,
%   obtained from J by the inverse recurrence (prepare_mc_chi_targets).
%   It needs J(omega).
% Both fits add a weak out-of-band penalty with the chain-extended BSE closure
% as its target, since the inverse recurrence is unstable outside the band.
% start: starting network in reference units; default universal.net.
assert(N>=1 && N==floor(N));
if nargin<7||isempty(start),start=universal.net;end
assert(numel(start.h)==cfg.M,'The starting network must have M=%d modes.',cfg.M);
q=chain_only_bse_parameters(N,e.omega(1:N+1),e.kappa(1:N),e.omegaInf,e.kappaInf);
link=mc_attachment_coupling(N,e);o=cfg.chiFit;f=cfg.chiBandFraction;
% Reference units: x=(w-omegaBar)/kappaBar. A network response R(x) gives the
% closure (link^2/kappaBar)*R(x), so the targets are scaled by kappaBar/link^2.
scale=q.kappaBar/link^2;
u=linspace(o.outOfBandRange(1),o.outOfBandRange(2),o.gridOut).';xOut=[-flipud(u);u];
target=@(x)scale*bse_boundary_self_energy(-1i*(q.omegaBar+q.kappaBar*x),N,e.omega,e.kappa,q.omegaBar,q.kappaBar);
switch kind
    case 'extended'
        xIn=linspace(-2*f,2*f,o.gridIn).';tIn=target(xIn);
        source='chi fit to the chain-extended BSE closure at site N';
    case 'exact'
        row=find(targets.depths==N,1);assert(~isempty(row),'No exact tail at N=%d.',N);
        xIn=(targets.wIn-q.omegaBar)/q.kappaBar;tIn=scale*targets.tail(row,:).';
        source='chi fit to the exact tail at site N';
    otherwise
        error('Unknown chi-fitted MC kind "%s".',kind);
end
opts=struct('fixOnsite',false,'outOfBandWeight',o.outOfBandWeight,'sgMin',o.sgMin, ...
    'maxIter',o.maxIterLocal,'tol',o.tol);
[net,info]=fit_mc_network_chi(xIn,tIn,xOut,target(xOut),start,opts);
mc=mc_model_from_network(net,q.omegaBar,q.kappaBar,link,source);
mc.fitInfo=info;mc.extension=q;mc.kind=kind;mc.attachmentSite=N;
end
