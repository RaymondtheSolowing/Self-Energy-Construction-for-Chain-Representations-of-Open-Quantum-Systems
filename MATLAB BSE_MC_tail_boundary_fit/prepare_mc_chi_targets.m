function targets=prepare_mc_chi_targets(cfg,e)
%PREPARE_MC_CHI_TARGETS Exact tail self-energies on the real axis, in the band.
% Row r of targets.tail is Sigma_exact_(N+1 -> N)(0^+-i*w) for N=targets.depths(r),
% obtained from the exact system-port response J+1i*Delta by the inverse
% recurrence of exact_tail_self_energies. Inside the band each inverse step
% preserves chi, so the transport is stable. Outside the band it amplifies
% rounding errors exponentially with depth, so no exact tail is formed there.
center=e.omegaInf;width=2*e.kappaInf;f=cfg.chiBandFraction;
wIn=linspace(center-f*width,center+f*width,cfg.chiGridN).';
assert(min(wIn)>cfg.bath.omegaL && max(wIn)<cfg.bath.omegaR);
depths=unique(cfg.Nscan(:));assert(all(depths>=1));
port=continuum_self_energy_axis(wIn,cfg.bath,cfg.chiPVN);
tails=exact_tail_self_energies(-1i*wIn,port,e.g,e.omega,e.kappa,max(depths));
targets=struct('wIn',wIn,'nIn',numel(wIn),'depths',depths,'tail',tails(depths+1,:),'port',port, ...
    'bandFraction',f,'source','Exact continuum port response J+i*Delta on the real axis, inverse Mobius recurrence');
fprintf('Cut-invariant targets: %d tail depths on %d band-interior frequencies.\n',numel(depths),numel(wIn));
end
