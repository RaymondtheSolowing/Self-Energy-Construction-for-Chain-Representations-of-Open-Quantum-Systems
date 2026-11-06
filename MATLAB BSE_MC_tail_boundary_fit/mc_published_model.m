function mc=mc_published_model(N,M,e)
%MC_PUBLISHED_MODEL Published (semicircle-fit) MC attached to site N.
% The published universal coefficients are shifted and rescaled to the bath:
% internal frequencies, hoppings and losses use omegaInf and kappaInf, and the
% couplings to site N are mc_attachment_coupling(N,e)*c. With the default
% 'kappaN', the closure approximates kappa_N^2/(p+1i*omegaInf+Sigma_inf(p)).
mc=mc_parameters(M,e.omegaInf,e.kappaInf,mc_attachment_coupling(N,e));
mode='kappaN';if isfield(e,'mcPortLink'),mode=e.mcPortLink;end
mc.portLinkMode=mode;mc.attachmentSite=N;mc.kind='published';
end
