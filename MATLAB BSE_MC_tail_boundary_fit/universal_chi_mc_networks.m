function universal=universal_chi_mc_networks(cfg,root)
%UNIVERSAL_CHI_MC_NETWORKS MC-form network fitted to the semicircle by minimizing chi.
% The published MC with M=cfg.M modes is the starting point and keeps its
% structure: all mode frequencies at omegaInf (h=0), nearest-neighbour
% hoppings, local losses, and couplings of every mode to the attachment site.
% Only the fit metric changes, from the time-domain correlation error to the
% band-interior rms of the cut-invariant error. The fit is done once, in units
% of kappaInf, and does not use J(omega). It is the starting point of the
% fits at each site N in chi_mc_model.
if nargin<2,root=fileparts(mfilename('fullpath'));end
path=fullfile(root,'cache','chi_mc_universal.mat');o=cfg.chiFit;
key=struct('M',cfg.M,'bandFraction',cfg.chiBandFraction,'outOfBandWeight',o.outOfBandWeight, ...
    'outOfBandRange',o.outOfBandRange,'gridIn',o.gridIn,'gridOut',o.gridOut, ...
    'maxIter',o.maxIterUniversal,'tol',o.tol,'sgMin',o.sgMin);
if cfg.reuseChiFitCache && isfile(path)
    try
        s=load(path);
        if isfield(s,'key') && isequaln(s.key,key)
            universal=s.universal;fprintf('Reused the universal chi-fitted MC network (unchanged settings).\n');return
        end
    catch
        fprintf('Could not read %s; refitting the universal network.\n',path);
    end
end
f=cfg.chiBandFraction;M=cfg.M;
xIn=linspace(-2*f,2*f,o.gridIn).';u=linspace(o.outOfBandRange(1),o.outOfBandRange(2),o.gridOut).';
xOut=[-flipud(u);u];
sIn=sigma_uniform_tail(-1i*xIn,0,1);sOut=sigma_uniform_tail(-1i*xOut,0,1);
opts=struct('fixOnsite',true,'outOfBandWeight',o.outOfBandWeight,'sgMin',o.sgMin, ...
    'maxIter',o.maxIterUniversal,'tol',o.tol);
net0=mc_published_network(M);chi0=cut_invariant_error(mc_network_response(xIn,net0),sIn);
[net,info]=fit_mc_network_chi(xIn,sIn,xOut,sOut,net0,opts);
universal=struct('M',M,'net',net,'chiRms',info.chiRms,'chiMax',info.chiMax, ...
    'publishedChiRms',sqrt(mean(chi0.^2)),'publishedChiMax',max(chi0),'info',info);
fprintf('Universal chi fit M=%d: chi_rms %.4g (published %.4g), %d iterations.\n', ...
    M,info.chiRms,universal.publishedChiRms,info.iterations);
if ~isfolder(fileparts(path)),mkdir(fileparts(path));end
save(path,'universal','key','-v7');
end
