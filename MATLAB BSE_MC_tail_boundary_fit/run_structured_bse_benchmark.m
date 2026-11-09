function results=run_structured_bse_benchmark(stage,overrides)
%RUN_STRUCTURED_BSE_BENCHMARK Five closures of a chain with N explicit modes.
% Closures at site N: asymptotic BSE, chain-extended BSE, semicircle-fit
% (published) MC, MC chi fit to the chain-extended BSE closure, and MC chi fit
% to the exact tail. Every MC closure keeps the exact chain up to site N and
% adds cfg.M damped modes coupled to site N. N runs over cfg.Nscan.
% The hard cutoff at N=cfg.Nmain appears only in Figure 1.
% Stages: 'all'; 'fit' (chain coefficients and MC fits); 'evaluate' (uses the
% saved fits); 'plot' (figures from the saved results).
if nargin<1,stage='all';end
if nargin<2,overrides=struct();end
assert(ismember(stage,{'all','fit','evaluate','plot'}),'Unknown stage "%s".',stage);
root=fileparts(mfilename('fullpath'));addpath(root,'-begin');cfg=benchmark_settings(overrides);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
fitPath=fullfile(cfg.outputDir,'fitted_environments.mat');resultPath=fullfile(cfg.outputDir,'benchmark_results.mat');
if strcmp(stage,'plot')
    saved=load(resultPath,'results');results=saved.results;
    results.cfg.figureVisible=cfg.figureVisible;results.cfg.outputDir=cfg.outputDir;
    % Refresh the legend labels of the saved results, so they match closure_labels.
    [results.methodNames,results.closureComparison.targetNames,results.closureComparison.networkNames]=closure_labels();
    save(resultPath,'results','-v7');plot_benchmark_figures(results);return
end
clockStart=tic;
if ismember(stage,{'all','fit'})
    fit=fit_closures(cfg,root);fit.fitSeconds=toc(clockStart);save(fitPath,'fit','-v7');
    if strcmp(stage,'fit'),results=fit;return;end
else
    saved=load(fitPath,'fit');fit=saved.fit;check=cfg;
    check.outputDir=fit.cfg.outputDir;check.figureVisible=fit.cfg.figureVisible;
    check.reuseChiFitCache=fit.cfg.reuseChiFitCache;
    assert(isequaln(check,fit.cfg),'Settings differ from the saved fits. Run fit or all.');fit.cfg=cfg;
end
results=evaluate_closures(fit,cfg);
results.closureComparison=closure_comparison(fit,cfg);
results.totalSeconds=toc(clockStart);
save(resultPath,'results','-v7');export_benchmark_tables(results);plot_benchmark_figures(results);
fprintf('Finished: %s\n',cfg.outputDir);
end

function fit=fit_closures(cfg,root)
% Chain coefficients, exact references, chain extensions, and the three MC closures at every N.
omegaInf=(cfg.bath.omegaL+cfg.bath.omegaR)/2;kappaInf=(cfg.bath.omegaR-cfg.bath.omegaL)/4;
[g,omega,kappa,chainInfo]=structured_chain_coeffs(max(cfg.Nscan)+1,cfg.bath,cfg.chainQuadN);
e=struct('g',g,'omega',omega,'kappa',kappa,'omegaInf',omegaInf,'kappaInf',kappaInf,'mcPortLink',cfg.mcPortLink);
w=linspace(cfg.bath.omegaL,cfg.bath.omegaR,cfg.spectralGridN);
exact=continuum_self_energy_regularized(cfg.etaSpec-1i*w,cfg.bath,cfg.continuumSpecQuadN);
fprintf('Closures at N=%d,...,%d chain modes; MC closures with M=%d damped modes at site N.\n',min(cfg.Nscan),max(cfg.Nscan),cfg.M);
targets=prepare_mc_chi_targets(cfg,e);universal=universal_chi_mc_networks(cfg,root);
n=numel(cfg.Nscan);extension=cell(n,1);models=cell(n,3);
for j=1:n
    N=cfg.Nscan(j);
    extension{j}=chain_only_bse_parameters(N,omega(1:N+1),kappa(1:N),omegaInf,kappaInf);
    assert(extension{j}.nonnegativeBath,'The chain-extended BSE has a nonpositive bath at N=%d.',N);
    models{j,1}=mc_published_model(N,cfg.M,e);
    models{j,2}=chi_mc_model('extended',N,e,universal,targets,cfg);
    models{j,3}=chi_mc_model('exact',N,e,universal,targets,cfg,models{j,2}.network);
    fprintf('N=%d: chi_rms of the fit to the chain-extended BSE closure %.4g, to the exact tail %.4g\n', ...
        N,models{j,2}.fitInfo.chiRms,models{j,3}.fitInfo.chiRms);
end
fit=e;fit.cfg=cfg;fit.chainInfo=chainInfo;fit.omegaSpec=w;fit.SigmaExact=exact;
fit.chiTargets=targets;fit.universalChiMC=universal;fit.extension=extension;fit.mcModels=models;
end

function r=evaluate_closures(fit,cfg)
% Populations, broadened port spectra and cut-invariant errors of the five closures.
g=fit.g;omega=fit.omega;kappa=fit.kappa;T=fit.chiTargets;pc=-1i*T.wIn;
w=fit.omegaSpec;p=cfg.etaSpec-1i*w;t=0:cfg.dt:cfg.tMax;op=cfg.compositeOpts;
fprintf('Computing the exact continuum reference and the closure dynamics...\n');
[Aexact,exactInfo]=continuum_amplitude_composite(t,cfg.epsilon,cfg.bath,op.Nquad,op.Npv,op.edgeOrder);
Aexact=Aexact(:);Pexact=abs(Aexact).^2;Jexact=real(fit.SigmaExact(:));
normP=trapz(t(:),Pexact.^2);normA=trapz(t(:),Pexact);
% Spectral error over the full band, as the published MC is compared with the
% semicircle. It stays finite at the band edges, where chi=1 for any closure, so
% only the cut-invariant error uses the band interior.
normJ=trapz(w(:),Jexact.^2);
n=numel(cfg.Nscan);K=5;curveN=unique([cfg.Nmain cfg.Ntable]);
populationError=zeros(n,K);populationMaxError=populationError;spectralError=populationError;
populationBound=populationError;amplitudeError=populationError;
chiRms=populationError;evaluationSeconds=populationError;closureRms=zeros(n,3);publishedFitRms=zeros(n,1);
bseInfo=cell(n,2);P=zeros(numel(t),numel(curveN),K);Sigma=zeros(numel(w),numel(curveN),K);
for j=1:n
    N=cfg.Nscan(j);q=fit.extension{j};exactTail=T.tail(T.depths==N,:).';
    parameters=[fit.omegaInf fit.kappaInf;q.omegaBar q.kappaBar];
    amplitudes=cell(1,K);responses=amplitudes;closures=amplitudes;
    for k=1:2
        om=parameters(k,1);kap=parameters(k,2);started=tic;
        [amplitudes{k},bseInfo{j,k}]=bse_amplitude_kquad(t,N,cfg.epsilon,g,omega,kappa,om,kap,om-2*kap,om+2*kap,cfg.bseOpts);
        evaluationSeconds(j,k)=toc(started);
        responses{k}=bse_port_self_energy(p,N,g,omega,kappa,om,kap);
        closures{k}=bse_boundary_self_energy(pc,N,omega,kappa,om,kap);
    end
    for k=3:5
        model=fit.mcModels{j,k-2};started=tic;
        [H,loss,v]=mc_environment_matrix(N,g,omega,kappa,model);
        amplitudes{k}=mc_population_amplitude(t,cfg.epsilon,H,loss,v);evaluationSeconds(j,k)=toc(started);
        responses{k}=mc_port_self_energy(p,N,g,omega,kappa,model);
        closures{k}=mc_tail_self_energy(pc,model);
    end
    c=find(curveN==N,1);
    for k=1:K
        Pk=abs(amplitudes{k}(:)).^2;dA=abs(amplitudes{k}(:)-Aexact);
        populationError(j,k)=sqrt(trapz(t(:),(Pk-Pexact).^2)/normP);
        populationMaxError(j,k)=max(abs(Pk-Pexact));
        % |P_X-P|<=2|A||dA|+|dA|^2 for any phase of the error amplitude dA.
        populationBound(j,k)=sqrt(trapz(t(:),(2*abs(Aexact).*dA+dA.^2).^2)/normP);
        amplitudeError(j,k)=sqrt(trapz(t(:),dA.^2)/normA);
        Jk=real(responses{k}(:));spectralError(j,k)=sqrt(trapz(w(:),(Jk-Jexact).^2)/normJ);
        chiRms(j,k)=chi_rms(T.wIn,cut_invariant_error(closures{k}(:),exactTail));
        if ~isempty(c),P(:,c,k)=Pk;Sigma(:,c,k)=responses{k}(:);end
    end
    % Closure error of each MC closure against its own target at site N.
    closureRms(j,:)=[chi_rms(T.wIn,cut_invariant_error(closures{3}(:),closures{1}(:))), ...
        chi_rms(T.wIn,cut_invariant_error(closures{4}(:),closures{2}(:))),chiRms(j,5)];
    % Published MC against the response it approximates, (link/kappaInf)^2*Sigma_inf.
    link=fit.mcModels{j,1}.portLink;
    publishedFitRms(j)=chi_rms(T.wIn,cut_invariant_error(closures{3}(:), ...
        (link/fit.kappaInf)^2*sigma_uniform_tail(pc(:),fit.omegaInf,fit.kappaInf)));
    fprintf('N=%d relative population errors [BSE inf, BSE bar, MC semicircle, MC chi bar, MC chi exact]: %s\n', ...
        N,mat2str(populationError(j,:),4));
end
Ahard=finite_chain_amplitude_spectral(cfg.Nmain,t,cfg.epsilon,g,omega,kappa,zeros(cfg.Nmain,1));
r=struct('cfg',cfg,'fit',fit,'t',t,'Pexact',Pexact,'exactInfo',exactInfo,'Phard',abs(Ahard(:)).^2, ...
    'methodNames',{closure_labels()}, ...
    'methodTags',{{'BSE_asymptotic','BSE_chain_extended','MC_semicircle','MC_chi_extended','MC_chi_exact'}}, ...
    'populationError',populationError,'populationMaxError',populationMaxError,'populationBound',populationBound, ...
    'amplitudeError',amplitudeError,'spectralError',spectralError, ...
    'chiRms',chiRms,'closureRms',closureRms,'publishedFitRms',publishedFitRms,'evaluationSeconds',evaluationSeconds, ...
    'bseInfo',{bseInfo},'curveN',curveN,'P',P,'Sigma',Sigma);
end

function q=closure_comparison(fit,cfg)
% Figures 8 and 9 at N=cfg.Ncut. Targets at site N: the asymptotic BSE closure
% kappa_N^2/(p+1i*omega_(N+1)+Sigma_inf), the chain-extended BSE closure with
% bar{Sigma}_N, and the exact tail Sigma_exact_(N+1 -> N). Closure error: each MC
% closure against its own target. Depth error: each BSE closure against the
% exact tail. Total error: each MC closure against the exact tail.
N=cfg.Ncut;j=find(cfg.Nscan==N,1);ext=fit.extension{j};models=fit.mcModels(j,:);
w=linspace(cfg.bath.omegaL,cfg.bath.omegaR,cfg.frequencyGridN).';w=w(2:end-1);p=-1i*w;
port=continuum_self_energy_axis(w,cfg.bath,cfg.chiPVN);
tails=exact_tail_self_energies(p,port,fit.g,fit.omega,fit.kappa,N);
targets=[bse_boundary_self_energy(p,N,fit.omega,fit.kappa,fit.omegaInf,fit.kappaInf), ...
    bse_boundary_self_energy(p,N,fit.omega,fit.kappa,ext.omegaBar,ext.kappaBar),tails(N+1,:).'];
networks=[mc_tail_self_energy(p,models{1}),mc_tail_self_energy(p,models{2}),mc_tail_self_energy(p,models{3})];
exactTail=repmat(targets(:,3),1,3);
closure=cut_invariant_error(networks,targets);depth=cut_invariant_error(targets(:,1:2),exactTail(:,1:2));
total=cut_invariant_error(networks,exactTail);
f=cfg.chiBandFraction;inside=abs(w-fit.omegaInf)<=f*2*fit.kappaInf;
bandRms=@(x)arrayfun(@(k)chi_rms(w(inside),x(inside,k)),1:size(x,2));
[~,targetNames,networkNames]=closure_labels();
q=struct('M',cfg.M,'N',N,'omega',w,'targets',targets,'networks',networks,'closure',closure,'depth',depth, ...
    'total',total,'closureRms',bandRms(closure),'depthRms',bandRms(depth),'totalRms',bandRms(total),'bandFraction',f, ...
    'targetNames',{targetNames},'networkNames',{networkNames});
fprintf('N=%d, M=%d: closure error chi_rms %s; depth error %s; total %s\n',N,cfg.M, ...
    mat2str(q.closureRms,4),mat2str(q.depthRms,4),mat2str(q.totalRms,4));
end

function [methods,targets,networks]=closure_labels()
%CLOSURE_LABELS Legend labels of the five closures, of the MC targets and of the MC networks.
methods={'Asymptotic BSE $\Sigma_\infty$','Chain-extended BSE $\bar{\Sigma}_N$','Semicircle-fit MC', ...
    'MC: $\chi$ fit to $\bar{\Sigma}_N$','MC: $\chi$ fit to exact tail'};
targets={'Asymptotic BSE $\Sigma_\infty$','Chain-extended BSE $\bar{\Sigma}_N$','Exact tail'};
networks={'Semicircle-fit MC','MC: $\chi$ fit to $\bar{\Sigma}_N$','MC: $\chi$ fit to exact tail'};
end
