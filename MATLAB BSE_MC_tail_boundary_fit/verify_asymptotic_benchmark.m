function report=verify_asymptotic_benchmark
%VERIFY_ASYMPTOTIC_BENCHMARK Independent checks of the five-closure benchmark.
% MC closure stability and passivity, matrix-exponential propagation, refined
% time, frequency and real-axis quadratures for every N and closure, refined
% chain coefficients, cut invariance, the chi budget, and the figure data.
root=fileparts(mfilename('fullpath'));addpath(root);c=benchmark_settings();
s=load(fullfile(c.outputDir,'benchmark_results.mat'),'results');r=s.results;f=r.fit;n=numel(c.Nscan);
assert(isequal(r.methodTags,{'BSE_asymptotic','BSE_chain_extended','MC_semicircle','MC_chi_extended','MC_chi_exact'}));
report=struct('MATLAB',version,'protocolVersion',c.protocolVersion);
% 1. MC closures: Hermitian mode matrix, nonnegative losses, stable poles, passive
% response, and the single-excitation amplitude against a matrix exponential.
maxAmplitude=0;maxPole=-Inf;minReal=Inf;times=[0 1 10 50 200];
for j=1:n
    N=c.Nscan(j);
    for k=1:3
        m=f.mcModels{j,k};assert(m.M==c.M && norm(m.h-m.h','fro')<1e-12 && all(m.gamma>=0));
        maxPole=max(maxPole,max(real(eig(-1i*m.h-diag(m.gamma/2)))));
        minReal=min(minReal,min(real(mc_tail_self_energy(1e-7-1i*linspace(-4,5,4001),m))));
        [H,loss,v]=mc_environment_matrix(N,f.g,f.omega,f.kappa,m);
        [A,info]=mc_population_amplitude(times,c.epsilon,H,loss,v);
        for q=1:numel(times),U=expm(-1i*info.Hfull*times(q));maxAmplitude=max(maxAmplitude,abs(A(q)-U(1,1)));end
    end
end
assert(maxAmplitude<1e-10 && maxPole<1e-10 && minReal>-1e-10);
report.MaxExpmAmplitudeDifference=maxAmplitude;report.MaxPoleRealPart=maxPole;report.MinRealClosureSelfEnergy=minReal;
% 2. Exact continuum reference: finer quadrature and a twice finer time grid.
t=0:c.dt/2:c.tMax;op=c.compositeOpts;
[A,info]=continuum_amplitude_composite(t,c.epsilon,c.bath,2*op.Nquad,2*op.Npv,op.edgeOrder+8);
Pref=abs(A(:)).^2;report.ReferencePopulationChange=max(abs(Pref(1:2:end)-r.Pexact(:)));
report.RefinedReferenceNormalization=abs(info.totalWeight-1);assert(report.ReferencePopulationChange<1e-10);
% 3. Population errors of all closures at every N on the refined time grid, with
% a twice finer BSE spectral inversion.
opts=c.bseOpts;opts.NkMin=2*c.bseOpts.NkMin-1;opts.phaseStepMax=c.bseOpts.phaseStepMax/2;
normP=trapz(t(:),Pref.^2);relative=zeros(n,5);maximum=relative;bound=relative;
for j=1:n
    N=c.Nscan(j);q=f.extension{j};parameters=[f.omegaInf f.kappaInf;q.omegaBar q.kappaBar];
    for k=1:5
        if k<=2
            om=parameters(k,1);kap=parameters(k,2);
            Ak=bse_amplitude_kquad(t,N,c.epsilon,f.g,f.omega,f.kappa,om,kap,om-2*kap,om+2*kap,opts);
        else
            [H,loss,v]=mc_environment_matrix(N,f.g,f.omega,f.kappa,f.mcModels{j,k-2});
            Ak=mc_population_amplitude(t,c.epsilon,H,loss,v);
        end
        Pk=abs(Ak(:)).^2;relative(j,k)=sqrt(trapz(t(:),(Pk-Pref).^2)/normP);maximum(j,k)=max(abs(Pk-Pref));
        dA=abs(Ak(:)-A(:));bound(j,k)=sqrt(trapz(t(:),(2*abs(A(:)).*dA+dA.^2).^2)/normP);
    end
    if mod(j,10)==0,fprintf('Verified the populations through N=%d.\n',N);end
end
report.RefinedRelativePopulationErrors=relative;
report.MaxRelativePopulationErrorChange=max(abs(relative-r.populationError)./r.populationError,[],'all');
report.MaxMaxPopulationErrorChange=max(abs(maximum-r.populationMaxError)./r.populationMaxError,[],'all');
report.MaxPopulationBoundChange=max(abs(bound-r.populationBound)./r.populationBound,[],'all');
assert(report.MaxRelativePopulationErrorChange<.02 && report.MaxMaxPopulationErrorChange<.02 && report.MaxPopulationBoundChange<.02);
% The bound holds pointwise in time, so it bounds the relative population error.
assert(all(r.populationBound>=r.populationError,'all'));
% 4. Spectral errors with a twice finer frequency grid and continuum quadrature.
% The spectral error uses the full band; the cut-invariant error uses the band interior.
w=linspace(c.bath.omegaL,c.bath.omegaR,2*c.spectralGridN-1);p=c.etaSpec-1i*w;
Jref=real(continuum_self_energy_regularized(p,c.bath,2*c.continuumSpecQuadN));Jref=Jref(:);
normJ=trapz(w(:),Jref.^2);spectral=zeros(n,5);
for j=1:n
    N=c.Nscan(j);q=f.extension{j};
    responses={bse_port_self_energy(p,N,f.g,f.omega,f.kappa,f.omegaInf,f.kappaInf), ...
        bse_port_self_energy(p,N,f.g,f.omega,f.kappa,q.omegaBar,q.kappaBar)};
    for k=3:5,responses{k}=mc_port_self_energy(p,N,f.g,f.omega,f.kappa,f.mcModels{j,k-2});end
    for k=1:5,spectral(j,k)=sqrt(trapz(w(:),(real(responses{k}(:))-Jref).^2)/normJ);end
end
report.MaxSpectralErrorChange=max(abs(spectral-r.spectralError)./r.spectralError,[],'all');
assert(report.MaxSpectralErrorChange<.02);
% 5. Chain coefficients from a twice finer quadrature of the spectral measure.
[g2,om2,kap2]=structured_chain_coeffs(max(c.Nscan)+1,c.bath,2*c.chainQuadN);
report.ChainCoefficientChange=max([abs(g2-f.g);abs(om2-f.omega);abs(kap2-f.kappa)]);
assert(report.ChainCoefficientChange<1e-10);
% 6. Cut-invariant errors on a twice finer real-axis grid with a finer principal-value rule.
wf=linspace(min(f.chiTargets.wIn),max(f.chiTargets.wIn),2*c.chiGridN-1).';pf=-1i*wf;
tails=exact_tail_self_energies(pf,continuum_self_energy_axis(wf,c.bath,2*c.chiPVN),f.g,f.omega,f.kappa,max(c.Nscan));
chi=zeros(n,5);
for j=1:n
    N=c.Nscan(j);q=f.extension{j};reference=tails(N+1,:).';
    closures={bse_boundary_self_energy(pf,N,f.omega,f.kappa,f.omegaInf,f.kappaInf), ...
        bse_boundary_self_energy(pf,N,f.omega,f.kappa,q.omegaBar,q.kappaBar)};
    for k=3:5,closures{k}=mc_tail_self_energy(pf,f.mcModels{j,k-2});end
    for k=1:5,chi(j,k)=chi_rms(wf,cut_invariant_error(closures{k}(:),reference));end
end
report.MaxChiRmsChange=max(abs(chi-r.chiRms)./r.chiRms,[],'all');assert(report.MaxChiRmsChange<.02);
% 7. Cut invariance at N=Ncut: chi between a closure and the exact tail at the cut
% equals chi between the two responses at the system port.
q=r.closureComparison;inside=abs(q.omega-f.omegaInf)<=q.bandFraction*2*f.kappaInf;wq=q.omega(inside);
portExact=continuum_self_energy_axis(wq,c.bath,c.chiPVN);invariance=0;
closures=[q.targets(inside,1:2),q.networks(inside,:)];
for k=1:5
    atCut=cut_invariant_error(closures(:,k),q.targets(inside,3));
    port=boundary_self_energy_chain(-1i*wq,q.N,f.g,f.omega,f.kappa,zeros(q.N,1),closures(:,k));
    invariance=max(invariance,max(abs(cut_invariant_error(port,portExact)./atCut-1)));
end
report.ChiCutPortInvarianceMaxRelativeDifference=invariance;assert(invariance<1e-6);
% 8. Chi budget at N=Ncut for the semicircle-fit MC and the chi fit to bar{Sigma}, at every
% band-interior frequency: |chi_M-chi_N|/(1-chi_M*chi_N)<=chi_tot<=(chi_M+chi_N)/(1+chi_M*chi_N).
violation=-Inf;
for k=1:2
    cm=q.closure(inside,k);cn=q.depth(inside,k);ct=q.total(inside,k);
    violation=max([violation;abs(cm-cn)./(1-cm.*cn)-ct;ct-(cm+cn)./(1+cm.*cn)]);
end
report.ChiBudgetMaxViolation=violation;assert(violation<1e-12);
% 9. The published MC against the response it approximates, (link/kappaInf)^2*Sigma_inf.
% Chi is unchanged by the rescaling, so this error is the same at every N.
report.PublishedUniversalFitChiRms=r.publishedFitRms(1);
report.PublishedUniversalFitSpreadOverN=(max(r.publishedFitRms)-min(r.publishedFitRms))/min(r.publishedFitRms);
assert(report.PublishedUniversalFitSpreadOverN<1e-10);
% 10. Figures and the data of Figure 3.
expected={'Figure1_time_domain','Figure2_effective_spectrum','Figure3_population_error_vs_N', ...
    'Figure4_spectral_error_vs_N','Figure5_total_error_vs_N','Figure6_depth_error_vs_N', ...
    'Figure7_closure_error_vs_N','Figure8_self_energy_at_cut','Figure9_total_error_at_cut'};
found=cellfun(@(x)isfile(fullfile(c.outputDir,[x '.pdf'])),expected);
report.FigureCount=sum(found);assert(all(found),'Missing figures; run the plot stage.');
fig=openfig(fullfile(c.outputDir,'Figure3_population_error_vs_N.fig'),'invisible');
lines=findobj(fig,'Type','line');matched=0;
for k=1:numel(lines)
    column=find(strcmp(lines(k).DisplayName,r.methodNames),1);
    if ~isempty(column)
        assert(isequal(lines(k).XData(:),c.Nscan(:)) && isequal(lines(k).YData(:),r.populationError(:,column)));
        matched=matched+1;
    end
end
close(fig);assert(matched==5);
report.Passed=true;
fid=fopen(fullfile(c.outputDir,'verification_asymptotic.json'),'w');fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));fclose(fid);
disp(report);
end
