function export_benchmark_tables(r)
%EXPORT_BENCHMARK_TABLES CSV data of Figures 1-9, the MC closures, Table I and the settings.
c=r.cfg;out=c.outputDir;tags=r.methodTags;N=c.Nscan(:);ext=r.fit.extension;
% Figures 3-5: errors of the five closures versus N (ChiRms is the total cut-invariant error
% against the exact tail, Figure 5), with the population error bound and the relative
% amplitude error.
s=table(N,cellfun(@(x)x.omegaBar,ext),cellfun(@(x)x.kappaBar,ext),cellfun(@(x)x.a,ext), ...
    cellfun(@(x)x.b,ext),cellfun(@(x)x.sampleCount,ext),cellfun(@(x)x.nonnegativeBath,ext), ...
    'VariableNames',{'N','OmegaBar','KappaBar','aCorrection','bCorrection','FitSampleCount','NonnegativeBath'});
for k=1:numel(tags)
    s.([tags{k} '_RelativePopulationError'])=r.populationError(:,k);
    s.([tags{k} '_PopulationErrorBound'])=r.populationBound(:,k);
    s.([tags{k} '_MaxPopulationError'])=r.populationMaxError(:,k);
    s.([tags{k} '_RelativeAmplitudeError'])=r.amplitudeError(:,k);
    s.([tags{k} '_RelativeSpectralError'])=r.spectralError(:,k);
    s.([tags{k} '_ChiRms'])=r.chiRms(:,k);
    s.([tags{k} '_EvaluationSeconds'])=r.evaluationSeconds(:,k);
end
writetable(s,fullfile(out,'closure_errors_vs_N.csv'));
% Figures 6-7 and the error budget at every N: closure errors of the MC closures against their own
% targets, depth errors of the BSE closures, and total errors against the exact tail.
b=table(N,r.closureRms(:,1),r.closureRms(:,2),r.closureRms(:,3),r.chiRms(:,1),r.chiRms(:,2), ...
    r.chiRms(:,3),r.chiRms(:,4),r.chiRms(:,5),r.publishedFitRms, ...
    'VariableNames',{'N','MC_semicircle_ClosureChiRms','MC_chi_extended_ClosureChiRms','MC_chi_exact_ClosureChiRms', ...
    'BSE_asymptotic_DepthChiRms','BSE_chain_extended_DepthChiRms','MC_semicircle_TotalChiRms', ...
    'MC_chi_extended_TotalChiRms','MC_chi_exact_TotalChiRms','MC_semicircle_UniversalFitChiRms'});
writetable(b,fullfile(out,'chi_budget_vs_N.csv'));
% Figures 1-2 and the curves at the other Table I values of N.
t=table(r.t(:),r.Pexact(:),r.Phard(:),'VariableNames',{'t','Reference',sprintf('HardCutoff_N%d',c.Nmain)});
w=r.fit.omegaSpec(:);sp=table(w,real(r.fit.SigmaExact(:)),imag(r.fit.SigmaExact(:)), ...
    'VariableNames',{'omega','Reference_Jeta','Reference_ImagSigma'});
for j=1:numel(r.curveN)
    for k=1:numel(tags)
        name=sprintf('%s_N%d',tags{k},r.curveN(j));t.(name)=r.P(:,j,k);
        sp.([name '_Jeta'])=real(r.Sigma(:,j,k));sp.([name '_ImagSigma'])=imag(r.Sigma(:,j,k));
    end
end
writetable(t,fullfile(out,'time_domain_curves.csv'));writetable(sp,fullfile(out,'spectral_curves.csv'));
% Figures 8-9. Figure 9 plots the TotalError columns of all five closures.
q=r.closureComparison;tn={'BSE_asymptotic','BSE_chain_extended','ExactTail'};nn=tags(3:5);
s=table(q.omega,'VariableNames',{'omega'});
for k=1:3
    s.([tn{k} '_Real'])=real(q.targets(:,k));s.([tn{k} '_Imag'])=imag(q.targets(:,k));
end
for k=1:3
    s.([nn{k} '_Real'])=real(q.networks(:,k));s.([nn{k} '_Imag'])=imag(q.networks(:,k));
    s.([nn{k} '_ClosureError'])=q.closure(:,k);s.([nn{k} '_TotalError'])=q.total(:,k);
end
s.BSE_asymptotic_DepthError=q.depth(:,1);s.BSE_chain_extended_DepthError=q.depth(:,2);
% For the BSE closures the total error against the exact tail equals the depth error.
s.BSE_asymptotic_TotalError=q.depth(:,1);s.BSE_chain_extended_TotalError=q.depth(:,2);
writetable(s,fullfile(out,sprintf('frequency_comparison_at_N%d.csv',q.N)));
% Physical parameters of the three MC closures at every N: mode frequencies,
% loss rates, hoppings to the next mode, and couplings to site N.
rows=[];
for j=1:numel(N)
    for k=1:3
        m=r.fit.mcModels{j,k};M=m.M;
        rows=[rows;[repmat([N(j) k],M,1),(1:M).',m.omega(:),m.gamma(:),[m.hopping(:);NaN], ...
            real(m.c(:)),imag(m.c(:)),repmat(m.portLink,M,1)]]; %#ok<AGROW>
    end
end
writetable(array2table(rows,'VariableNames',{'N','Closure','Mode','Omega','Gamma','HoppingToNext', ...
    'ReCoupling','ImCoupling','AttachmentCoupling'}),fullfile(out,'mc_closure_networks.csv'));
% Universal chi fit to the semicircle, in units of kappaInf with frequencies from omegaInf.
u=r.fit.universalChiMC;n=u.net;M=u.M;
writetable(array2table([(1:M).',n.h(:),[n.t(:);NaN],n.gam(:),real(n.c(:)),imag(n.c(:)), ...
    repmat([u.chiRms u.publishedChiRms],M,1)],'VariableNames',{'Mode','h','t','gamma','ReC','ImC', ...
    'ChiRms','PublishedChiRms'}),fullfile(out,'chi_mc_universal_network.csv'));
fid=fopen(fullfile(out,'settings.json'),'w');fprintf(fid,'%s\n',jsonencode(c,PrettyPrint=true));fclose(fid);
export_table_I(r);
end
