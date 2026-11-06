function export_table_I(r)
%EXPORT_TABLE_I Chain-extended BSE parameters and relative population errors at N=cfg.Ntable.
c=r.cfg;[~,rows]=ismember(c.Ntable,c.Nscan);
fid=fopen(fullfile(c.outputDir,'Table_I.tex'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
line('\begin{table}[t]');line('\centering');line('\small');line('\setlength{\tabcolsep}{3pt}');
line('\renewcommand{\arraystretch}{1.3}');line('');
line('\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}|c|c|c|}');line('\hline');
line('\multicolumn{3}{|c|}{\rule[-0.7ex]{0pt}{3.2ex}BSE: chain extension $\bar{\Sigma}_N$} \\');line('\hline');
line('$N$ & $\bar{\omega}_N$ & $\bar{\kappa}_N$ \\');line('\hline');
for j=rows
    q=r.fit.extension{j};line([sprintf('%d & $%.10f$ & $%.10f$',c.Nscan(j),q.omegaBar,q.kappaBar) ' \\']);
end
line('\hline');line('\end{tabular*}');line('');line('\medskip');
line('\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}|c|c|c|c|c|c|}');line('\hline');
line(['\multicolumn{6}{|c|}{\rule[-0.7ex]{0pt}{3.2ex}Relative population error $\varepsilon_X^{(P)}$, $T=' num2str(c.tMax) '$} \\']);
line('\hline');
line(['& \multicolumn{2}{c|}{BSE} & \multicolumn{3}{c|}{MC, $M=' num2str(c.M) '$} \\']);
line('$N$ & $\Sigma_\infty$ & $\bar{\Sigma}_N$ & semicircle & $\chi$: $\bar{\Sigma}$ & $\chi$: exact \\');line('\hline');
for j=rows
    entries=arrayfun(@sci,r.populationError(j,:),'UniformOutput',false);
    line([sprintf('%d',c.Nscan(j)) sprintf(' & %s',entries{:}) ' \\']);
end
line('\hline');line('\end{tabular*}');line('');
line(['\caption{Chain-extended BSE parameters and relative population errors of the five closures for $T=' num2str(c.tMax) '$.']);
line(['The MC closures use $M=' num2str(c.M) '$ damped modes coupled to site $N$.}']);
line('\label{tab:termination_parameters_errors}');line('\end{table}');
    function line(s),fprintf(fid,'%s\n',s);end
end
function s=sci(x)
if x==0,s='$0$';return;end
ex=floor(log10(abs(x)));m=round(10*x/10^ex)/10;
if abs(m)>=10,m=m/10;ex=ex+1;end
s=[sprintf('$%.1f',m) '\times10^{' num2str(ex) '}$'];
end
