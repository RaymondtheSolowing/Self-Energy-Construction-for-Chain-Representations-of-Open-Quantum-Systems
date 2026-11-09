function plot_benchmark_figures(r,figureNumbers)
%PLOT_BENCHMARK_FIGURES Figures 1-9 of the five-closure benchmark.
% 1: populations at N=Nmain with the hard cutoff and the exact continuum.
% 2: broadened port spectra at N=Nmain, with zoomed panels of the peak and the upper band edge.
% 3-4: relative population error and relative spectral error versus N.
% 5: total cut-invariant error of the five closures against the exact tail versus N.
% 6: depth error of the BSE closures versus N; 7: closure error of the MC closures versus N.
% 8: MC closures and their targets at N=Ncut.
% 9: total cut-invariant error of the five closures against the exact tail at N=Ncut.
% Legends use place_legend. In Figures 3-9, which are panels in the manuscript, the legends
% are centered and at most panel_legend_width of the axes wide, so the top-left corner stays
% free for the panel labels (a), (b), ...
if nargin<2,figureNumbers=1:9;end
c=r.cfg;out=c.outputDir;s=closure_styles(r.methodNames);main=find(r.curveN==c.Nmain,1);
if ismember(1,figureNumbers)
    f=make(c,'Time-domain comparison',[720 500]);ax=axes(f);hold(ax,'on');mask=r.t<=c.tFigure1Max;
    t=r.t(mask);t=t(:);Pexact=r.Pexact(mask);Phard=r.Phard(mask);P=reshape(r.P(mask,main,:),[],5);
    hard=plot(ax,t,Phard,'-','Color',[.72 .72 .72],'LineWidth',2.6,'DisplayName','Hard cutoff');
    h=plot(ax,t,Pexact,'k','LineWidth',2,'DisplayName','Exact continuum');h=[h hard];
    for k=1:5,h=[h plot_closure(ax,t,P(:,k),s,k)];end %#ok<AGROW>
    xlabel(ax,'t');ylabel(ax,'P_e(t)');xlim(ax,[0 c.tFigure1Max]);ylim(ax,[0 1.02]);finish(ax);
    % Top center: after the initial decay, no curve rises above P_e=0.4 below t=90.
    place_legend(ax,h,'north',1,t,[Pexact(:) Phard(:) P],.6);
    save_plot(f,out,'Figure1_time_domain');
end
if ismember(2,figureNumbers)
    % Main panel and two zoomed panels: the peak of J and the upper band edge.
    f=make(c,'Broadened effective spectrum',[720 720]);w=r.fit.omegaSpec(:);J=real(r.fit.SigmaExact(:));
    curves=real(reshape(r.Sigma(:,main,:),numel(w),5));
    ax=axes('Parent',f,'Position',[.115 .465 .85 .515]);hold(ax,'on');
    h=plot(ax,w,J,'k','LineWidth',2.4,'DisplayName',sprintf('Exact continuum, $\\eta=%.4g$',c.etaSpec));
    for k=1:5,h=[h plot_closure(ax,w,curves(:,k),s,k)];end %#ok<AGROW>
    xlabel(ax,'\omega');ylabel(ax,'$J_{\mathrm{eff}}^{(\eta)}(\omega)$','Interpreter','latex');
    xlim(ax,[c.bath.omegaL c.bath.omegaR]);ylim(ax,[0 1.1*max(J)]);finish(ax);
    % Upper left, narrower than the rise of the peak.
    place_legend(ax,h,'northwest',1,w,[J curves],.56);
    [~,peak]=max(J);windows=[w(peak)+[-.015 .015];c.bath.omegaR-.05 c.bath.omegaR];
    positions=[.115 .075 .375 .27;.59 .075 .375 .27];
    for z=1:2
        in=w>=windows(z,1) & w<=windows(z,2);values=[J(in) curves(in,:)];
        low=min(values(:));high=max(values(:));pad=.06*(high-low);
        plot(ax,windows(z,[1 2 2 1 1]),[low-pad low-pad high+pad high+pad low-pad],':','Color',[.45 .45 .45],'LineWidth',.8,'HandleVisibility','off');
        zx=axes('Parent',f,'Position',positions(z,:));hold(zx,'on');
        plot(zx,w(in),J(in),'k','LineWidth',2.4);
        for k=1:5,plot_closure(zx,w(in),curves(in,k),s,k);end
        xlim(zx,windows(z,:));ylim(zx,[low-pad high+pad]);xlabel(zx,'\omega');finish(zx);
        xticks(zx,ceil(100*windows(z,1))/100:.01:windows(z,2));set(zx,'FontSize',12);
        if z==1,ylabel(zx,'$J_{\mathrm{eff}}^{(\eta)}(\omega)$','Interpreter','latex','FontSize',14);end
    end
    save_plot(f,out,'Figure2_effective_spectrum');
end
if ismember(3,figureNumbers)
    f=errors_vs_N(r,r.populationError,'Relative population error $\varepsilon_X^{(P)}$',s);
    save_plot(f,out,'Figure3_population_error_vs_N');
end
if ismember(4,figureNumbers)
    f=errors_vs_N(r,r.spectralError,'Relative spectral error $\varepsilon_X^{(J)}$',s);
    save_plot(f,out,'Figure4_spectral_error_vs_N');
end
if ismember(5,figureNumbers)
    % Total error: each closure against the exact tail at site N. For the BSE closures it
    % equals the depth error of Figure 6.
    f=errors_vs_N(r,r.chiRms,'Cut-invariant error $\chi_{\rm tot}^{\rm rms}$',s);
    save_plot(f,out,'Figure5_total_error_vs_N');
end
if ismember(6,figureNumbers)
    % Depth error: each BSE closure against the exact tail at site N.
    f=errors_vs_N(r,r.chiRms(:,1:2),'Depth error $\chi_N^{\rm rms}$',s,1:2);
    save_plot(f,out,'Figure6_depth_error_vs_N');
end
if ismember(7,figureNumbers)
    % Closure error: each MC closure against its own target at site N.
    f=errors_vs_N(r,r.closureRms,'Closure error $\chi_M^{\rm rms}$',s,3:5);
    save_plot(f,out,'Figure7_closure_error_vs_N');
end
% Figure 8: targets use the BSE styles and the exact tail is black; MC closures use their own styles.
targetColors=[s.colors(1:2,:);0 0 0];targetStyles=[s.styles(1:2) {'-'}];targetWidths=[s.widths(1:2) 2.4];
if ismember(8,figureNumbers)
    % Each MC closure with its own target at site N.
    q=r.closureComparison;f=make(c,'MC closures and their targets',[750 500]);ax=axes(f);hold(ax,'on');h=[];
    for k=[3 1 2]
        h=[h plot(ax,q.omega,real(q.targets(:,k)),targetStyles{k},'Color',targetColors(k,:),'LineWidth',targetWidths(k),'DisplayName',q.targetNames{k})]; %#ok<AGROW>
    end
    for k=1:3,h=[h plot_closure(ax,q.omega,real(q.networks(:,k)),s,k+2,q.networkNames{k})];end %#ok<AGROW>
    top=max(real(q.targets(:,1)));
    xlabel(ax,'$\omega$','Interpreter','latex');
    ylabel(ax,'$\mathrm{Re}\,\Sigma_{N+1\to N}(0^+-i\omega)$','Interpreter','latex');
    xlim(ax,[c.bath.omegaL c.bath.omegaR]);ylim(ax,[0 1.15*top]);finish(ax);
    place_legend(ax,h,'north',2,q.omega,real([q.targets q.networks]),panel_legend_width,sprintf('$M=%d$, $N=%d$',q.M,q.N));
    % Inset: the upper band edge, where the closures differ most.
    nearEdge=q.omega>=1-.06*(c.bath.omegaR-c.bath.omegaL);
    ins=axes('Parent',f,'Position',[.39 .17 .26 .24]);hold(ins,'on');
    for k=[3 1 2],plot(ins,q.omega(nearEdge),real(q.targets(nearEdge,k)),targetStyles{k},'Color',targetColors(k,:),'LineWidth',targetWidths(k));end
    for k=1:3,plot_closure(ins,q.omega(nearEdge),real(q.networks(nearEdge,k)),s,k+2);end
    xlim(ins,[min(q.omega(nearEdge)) c.bath.omegaR]);set(ins,'FontSize',10);box(ins,'on');
    save_plot(f,out,'Figure8_self_energy_at_cut');
end
if ismember(9,figureNumbers)
    % Total error of every closure against the exact tail, the reference. For the BSE
    % closures it equals the depth error; for the chi fit to the exact tail, the closure error.
    q=r.closureComparison;f=make(c,'Total cut-invariant error',[750 500]);ax=axes(f);hold(ax,'on');set(ax,'YScale','log');
    values=[q.depth(:,1:2),q.total(:,1:3)];h=[];
    for k=[1 2 0 3 4 5]
        if k==0
            % Empty legend slot, so the MC closures fill the second column as in Figure 8.
            h=[h plot(ax,NaN,NaN,'LineStyle','none','Marker','none','DisplayName',' ')]; %#ok<AGROW>
        else
            h=[h plot_closure(ax,q.omega,values(:,k),s,k)]; %#ok<AGROW>
        end
    end
    inside=abs(q.omega-r.fit.omegaInf)<=q.bandFraction*2*r.fit.kappaInf;
    low=10^floor(log10(min(min(values(inside,:)))));
    xlabel(ax,'$\omega$','Interpreter','latex');ylabel(ax,'Cut-invariant error $\chi_{\rm tot}(\omega)$','Interpreter','latex');
    xlim(ax,[c.bath.omegaL c.bath.omegaR]);ylim(ax,[low 10^(.3+log10(max(values(:))))]);finish(ax);
    place_legend(ax,h,'north',2,q.omega,values,panel_legend_width,sprintf('$M=%d$, $N=%d$',q.M,q.N));
    save_plot(f,out,'Figure9_total_error_at_cut');
end
end

function s=closure_styles(names)
% One style per closure, shared by all figures.
s.names=names;
s.colors=[.25 .25 .25;.466 .674 .188;.494 .184 .556;.85 .325 .098;0 .447 .741];
s.styles={'--',':','-','-.','--'};s.widths=[1.5 2 1.8 1.5 1.3];
end

function h=plot_closure(ax,x,y,s,k,name)
if nargin<6,name=s.names{k};end
h=plot(ax,x,y,s.styles{k},'Color',s.colors(k,:),'LineWidth',s.widths(k),'DisplayName',name);
end

function f=errors_vs_N(r,values,label,s,closures)
% Error versus the number of explicit chain modes, with the legend above the curves.
% Column j of values belongs to closure closures(j). Lines only, with the closure styles
% of Figures 8 and 9. With all five closures, the legend has the BSE closures in the
% first column and the MC closures in the second.
if nargin<5,closures=1:5;end
c=r.cfg;N=c.Nscan(:);columns=1;order=1:numel(closures);
if isequal(closures,1:5),columns=2;order=[1 2 0 3 4 5];end
f=make(c,'Error versus N',[720 530]);ax=axes(f);hold(ax,'on');set(ax,'YScale','log');h=[];
for j=order
    if j==0
        % Empty legend slot, so the MC closures fill the second column.
        h=[h plot(ax,NaN,NaN,'LineStyle','none','Marker','none','DisplayName',' ')]; %#ok<AGROW>
    else
        h=[h plot_closure(ax,N,values(:,j),s,closures(j))]; %#ok<AGROW>
    end
end
xlabel(ax,'Explicit chain modes $N$','Interpreter','latex');ylabel(ax,label,'Interpreter','latex');
xlim(ax,[N(1) N(end)]);
low=10^floor(log10(min(values(:))/1.4));top=max(values(:));
ylim(ax,[low 10^(log10(top)+.2)]);finish(ax);
place_legend(ax,h,'north',columns,N,values,panel_legend_width);
end

function lg=place_legend(ax,h,location,columns,x,Y,maxWidth,legendTitle)
% Legend with font size 17. The figures are reduced to about half their width in the
% manuscript, so this is close to the caption size there. The font is reduced until the
% legend is at most maxWidth times the axes wide. The upper y-limit is then raised
% until the curves in the horizontal range of the legend lie below it.
lg=legend(ax,h,'Location',location,'NumColumns',columns,'Interpreter','latex','FontSize',17);
if nargin>=8,title(lg,legendTitle,'Interpreter','latex');end
drawnow;
while width_of(lg)>maxWidth*width_of(ax) && get(lg,'FontSize')>11
    set(lg,'FontSize',get(lg,'FontSize')-1);
    if nargin>=8 && isprop(lg,'Title'),set(get(lg,'Title'),'FontSize',get(lg,'FontSize'));end
    drawnow;
end
% Two passes, since new tick labels can move the axes and the legend with them.
for pass=1:2,clear_legend(ax,lg,x,Y);end
end

function w=panel_legend_width
% Largest legend width of the panel figures, as a fraction of the axes width. A centered
% legend then leaves 12.5% of the axes width free on the left for the panel label, which
% reaches about 9% in Figure 8.
w=.75;
end

function w=width_of(obj)
p=get(obj,'Position');w=p(3);
end

function clear_legend(ax,lg,x,Y)
% Raise the upper y-limit until the curves in the horizontal range of the legend lie below
% its lower edge, with small margins. The y-limit is never lowered.
drawnow;a=get(ax,'Position');p=get(lg,'Position');xl=xlim(ax);yl=ylim(ax);
xr=xl(1)+([p(1)-.03*a(3) p(1)+p(3)+.03*a(3)]-a(1))/a(3)*diff(xl);
bottom=(p(2)-a(2))/a(4)-.03;
v=Y(x>=xr(1) & x<=xr(2),:);v=v(isfinite(v));
isLog=strcmp(get(ax,'YScale'),'log');
if isLog,v=v(v>0);end
if isempty(v) || bottom<=.1,return;end
if isLog
    need=10^(log10(yl(1))+(log10(max(v))-log10(yl(1)))/bottom);
else
    need=yl(1)+(max(v)-yl(1))/bottom;
end
if need>yl(2),ylim(ax,[yl(1) need]);end
end

function f=make(c,name,sz)
f=figure('Color','w','Name',name,'Visible',c.figureVisible,'Position',[100 100 sz]);
set(f,'DefaultAxesFontSize',14,'DefaultTextFontSize',14,'DefaultAxesLineWidth',.8);
end

function finish(ax)
box(ax,'on');ax.XMinorTick='on';ax.YMinorTick='on';
end

function save_plot(f,out,name)
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector');
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',300);savefig(f,fullfile(out,[name '.fig']));
if strcmp(f.Visible,'off'),close(f);end
end
