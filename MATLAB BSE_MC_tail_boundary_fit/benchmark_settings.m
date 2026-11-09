function cfg=benchmark_settings(overrides)
%BENCHMARK_SETTINGS Settings of the five-closure benchmark.
% Closures of a chain with N explicit modes: asymptotic BSE, chain-extended
% BSE, semicircle-fit (published) MC, MC chi fit to the chain-extended BSE
% closure, and MC chi fit to the exact tail. Every MC closure adds M damped
% modes coupled to site N.
if nargin<1,overrides=struct();end
cfg.protocolVersion='BSE-MC-closures-cut-invariant-error-2026-11-02-v4';
cfg.bath=struct('omegaL',0,'omegaR',1,'alpha',.02,'bumpAmplitude',4,'bumpCenter',.70,'bumpWidth',.07);
cfg.epsilon=.5;cfg.tMax=2000;cfg.etaSpec=[];cfg.dt=.1;
% Nscan: chain modes for Figures 3-7; Nmain: Figures 1-2; Ncut: Figures 8-9;
% Ntable: Table I and the stored curves. M: damped modes of every MC closure.
cfg.tFigure1Max=100;cfg.Nmain=15;cfg.Nscan=10:60;cfg.Ncut=50;cfg.Ntable=[15 30 45 60];cfg.M=10;
% Couplings of the MC modes to site N: 'kappaN' scales the published couplings
% by kappa_N, the attachment coupling of the retained chain; 'kappaInf' uses
% kappaInf, the rescaling of the published MC. Only the semicircle-fit MC
% depends on this choice; the chi fits absorb the scale.
cfg.mcPortLink='kappaN';
% Cut-invariant error on the real axis, band interior |w-omegaInf|<=chiBandFraction*2*kappaInf,
% here 0.01<=w<=0.99; the chi fits use the same window. chi=1 at the band edges for any
% closure. The relative spectral error uses the full band 0<=w<=1.
cfg.chiBandFraction=.98;cfg.chiGridN=2001;cfg.chiPVN=256;
% Chi fits of MC-form networks: in-band chi plus a weak absolute penalty for
% outOfBandRange(1)<=|w-omegaRef|/kappaRef<=outOfBandRange(2). The universal
% fit to the semicircle is done once; the fits at each N start from it.
cfg.chiFit=struct('outOfBandWeight',.3,'outOfBandRange',[2.05 8],'gridIn',1201,'gridOut',150, ...
    'maxIterUniversal',20000,'maxIterLocal',2000,'tol',1e-13,'sgMin',1e-2);
cfg.reuseChiFitCache=true;
cfg.frequencyGridN=4001;cfg.chainQuadN=4000;cfg.continuumSpecQuadN=1200;cfg.spectralGridN=8001;
cfg.compositeOpts=struct('Nquad',4800,'Npv',256,'edgeOrder',24);
cfg.bseOpts=struct('eta',0,'normalizationTolerance',1e-10,'phaseStepMax',pi/8, ...
    'NkMin',8001,'findBoundStates',true,'boundEmax',10,'maxRefinements',3);
cfg.outputDir=fullfile(fileparts(mfilename('fullpath')),'results');
cfg.figureVisible='off';
names=fieldnames(overrides);
for k=1:numel(names)
    assert(isfield(cfg,names{k}),'Unknown setting %s',names{k});cfg.(names{k})=overrides.(names{k});
end
if isempty(cfg.etaSpec),cfg.etaSpec=1/cfg.tMax;end
assert(cfg.etaSpec>0 && cfg.tMax>0 && cfg.dt>0 && cfg.tFigure1Max<=cfg.tMax);
N=cfg.Nscan;
assert(isequal(N,unique(round(N))) && min(N)>=1,'Nscan must be increasing positive integers.');
assert(ismember(cfg.Nmain,N) && ismember(cfg.Ncut,N) && all(ismember(cfg.Ntable,N)));
assert(ismember(cfg.M,[6 8 10]),'Published MC coefficients exist for M=6, 8, 10.');
assert(ismember(cfg.mcPortLink,{'kappaN','kappaInf'}));
assert(cfg.chiBandFraction>0 && cfg.chiBandFraction<1 && cfg.spectralGridN>=2001);
end
