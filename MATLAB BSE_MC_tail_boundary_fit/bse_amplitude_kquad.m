function [A0,info] = bse_amplitude_kquad( ...
    t,N,epsilon,g,omega,kappa,omegaBar,kappaBar,omegaL,omegaR,opts)
%BSE_AMPLITUDE_KQUAD Homogeneous-continuation BSE amplitude on its tail band.
%
% For the time-domain spectral representation the physical boundary value
% on the continuum is used directly (eta = 0^+).  sigma_uniform_tail.m
% implements the retarded branch on the imaginary axis, so no artificial
% positive broadening is required here.  A finite eta would move a small
% amount of Lorentzian spectral weight outside the finite band and should
% therefore be reserved for frequency-domain plotting only.

if nargin < 11
    opts = struct();
end

span = omegaR-omegaL;
assert(max(abs([omegaL-(omegaBar-2*kappaBar),omegaR-(omegaBar+2*kappaBar)])) ...
    <1e-12*max(1,abs(omegaBar)+2*kappaBar),'Use the fitted tail band edges for BSE inversion.');
eta = get_opt(opts,'eta',0);
if eta < 0
    error('opts.eta must be nonnegative.');
end
normalizationTolerance = get_opt(opts,'normalizationTolerance',1e-8);
maxRefinements = get_opt(opts,'maxRefinements',2);
phaseStepMax = get_opt(opts,'phaseStepMax',pi/8);
NkMin = get_opt(opts,'NkMin',2001);
findBoundStates = get_opt(opts,'findBoundStates',true);
boundEmax = get_opt(opts,'boundEmax',10*max(1,span));

t = t(:).';
tMax = max(t);
vMax = 2*kappaBar;

NkPhase = ceil(pi*vMax*max(tMax,1)/phaseStepMax)+1;
Nk = ceil(max(NkMin,NkPhase));
if mod(Nk,2)==0
    Nk = Nk+1;
end

k = linspace(0,pi,Nk).';
dk = pi/(Nk-1);
omegaK = omegaBar-2*kappaBar*cos(k);
jacobian = 2*kappaBar*sin(k);

Gfun = @(p) green_bse( ...
    p,N,epsilon,g,omega,kappa,omegaBar,kappaBar);
G = Gfun(eta-1i*omegaK);
rho = real(G)/pi;

smallNegative = rho < 0 & rho > -1e-11;
rho(smallNegative) = 0;
if any(~isfinite(rho)) || any(rho < -1e-8)
    error('Invalid BSE continuum density.');
end

simpson = ones(Nk,1);
simpson(2:2:end-1) = 4;
simpson(3:2:end-2) = 2;
weight = (dk/3).*simpson.*jacobian.*rho;

Eb = [];
Zb = [];
if findBoundStates
    [Eb,Zb] = find_bound_states_from_green( ...
        Gfun,omegaL,omegaR,boundEmax,true);
end

info.Nk = Nk;
info.dk = dk;
info.phaseStepAtTmax = tMax*vMax*dk;
info.continuumWeight = sum(weight);
info.boundEnergies = Eb;
info.boundWeights = Zb;
info.totalWeight = info.continuumWeight+sum(Zb);
info.normalizationError = abs(info.totalWeight-1);
info.eta = eta;
info.refinementCount = 0;

if info.normalizationError > normalizationTolerance
    if maxRefinements > 0 && findBoundStates
        refined = opts;
        refined.NkMin = 2*(Nk-1)+1;
        refined.maxRefinements = maxRefinements-1;
        [A0,info] = bse_amplitude_kquad( ...
            t,N,epsilon,g,omega,kappa,omegaBar,kappaBar, ...
            omegaL,omegaR,refined);
        info.refinementCount = info.refinementCount+1;
        return;
    end
    error('BSE spectral weight %.12g differs from unity.',info.totalWeight);
end

Acont = zeros(size(t));
if numel(t)>1 && abs(t(1)) < 100*eps
    dt = t(2)-t(1);
    if max(abs(diff(t)-dt)) < 1e-10*max(1,abs(dt))
        phase = ones(Nk,1);
        phaseStep = exp(-1i*omegaK*dt);
        for m = 1:numel(t)
            Acont(m) = weight.'*phase;
            phase = phase.*phaseStep;
        end
    else
        Acont = nonuniform_sum(t,omegaK,weight);
    end
elseif numel(t)==1
    Acont(1) = weight.'*exp(-1i*omegaK*t(1));
else
    Acont = nonuniform_sum(t,omegaK,weight);
end

Abound = zeros(size(t));
for j = 1:numel(Eb)
    Abound = Abound+Zb(j)*exp(-1i*Eb(j)*t);
end
A0 = Acont+Abound;
end

function A = nonuniform_sum(t,omegaK,weight)
A = zeros(size(t));
chunkSize = 200;
for i1 = 1:chunkSize:numel(t)
    i2 = min(numel(t),i1+chunkSize-1);
    A(i1:i2) = weight.'*exp(-1i*(omegaK*t(i1:i2)));
end
end

function value = get_opt(opts,name,default)
if isfield(opts,name)
    value = opts.(name);
else
    value = default;
end
end
