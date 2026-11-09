function [Eb,Zb,info] = continuum_bound_states(epsilon,bath,Emax)
%CONTINUUM_BOUND_STATES Bound poles of the exact continuum qubit Green function.
%
% Outside the band the pole equation is
%   E - epsilon - (1/pi) int J(nu)/(E-nu) dnu = 0.
% If J(edge) > 0 the integral diverges logarithmically at that edge, so a
% pole can lie exponentially close to the band.  The search therefore uses
% the edge distance d = |E-edge| on a logarithmic grid, and the edge value
% J(edge) is subtracted analytically from both integrals.

if nargin < 3 || isempty(Emax)
    Emax = 10;
end

L = bath.omegaL;
R = bath.omegaR;
span = R-L;
[JL,dJL] = structured_spectral_density(L,bath);
[JR,dJR] = structured_spectral_density(R,bath);
q = {'RelTol',1e-12,'AbsTol',1e-16};

Eb = [];
Zb = [];
info.edgeDistance = [];

for side = [-1 1]
    if side > 0
        edge = R; Je = JR; dJe = dJR; s = @(nu) R-nu;
    else
        edge = L; Je = JL; dJe = dJL; s = @(nu) nu-L;
    end
    lin = @(nu) dJe.*(nu-edge);

    % S(E) = (1/pi) int J/(E-nu), with E = edge + side*d
    Sfun = @(d) side*( integral(@(nu) ...
        (structured_spectral_density(nu,bath)-Je)./(d+s(nu)),L,R,q{:}) ...
        + Je*log((d+span)/d) )/pi;
    % T(E) = (1/pi) int J/(E-nu)^2.  Near the edge the analytic terms
    % dominate, so the absolute tolerance of the remainder is scaled by them.
    Tsing = @(d) Je*(1/d-1/(d+span)) ...
        - side*dJe*(log((d+span)/d)+d/(d+span)-1);
    Tfun = @(d) ( integral(@(nu) ...
        (structured_spectral_density(nu,bath)-Je-lin(nu))./(d+s(nu)).^2,L,R, ...
        'RelTol',1e-12,'AbsTol',1e-14*max(1,abs(Tsing(d)))) + Tsing(d) )/pi;

    Ffun = @(u) (edge+side*exp(u)) - epsilon - Sfun(exp(u));
    uGrid = linspace(log(1e-300*span),log(Emax*span),400);
    Fg = arrayfun(Ffun,uGrid);
    k = find(isfinite(Fg(1:end-1)) & isfinite(Fg(2:end)) & ...
        sign(Fg(1:end-1)).*sign(Fg(2:end)) <= 0,1,'first');
    if isempty(k)
        continue;
    end
    d = exp(fzero(Ffun,uGrid(k:k+1),optimset('TolX',1e-13,'Display','off')));
    Z = 1/(1+Tfun(d));
    if Z > 0 && Z <= 1+1e-8
        Eb(end+1) = edge+side*d; %#ok<AGROW>
        Zb(end+1) = real(Z); %#ok<AGROW>
        info.edgeDistance(end+1) = d; %#ok<AGROW>
    end
end

[Eb,order] = sort(Eb);
Zb = Zb(order);
if ~isempty(info.edgeDistance)
    info.edgeDistance = info.edgeDistance(order);
end
end
