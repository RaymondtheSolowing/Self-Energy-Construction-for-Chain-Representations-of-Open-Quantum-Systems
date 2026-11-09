function q=chain_only_bse_parameters(N,omega,kappa,omegaInf,kappaInf)
% Least-squares 1/n^2 continuation from last third of AVAILABLE coefficients.
% Inputs needed: omega(1:N+1), kappa(1:N), and band-edge asymptotes only.
% omega(n+1)-omegaInf ~= a/n^2; kappa(n)-kappaInf ~= b/n^2.
assert(numel(omega)>=N+1 && numel(kappa)>=N);
% Take the last ceil(length/3) values of EACH available coefficient series.
% Exclude omega_1 because its corresponding n=0 cannot enter 1/n^2.
nw=(max(1,N+1-ceil((N+1)/3)):N)';nk=(N-ceil(N/3)+1:N)';
xw=1./nw.^2;xk=1./nk.^2;
a=xw\(omega(nw+1)-omegaInf);b=xk\(kappa(nk)-kappaInf);
q=struct('N',N,'a',a,'b',b,'omegaSampleIndices',nw+1,'kappaSampleIndices',nk,'sampleCount',min(numel(nw),numel(nk)), ...
    'omegaBar',omegaInf+a/(N+1)^2,'kappaBar',kappaInf+b/(N+1)^2, ...
    'method','Last-third coefficient fit to 1/n^2');
q.bandLower=q.omegaBar-2*q.kappaBar;q.bandUpper=q.omegaBar+2*q.kappaBar;
% Report validity; never silently clip the proposed rule to improve its error.
pivot=omega(1);
for j=2:N+1, pivot=omega(j)-kappa(j-1)^2/pivot; end
q.prefixZeroEnergyPivot=pivot;
q.nonnegativeBath=q.kappaBar>0 && q.omegaBar>=2*q.kappaBar+max(q.kappaBar-pivot,0)^2/pivot;
end
