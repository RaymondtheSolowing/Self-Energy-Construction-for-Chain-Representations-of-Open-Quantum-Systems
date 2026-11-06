function tails=exact_tail_self_energies(p,SigmaPort,g,omega,kappa,maxDepth)
%EXACT_TAIL_SELF_ENERGIES Strip an exact prefix by inverse Mobius steps.
% Row n+1 is Sigma_exact_(n+1 -> n), including n=0 (the system port).
% No uniform-tail approximation is used to construct these targets.
assert(maxDepth>=0 && maxDepth==floor(maxDepth));
assert(numel(omega)>=maxDepth && numel(kappa)>=maxDepth-1);
p=p(:).'; response=SigmaPort(:).';
assert(numel(p)==numel(response));
tails=zeros(maxDepth+1,numel(p)); tails(1,:)=response;
for n=1:maxDepth
    if n==1, link=g; else, link=kappa(n-1); end
    response=link^2./response-p-1i*omega(n);
    tails(n+1,:)=response;
end
assert(all(isfinite(tails),'all'),'Inverse Mobius transformation is nonfinite.');
assert(min(real(tails),[],'all')>=-1e-10, ...
    'Tail reference lost passivity; refine the continuum/chain calculation.');
end
