function net=mc_published_network(M)
%MC_PUBLISHED_NETWORK Published universal MC6/8/10 in reference units.
% Same coefficients as mc_parameters, in units of kappaInf with frequencies
% measured from omegaInf: h=-2*alpha_I, t=-2*beta_I, gam=-4*alpha_R, c=coupls.
assert(ismember(M,[6 8 10]));
folder=fullfile(fileparts(mfilename('fullpath')),'mc_coefficients');
a=load(fullfile(folder,sprintf('alphas_%d.dat',M)));
b=load(fullfile(folder,sprintf('betas_%d.dat',M)));
c=load(fullfile(folder,sprintf('coupls_%d.dat',M)));
net=struct('h',-2*a(:,2),'t',-2*b(:,2),'gam',-4*a(:,1),'c',c(:,1)+1i*c(:,2));
assert(all(net.gam>=0));
end
