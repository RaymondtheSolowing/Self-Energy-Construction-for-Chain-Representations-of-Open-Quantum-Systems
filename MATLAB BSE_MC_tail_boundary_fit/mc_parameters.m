function mc=mc_parameters(M,omegaInf,kappaInf,portLink)
%MC_PARAMETERS Published universal MC6/8/10 coefficients, without refitting.
% H_int = sum_j (c_j a_N^dagger d_j + conj(c_j) d_j^dagger a_N).
% The normalized coupling vector is scaled by portLink, the coupling to the
% attachment site N (mc_published_model chooses kappa_N or kappaInf).
% Internal frequencies, hoppings and losses use omegaInf and kappaInf.
% No normalization or rounding of the published data.
if nargin<4, portLink=kappaInf; end
assert(ismember(M,[6 8 10]) && kappaInf>0 && portLink>0);
folder=fullfile(fileparts(mfilename('fullpath')),'mc_coefficients');
a=load(fullfile(folder,sprintf('alphas_%d.dat',M)));
b=load(fullfile(folder,sprintf('betas_%d.dat',M)));
c=load(fullfile(folder,sprintf('coupls_%d.dat',M)));
assert(isequal(size(a),[M 2]) && isequal(size(b),[M-1 2]) && isequal(size(c),[M 2]));
assert(all(b(:,1)==0));
mc.M=M;
mc.omega=omegaInf-2*kappaInf*a(:,2);
mc.gamma=-4*kappaInf*a(:,1);
mc.hopping=-2*kappaInf*b(:,2);
mc.c=portLink*(c(:,1)+1i*c(:,2));
mc.h=diag(mc.omega)+diag(mc.hopping,1)+diag(mc.hopping,-1);
mc.portLink=portLink;
assert(all(mc.gamma>=0));
mc.source='phaerrax/markovian_closure_fermions@65d476f0eda2b6c36ecd0f8d3785b481d2526ab1';
end
