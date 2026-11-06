function mc=mc_model_from_network(net,omegaRef,kappaRef,portLink,source)
%MC_MODEL_FROM_NETWORK Physical closure modes from a network in reference units.
% Mode frequencies omegaRef+kappaRef*h, hoppings kappaRef*t, Lindblad rates
% kappaRef*gam, and couplings portLink*c to the attachment site. The fields
% match mc_parameters, so mc_environment_matrix can use the result directly.
if nargin<5,source='network';end
assert(kappaRef>0 && portLink>0);
M=numel(net.h);
mc.M=M;
mc.omega=omegaRef+kappaRef*net.h(:);
mc.gamma=kappaRef*net.gam(:);
mc.hopping=kappaRef*net.t(:);
mc.c=portLink*net.c(:);
mc.h=diag(mc.omega)+diag(mc.hopping,1)+diag(mc.hopping,-1);
mc.portLink=portLink;
mc.omegaRef=omegaRef;mc.kappaRef=kappaRef;mc.network=net;
assert(all(mc.gamma>=0));
mc.source=source;
end
