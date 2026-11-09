function [q,info]=lm_least_squares(fun,q,opts)
%LM_LEAST_SQUARES Levenberg-Marquardt minimization of sum(r.^2), base MATLAB.
% [r,J]=fun(q) returns the real residual vector and its Jacobian.
% Damping update of Nielsen (1999) with More's monotone diagonal scaling.
if nargin<3,opts=struct();end
maxIter=get_opt(opts,'maxIter',500);tol=get_opt(opts,'tol',1e-13);
q=q(:);[r,J]=fun(q);cost=r.'*r;
info=struct('initialCost',cost,'iterations',0,'evaluations',1);
A=J.'*J;g=J.'*r;D=max(diag(A),eps);
lambda=get_opt(opts,'lambda',1e-3)*max(diag(A)./D);nu=2;stall=0;
for it=1:maxIter
    D=max(D,diag(A));
    step=-(A+lambda*diag(D))\g;
    predicted=-g.'*step+lambda*(step.'*(D.*step));
    qn=q+step;[rn,Jn]=fun(qn);info.evaluations=info.evaluations+1;costn=rn.'*rn;
    rho=(cost-costn)/max(predicted,realmin);
    if isfinite(costn) && rho>0
        gain=(cost-costn)/cost;q=qn;r=rn;J=Jn;cost=costn;A=J.'*J;g=J.'*r;
        lambda=lambda*max(1/3,1-(2*rho-1)^3);nu=2;
        if gain<tol,stall=stall+1;else,stall=0;end
        if stall>=10 || norm(step)<tol*(1+norm(q)),break;end
    else
        lambda=lambda*nu;nu=2*nu;
        if lambda>1e20,break;end
    end
    info.iterations=it;
end
info.iterations=it;info.cost=cost;info.lambda=lambda;
end

function value=get_opt(opts,name,default)
if isfield(opts,name),value=opts.(name);else,value=default;end
end
