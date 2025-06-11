function [y] = proxCTV_x(x,prox_opts)
% Proximal update of Constrained TV regulation with FGP method
% Beck, A. & Teboulle, M. Fast Gradient-Based Algorithms for Constrained Total Variation Image Denoising and Deblurring Problems. IEEE Transactions on Image Processing 18, 2419-2434, doi:10.1109/tip.2009.2028250 (2009).
gamma=prox_opts.gamma;
lambda_TV_x=prox_opts.lambda_TV_x;
prox_iter_num=prox_opts.prox_iter_num;


[n1,n2] = size(x);
pq = zeros(n1,n2);
rs = zeros(n1,n2);
t=1;

for i=1:prox_iter_num
    t_prev=t;
    w_prev=pq;
    
    w_temp=1/8/gamma*grad(Proj_C(x-gamma*div(rs,lambda_TV_x)),lambda_TV_x);
    pq=Proj_P(rs+w_temp);
    
    t=(1+sqrt(1+4*t^2))/2;
    
    rs=pq+(t_prev-1)/t*(pq-w_prev);
end
y=Proj_C(x-gamma*div(pq,lambda_TV_x));
end


function [g] = grad(u,lambda_x)
%gradient
    g = cat(1,-diff(u, 1, 1), zeros(size(u(end,:))) )*lambda_x;
end

function [d] = div(g,lambda_x)
%gradient adjoint
    d=cat(1, zeros(size(-g(end,:))),diff(g, 1, 1))*lambda_x;
end



function [y] = Proj_C(x)
    y=max(x,0);
end

function [y] = Proj_P(x)
    y=x;
end
% 
% function [y] = Proj_C(x)
%     y=min(max(x,0),1);
% end
% 
% function [y] = Proj_P(x)
%     y=min(max(x,-1),1);
% end