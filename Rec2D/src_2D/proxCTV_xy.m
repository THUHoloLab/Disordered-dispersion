function [x] = proxCTV_xy(b,prox_opts)
% Proximal update of Constrained TV regulation with FGP method
% Beck, A. & Teboulle, M. Fast Gradient-Based Algorithms for Constrained Total Variation Image Denoising and Deblurring Problems. IEEE Transactions on Image Processing 18, 2419-2434, doi:10.1109/tip.2009.2028250 (2009).
gamma=prox_opts.gamma;
lambda_TV_xy=prox_opts.lambda_TV_xy;
prox_iter_num=prox_opts.prox_iter_num;


[n1,n2,n3] = size(b); %n3 is wavelength channel
pq = zeros(n1,n2,n3,2);
rs = zeros(n1,n2,n3,2);
t=1;

for i=1:prox_iter_num
    t_prev=t;
    w_prev=pq;
    
    w_temp=1/12/gamma*grad(Proj_C(b-gamma*div(rs,lambda_TV_xy)),lambda_TV_xy);
    pq=Proj_P(rs+w_temp);
    
    t=(1+sqrt(1+4*t^2))/2;
    
    rs=pq+(t_prev-1)/t*(pq-w_prev);
end
x=Proj_C(b-gamma*div(pq,lambda_TV_xy));
end


function [g] = grad(u,lambda_xy)
[n1,n2,n3] = size(u);
%gradient
    g(:, :, :, 1) = cat(1,   -diff(u, 1, 1), zeros(1,n2,n3)   )*lambda_xy;
    g(:, :, :, 2) = cat(2,   -diff(u, 1, 2), zeros(n1,1,n3)   )*lambda_xy;
end

function [d] = div(g,lambda_xy)
[n1,n2,n3,~] = size(g);
%gradient adjoint
    d1=cat(1,    zeros(1,n2,n3) , diff(g(:, :, :, 1), 1, 1)    )*lambda_xy;
    d2=cat(2,    zeros(n1,1,n3) , diff(g(:, :, :, 2), 1, 2)    )*lambda_xy;
    d=d1+d2;
end


function [y] = Proj_C(x)
    y=max(x,0);
end
function [y] = Proj_P(x)
    y=x;
end

% function [y] = Proj_C(x)
%     y=min(max(x,0),1);
% end
% 
% function [y] = Proj_P(x)
%     y=min(max(x,-1),1);
% end