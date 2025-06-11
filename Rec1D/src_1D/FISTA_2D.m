function [x_record,res_record,x]=FISTA_2D(b,A,A_adj,opts)
denoise_method=opts.denoise_method;
iter_num=opts.iter_num;
record_index=opts.record_index;
display_every=opts.display_every;

Lip=5;  % A Lipschitz constant of ∇f.

switch denoise_method
    case 'mixed_l1l2'
        prox_opts.gamma=1/Lip;
        prox_opts.lambda_TV_x=1e-3; % l1 norm
        prox_opts.lambda_TV_w=1e-0; % l2 norm
        prox_opts.prox_iter_num=5;
        prox = @(x) proxCTV_x(x,prox_opts);

        dF=@(x)  A_adj(A(x)-b) + 0.5*prox_opts.lambda_TV_w*div_w(grad_w(x));
end




figure(1)
%% FISTA MAIN
x=rand(size(A_adj(b)));
x_record=zeros([length(record_index),size(x)]);record_num=1;
res_record=zeros(1,iter_num);

y=x;
t=1;


for i=1:iter_num

x_prev=x;
t_prev=t;
x=prox(   y-1/Lip*dF(y)     );
t=(1+sqrt(1+4*t^2))/2;
y=x+(t_prev-1)/t*(x-x_prev);

res_record(i)=norm(A(x)-b);  %Fidelity error
fprintf(['iter= ' ,num2str(i), ' | residual= ' ,num2str(res_record(i)), '\n'])

if record_num<=length(record_index)
    if i==record_index(record_num)
        x_record(record_num,:,:)=x;
        record_num=record_num+1;
    end
end

if mod(i,display_every) == 0 
    subplot(1,2,1)
    x_display=gather(x);
    imagesc(flipud(x_display'));
    colorbar
    title('Rec')

    subplot(1,2,2),
    semilogy(res_record);
    title('Loss')
    drawnow
end
end





function [g] = grad_w(u)
%gradient
    g = cat(2,-diff(u, 1, 2), zeros(size(u(:,end))) );
end
function [d] = div_w(g)
%gradient adjoint
    d=cat(2,zeros(size(g(:,1))),diff(g, 1, 2));
end







end
