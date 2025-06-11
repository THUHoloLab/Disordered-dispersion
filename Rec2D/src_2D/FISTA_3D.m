function [x_record,res_record,x]=FISTA_3D(b,A,A_adj,opts)
denoise_method=opts.denoise_method;
iter_num=opts.iter_num;
record_index=opts.record_index;
display_every=opts.display_every;

Lip=2;  % A Lipschitz constant of ∇f.

switch denoise_method

    case 'mixed_l1l2'
        prox_opts.gamma=1/Lip;
        prox_opts.lambda_TV_xy=5e-4; % l1 norm
        prox_opts.lambda_TV_w=1e0; % l2 norm
        prox_opts.prox_iter_num=1;
        prox = @(x) proxCTV_xy(x,prox_opts);

        dF=@(x)  A_adj(A(x)-b) + 0.5*prox_opts.lambda_TV_w*div_w(grad_w(x));
end

figure()
%% FISTA MAIN
x=zeros(size(A_adj(b)));
x_record=zeros([length(record_index),size(x)]);record_num=1;
res_record=zeros(1,iter_num);
t=1;
y=x;
res_record=zeros(1,iter_num);

for i=1:iter_num

    x_prev=x;
    t_prev=t;
    
    x=prox(y-dF(y)/Lip);
    t=(1+sqrt(1+4*t^2))/2;
    y=x+(t_prev-1)/t*(x-x_prev);
    
    res_record(i)=norm(A(x)-b); %Fidelity error
    fprintf(['iter= ' ,num2str(i), ' | residual= ' ,num2str(res_record(i)), '\n'])
    
    if record_num<=length(record_index)
        if i==record_index(record_num)
            x_record(record_num,:,:,:)=x;
            record_num=record_num+1;
        end
    end
    
    if mod(i,display_every) == 0 
        plot_tensor_and_loss(gather(x), res_record, 3, 3)
    end

end



end




function [g] = grad_w(u)
%gradient
    g = cat(3,-diff(u, 1, 3), zeros(size(u(:,:,end))) );
end
function [d] = div_w(g)
%gradient adjoint
    d=cat(3, zeros(size(g(:,:,1))),  diff(g, 1, 3));
end






function plot_tensor_and_loss(x, res_record, num_slices, num_plots)
    clf
    % 获取张量的尺寸
    [dim_x, dim_y, dim_z] = size(x);
    num_slices=3;
    % 在 x 方向均匀截取 num_slices 个 slice
    x_indices = round(linspace(round(dim_x/6), round(dim_x/6*5), num_slices));
    for i = 1:num_slices
        subplot(3, 5, 1+(i-1)*5);
        x_slice = squeeze(x(x_indices(i), :, :));
        imagesc(x_slice);
        title(['X Slice ', num2str(i)]);
        colorbar
    end

    % 在 y 方向均匀截取 num_slices 个 slice
    y_indices = round(linspace(round(dim_y/6), round(dim_y/6*5), num_slices));
    for i = 1:num_slices
        subplot(3, 5, 2+(i-1)*5);
        y_slice = squeeze(x(:, y_indices(i), :));
        imagesc(y_slice);
        title(['Y Slice ', num2str(i)]);
        colorbar
    end

    % 在 z 方向均匀截取 num_slices 个 slice
    z_indices = round(linspace(round(dim_z/6), round(dim_z/6*5), num_slices));
    for i = 1:num_slices
        subplot(3, 5, 3+(i-1)*5);
        z_slice = squeeze(x(:, :, z_indices(i)));
        imagesc(z_slice);
        title(['W Slice ', num2str(i)]);
        colorbar
    end

    % 第四列，画多个 (x, y) 点处的 z 轴一维 plot
    subplot(3, 5, [4,9,14]);
    hold on;
    x_sample_indices = round(linspace(round(dim_x/6), round(dim_x/6*5), num_plots));
    y_sample_indices = round(linspace(round(dim_y/6), round(dim_y/6*5), num_plots));

    [xx_sample_indices,yy_sample_indices]=meshgrid(x_sample_indices,y_sample_indices);

    for i = 1:num_plots^2
        plot(squeeze(x(xx_sample_indices(i), yy_sample_indices(i), :)));
    end
    title('W-axis Plots');
    hold off;

    % 第五列，画 loss 曲线
    subplot(3, 5, [5,10,15]);
    semilogy(res_record);
    title('Loss');

    % 确保图形更新
    drawnow;
end
