function [x_record, res_record, x] = FISTA_3D(b, A, A_adj, opts)
%FISTA_3D Reconstruct a hyperspectral cube with mixed_l1l2 only.
%
% The GitHub demo intentionally keeps only the reconstruction method used in
% the manuscript demos:
%   mixed_l1l2 = L1-TV regularization in the two angular dimensions plus L2
%   smoothing along wavelength, with nonnegative projection.

denoise_method = get_option_field(opts, 'denoise_method', 'mixed_l1l2');
if ~strcmp(denoise_method, 'mixed_l1l2')
    error('FISTA_3D in this demo only supports opts.denoise_method = ''mixed_l1l2''.');
end

iter_num = opts.iter_num;
record_index = opts.record_index;
display_every = opts.display_every;

lambda_tv_xy = get_option_field(opts, 'lambda_tv_xy', get_option_field(opts, 'lambda_l1', 5e-4));
lambda_tv_w = get_option_field(opts, 'lambda_tv_w', get_option_field(opts, 'lambda_l2', 1e0));
Lip = get_option_field(opts, 'lip', 2);

prox_opts.gamma = 1 / Lip;
prox_opts.lambda_TV_xy = lambda_tv_xy;
prox_opts.lambda_TV_w = lambda_tv_w;
prox_opts.prox_iter_num = get_option_field(opts, 'prox_iter_num', 10);
prox = @(z) proxCTV_xy(z, prox_opts);
dF = @(z) A_adj(A(z) - b) + 0.5 * lambda_tv_w * div_w(grad_w(z));

figure();
initial_value = A_adj(b);
x = zeros(size(initial_value), 'like', initial_value);
x_record = zeros([length(record_index), size(x)], 'like', x);
res_record = zeros(1, iter_num, 'like', b);
record_num = 1;

y = x;
t = 1;

for iter_idx = 1:iter_num
    x_prev = x;
    t_prev = t;

    x = prox(y - dF(y) / Lip);
    t = (1 + sqrt(1 + 4 * t^2)) / 2;
    y = x + (t_prev - 1) / t * (x - x_prev);

    res_record(iter_idx) = norm(A(x) - b);
    fprintf('iter= %d | residual= %g\n', iter_idx, gather(res_record(iter_idx)));

    if record_num <= length(record_index) && iter_idx == record_index(record_num)
        x_record(record_num, :, :, :) = x;
        record_num = record_num + 1;
    end

    if isfinite(display_every) && display_every > 0 && mod(iter_idx, display_every) == 0
        plot_tensor_and_loss(gather(x), gather(res_record), 3, 3);
    end
end

end

function g = grad_w(u)
% Forward wavelength gradient with a zero Neumann boundary at the last band.
g = cat(3, -diff(u, 1, 3), zeros(size(u(:, :, end)), 'like', u));
end

function d = div_w(g)
% Adjoint of grad_w. The first wavelength band is handled explicitly.
d = cat(3, g(:, :, 1), diff(g, 1, 3));
end

function plot_tensor_and_loss(x, res_record, num_slices, num_plots)
% Display angular-wavelength slices, selected spectra, and the loss curve.
clf;
[dim_x, dim_y, dim_w] = size(x);

x_indices = round(linspace(round(dim_x / 6), round(dim_x / 6 * 5), num_slices));
y_indices = round(linspace(round(dim_y / 6), round(dim_y / 6 * 5), num_slices));
w_indices = round(linspace(round(dim_w / 6), round(dim_w / 6 * 5), num_slices));

for idx = 1:num_slices
    subplot(3, 5, 1 + (idx - 1) * 5);
    imagesc(squeeze(x(x_indices(idx), :, :)));
    title(['X slice ', num2str(idx)]);
    colorbar;

    subplot(3, 5, 2 + (idx - 1) * 5);
    imagesc(squeeze(x(:, y_indices(idx), :)));
    title(['Y slice ', num2str(idx)]);
    colorbar;

    subplot(3, 5, 3 + (idx - 1) * 5);
    imagesc(squeeze(x(:, :, w_indices(idx))));
    title(['W slice ', num2str(idx)]);
    colorbar;
end

subplot(3, 5, [4, 9, 14]);
hold on;
x_sample_indices = round(linspace(round(dim_x / 6), round(dim_x / 6 * 5), num_plots));
y_sample_indices = round(linspace(round(dim_y / 6), round(dim_y / 6 * 5), num_plots));
[xx_sample_indices, yy_sample_indices] = meshgrid(x_sample_indices, y_sample_indices);
for idx = 1:num_plots^2
    plot(squeeze(x(xx_sample_indices(idx), yy_sample_indices(idx), :)));
end
title('Spectral profiles');
hold off;

subplot(3, 5, [5, 10, 15]);
semilogy(res_record);
title('Loss');
drawnow;
end
