function [x_record, res_record, x] = FISTA_2D(b, A, A_adj, opts)
%FISTA_2D Reconstruct a spectral-angular 2D object with mixed_l1l2 only.
%
% The GitHub demo intentionally keeps only the reconstruction method used in
% the manuscript demos:
%   mixed_l1l2 = L1-TV regularization along theta_x plus L2 smoothing along
%   wavelength, with nonnegative projection.

denoise_method = get_option_field(opts, 'denoise_method', 'mixed_l1l2');
if ~strcmp(denoise_method, 'mixed_l1l2')
    error('FISTA_2D in this demo only supports opts.denoise_method = ''mixed_l1l2''.');
end

iter_num = opts.iter_num;
record_index = opts.record_index;
display_every = opts.display_every;

lambda_tv_x = get_option_field(opts, 'lambda_tv_x', get_option_field(opts, 'lambda_l1', 1e-3));
lambda_tv_w = get_option_field(opts, 'lambda_tv_w', get_option_field(opts, 'lambda_l2', 1e0));

Lip = get_option_field(opts, 'lip', max(5, 5 + 2 * lambda_tv_w));

prox_opts.gamma = 1 / Lip;
prox_opts.lambda_TV_x = lambda_tv_x;
prox_opts.lambda_TV_w = lambda_tv_w;
prox_opts.prox_iter_num = get_option_field(opts, 'prox_iter_num', 5);
prox = @(z) proxCTV_x(z, prox_opts);
dF = @(z) A_adj(A(z) - b) + 0.5 * lambda_tv_w * div_w(grad_w(z));

figure(1);
x = rand(size(A_adj(b)), 'like', A_adj(b));
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
        x_record(record_num, :, :) = x;
        record_num = record_num + 1;
    end

    if isfinite(display_every) && display_every > 0 && mod(iter_idx, display_every) == 0
        subplot(1, 2, 1);
        imagesc(flipud(gather(x)'));
        colorbar;
        title('Rec');

        subplot(1, 2, 2);
        semilogy(gather(res_record));
        title('Loss');
        drawnow;
    end
end

end

function g = grad_w(u)
% Forward wavelength gradient with a zero Neumann boundary at the last band.
g = cat(2, -diff(u, 1, 2), zeros(size(u(:, end)), 'like', u));
end

function d = div_w(g)
% Adjoint of grad_w. The first wavelength band is handled explicitly.
d = cat(2, g(:, 1), diff(g, 1, 2));
end
