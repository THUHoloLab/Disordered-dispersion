function x = proxCTV_xy(input_value, prox_opts)
% Proximal map for nonnegative anisotropic L1-TV in both angle axes.
%
% Solves
%   0.5 * ||x-input_value||_2^2
%   + gamma * lambda_TV_xy * (||D_x x||_1 + ||D_y x||_1)
%   + indicator(x >= 0)
% with a first-order primal-dual iteration.

gamma = prox_opts.gamma;
lambda_TV_xy = prox_opts.lambda_TV_xy;
prox_iter_num = prox_opts.prox_iter_num;
tv_weight = gamma * lambda_TV_xy;

% tau * sigma * ||[D_x; D_y]||^2 < 1, with the norm squared bounded by 8.
tau = 0.35;
sigma = 0.35;
extrapolation = 1;

[x_count, y_count, wavelength_count] = size(input_value);
x = max(input_value, 0);
x_bar = x;
dual_x = zeros(x_count - 1, y_count, wavelength_count, ...
    'like', input_value);
dual_y = zeros(x_count, y_count - 1, wavelength_count, ...
    'like', input_value);

for iter_idx = 1:prox_iter_num
    dual_x = dual_x + sigma * diff(x_bar, 1, 1);
    dual_y = dual_y + sigma * diff(x_bar, 1, 2);
    dual_x = min(max(dual_x, -tv_weight), tv_weight);
    dual_y = min(max(dual_y, -tv_weight), tv_weight);

    x_prev = x;
    divergence = angle_difference_adjoint_x(dual_x) + ...
        angle_difference_adjoint_y(dual_y);
    primal_step = x - tau * divergence;
    x = max((primal_step + tau * input_value) / (1 + tau), 0);
    x_bar = x + extrapolation * (x - x_prev);
end
end

function adjoint_value = angle_difference_adjoint_x(dual_x)
if size(dual_x, 1) == 0
    adjoint_value = zeros(1, size(dual_x, 2), size(dual_x, 3), ...
        'like', dual_x);
    return;
end

adjoint_value = cat(1, ...
    -dual_x(1, :, :), ...
    dual_x(1:end-1, :, :) - dual_x(2:end, :, :), ...
    dual_x(end, :, :));
end

function adjoint_value = angle_difference_adjoint_y(dual_y)
if size(dual_y, 2) == 0
    adjoint_value = zeros(size(dual_y, 1), 1, size(dual_y, 3), ...
        'like', dual_y);
    return;
end

adjoint_value = cat(2, ...
    -dual_y(:, 1, :), ...
    dual_y(:, 1:end-1, :) - dual_y(:, 2:end, :), ...
    dual_y(:, end, :));
end
