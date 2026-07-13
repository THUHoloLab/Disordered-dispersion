function y = proxCTV_x(input_value, prox_opts)
% Proximal map for nonnegative anisotropic L1-TV along the angle axis.
%
% Solves
%   0.5 * ||y-input_value||_2^2
%   + gamma * lambda_TV_x * ||diff(y,1,1)||_1
%   + indicator(y >= 0)
% with a first-order primal-dual iteration.

gamma = prox_opts.gamma;
lambda_TV_x = prox_opts.lambda_TV_x;
prox_iter_num = prox_opts.prox_iter_num;
tv_weight = gamma * lambda_TV_x;

tau = 0.49;
sigma = 0.49;
extrapolation = 1;

y = max(input_value, 0);
y_bar = y;
dual_x = zeros(size(input_value, 1) - 1, size(input_value, 2), ...
    'like', input_value);

for iter_idx = 1:prox_iter_num
    dual_x = dual_x + sigma * diff(y_bar, 1, 1);
    dual_x = min(max(dual_x, -tv_weight), tv_weight);

    y_prev = y;
    primal_step = y - tau * angle_difference_adjoint(dual_x);
    y = max((primal_step + tau * input_value) / (1 + tau), 0);
    y_bar = y + extrapolation * (y - y_prev);
end
end

function adjoint_value = angle_difference_adjoint(dual_x)
angle_count = size(dual_x, 1) + 1;

if angle_count == 1
    adjoint_value = zeros(1, size(dual_x, 2), 'like', dual_x);
    return;
end

adjoint_value = cat(1, ...
    -dual_x(1, :), ...
    dual_x(1:end-1, :) - dual_x(2:end, :), ...
    dual_x(end, :));
end
