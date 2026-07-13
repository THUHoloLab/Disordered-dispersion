function x = normalize_nonnegative(x)
%NORMALIZE_NONNEGATIVE Clip negative/NaN values and normalize by max.

x(x < 0) = 0;
x(isnan(x)) = 0;
max_value = max(x(:));
if max_value > 0
    x = x ./ max_value;
end
end
