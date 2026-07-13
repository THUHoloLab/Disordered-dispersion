clear; clc; close all;

% DEMO_2D_RECONSTRUCTION
% This demo reconstructs a 2D spectral-angular object x(theta_x, lambda).
% The encoder response is loaded from an experimentally calibrated
% transmission matrix, while the object and measurement are simulated.

%% Set paths and demo parameters
% Keep all paths relative to this GitHub demo folder so the code can be
% moved to another computer after download.
demo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(demo_root, 'src', 'common'));
addpath(fullfile(demo_root, 'src', 'rec2d'));

data_dir = fullfile(demo_root, 'data');
output_dir = fullfile(demo_root, 'outputs', '2D_demo');
ensure_dir(output_dir);

calibration_path = ensure_demo_calibration_file(data_dir);
% Use the original reconstruction grid used in the manuscript code. The
% calibration matrix is not downsampled; it is interpolated onto this grid.
wavelength_range = 400:0.5:700;
theta_x_range = -60:1:60;
rng(2);

%% Load and resample the experimentally calibrated 2D response
% The calibration file contains response_xyw(theta_x, theta_y, lambda).
% For this 2D demo we use the theta_y = 0 slice, yielding a sensing matrix
% T(theta_x, lambda).
load(calibration_path, 'response_xyw', 'theta_x_data', 'wavelength_data');
xw_response_data = squeeze(response_xyw(:, 1, :));

[XX_data, WW_data] = ndgrid(theta_x_data, wavelength_data);
[XX_samples, WW_samples] = ndgrid(theta_x_range, wavelength_range);
xw_sampling_matrix = interpn(XX_data, WW_data, xw_response_data, ...
    XX_samples, WW_samples, 'spline');
xw_sampling_matrix = normalize_nonnegative(xw_sampling_matrix);

%% Construct a simulated spectral-angular ground truth object
% The object contains one compact angular feature. Every angular position in
% the object shares the same smooth spectrum generated from three random
% control points in [0.5, 1]. The random points are resampled until their
% contrast is visually clear.
xw_gt = create_demo_2d_object(theta_x_range, wavelength_range);

fig = plot_2d_map(theta_x_range, wavelength_range, xw_gt, ...
    '\theta_x (deg)', 'Wavelength (nm)', '2D demo ground truth');
save_figure_png(fig, output_dir, '2D_ground_truth.png');

%% Simulate the encoded spectral measurement
% The forward model produces a single measured spectrum from the product of
% the object and the experimentally calibrated response.
measurement = Forward_model(xw_gt, xw_sampling_matrix);
measurement = normalize_nonnegative(measurement);

fig = figure('Color', 'w');
plot(wavelength_range, measurement, 'LineWidth', 1.4);
xlabel('Wavelength (nm)');
ylabel('Normalized intensity (a.u.)');
title('2D demo simulated measurement');
grid on; box on;
save_figure_png(fig, output_dir, '2D_measurement.png');

%% Reconstruct with FISTA and show the built-in live progress window
% FISTA_2D displays the reconstruction map and the residual curve every
% opts.display_every iterations. This animation is intentionally not saved;
% it is meant for interactive inspection while the demo runs.
opts.denoise_method = 'mixed_l1l2';
opts.iter_num = 2000;
opts.record_index = [];
opts.display_every = 100;
opts.lambda_l1 = 1e-3;
opts.lambda_l2 = 1e1;
opts.lambda_tv_x = opts.lambda_l1;
opts.lambda_tv_w = opts.lambda_l2;

if can_use_gpu()
    fprintf('[2D demo] Using GPU.\n');
    xw_sampling_matrix_solver = gpuArray(xw_sampling_matrix);
    measurement_solver = gpuArray(measurement);
else
    fprintf('[2D demo] Using CPU.\n');
    xw_sampling_matrix_solver = xw_sampling_matrix;
    measurement_solver = measurement;
end

A = @(x) Forward_model(x, xw_sampling_matrix_solver);
AT = @(y) Forward_model_transpose(y, xw_sampling_matrix_solver);

fprintf('[2D demo] Starting FISTA reconstruction.\n');
tic;
[~, res_record, rec_solver] = FISTA_2D(measurement_solver, A, AT, opts);
runtime_seconds = toc;
fprintf('[2D demo] Finished in %.1f s.\n', runtime_seconds);

rec = normalize_nonnegative(gather(rec_solver));
res_record = gather(res_record);

%% Save final data and visualization
% The final PNG summarizes the GT, reconstruction, measurement consistency,
% and mean spectral profile.
save(fullfile(output_dir, '2D_demo_result.mat'), ...
    'xw_gt', 'rec', 'measurement', 'xw_sampling_matrix', ...
    'theta_x_range', 'wavelength_range', 'opts', 'res_record', ...
    'runtime_seconds', '-v7.3');

fig = plot_2d_reconstruction_summary(theta_x_range, wavelength_range, ...
    xw_gt, rec, measurement, Forward_model(rec, xw_sampling_matrix), res_record);
save_figure_png(fig, output_dir, '2D_reconstruction_summary.png');

fprintf('[2D demo] Outputs saved to:\n  %s\n', output_dir);

function xw_gt = create_demo_2d_object(theta_x_range, wavelength_range)
    center_theta_deg = -10;
    support_width_deg = 15;
    half_width_deg = support_width_deg / 2;

    relative_theta = abs(theta_x_range(:) - center_theta_deg);
    angular_profile = zeros(numel(theta_x_range), 1);
    support_mask = relative_theta <= half_width_deg;
    angular_profile(support_mask) = 0.5 * ...
        (1 + cos(pi * relative_theta(support_mask) / half_width_deg));

    smooth_spectrum = create_random_smooth_spectrum(wavelength_range);
    xw_gt = angular_profile * smooth_spectrum(:).';
    xw_gt = normalize_nonnegative(xw_gt);
end

function spectrum = create_random_smooth_spectrum(wavelength_range)
    control_wavelengths = linspace(min(wavelength_range), max(wavelength_range), 3);
    control_values = random_high_contrast_control_values();
    spectrum = interp1(control_wavelengths, control_values, wavelength_range, 'pchip');
    spectrum = min(max(spectrum, 0), 1);
    spectrum = spectrum - min(spectrum);
    if max(spectrum) > 0
        spectrum = spectrum ./ max(spectrum);
    end
    spectrum = 0.35 + 0.65 * spectrum;
end

function control_values = random_high_contrast_control_values()
    min_contrast = 0.35;
    control_values = 0.5 + 0.5 * rand(1, 3);
    retry_count = 0;
    while range(control_values) < min_contrast && retry_count < 100
        control_values = 0.5 + 0.5 * rand(1, 3);
        retry_count = retry_count + 1;
    end
end

function fig = plot_2d_map(x_axis, y_axis, image_data, x_label, y_label, title_text)
    fig = figure('Color', 'w');
    imagesc(x_axis, y_axis, image_data');
    axis xy;
    xlabel(x_label);
    ylabel(y_label);
    title(title_text);
    colormap(parula);
    colorbar;
    box on;
end

function fig = plot_2d_reconstruction_summary(theta_x_range, wavelength_range, ...
    xw_gt, rec, measurement, measurement_rec, res_record)
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 22 14]);

    subplot(2, 3, 1);
    imagesc(theta_x_range, wavelength_range, xw_gt');
    axis xy; title('Ground truth'); xlabel('\theta_x (deg)'); ylabel('Wavelength (nm)');
    colormap(parula); colorbar; box on;

    subplot(2, 3, 2);
    imagesc(theta_x_range, wavelength_range, rec');
    axis xy; title('Reconstruction'); xlabel('\theta_x (deg)'); ylabel('Wavelength (nm)');
    colormap(parula); colorbar; box on;

    subplot(2, 3, 3);
    imagesc(theta_x_range, wavelength_range, abs(rec - xw_gt)');
    axis xy; title('Absolute error'); xlabel('\theta_x (deg)'); ylabel('Wavelength (nm)');
    colormap(parula); colorbar; box on;

    subplot(2, 3, 4);
    plot(wavelength_range, measurement, 'LineWidth', 1.2); hold on;
    plot(wavelength_range, normalize_nonnegative(measurement_rec), '--', 'LineWidth', 1.2);
    xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
    title('Measurement fit'); legend('Simulated', 'From rec', 'Box', 'off'); grid on; box on;

    subplot(2, 3, 5);
    foreground_mask = mean(xw_gt, 2) > 0.05 * max(mean(xw_gt, 2));
    gt_spectrum = normalize_curve(mean(xw_gt(foreground_mask, :), 1));
    rec_spectrum = normalize_curve(mean(rec(foreground_mask, :), 1));
    plot(wavelength_range, gt_spectrum, 'LineWidth', 1.2); hold on;
    plot(wavelength_range, rec_spectrum, '--', 'LineWidth', 1.2);
    xlabel('Wavelength (nm)'); ylabel('Normalized intensity');
    title('Foreground mean spectrum'); legend('GT', 'Rec', 'Box', 'off'); grid on; box on;

    subplot(2, 3, 6);
    semilogy(res_record, 'LineWidth', 1.2);
    xlabel('Iteration'); ylabel('Residual');
    title('FISTA loss'); grid on; box on;
end

function y = normalize_curve(y)
    y = y(:).';
    y(~isfinite(y)) = 0;
    y(y < 0) = 0;
    max_value = max(y);
    if max_value > 0
        y = y ./ max_value;
    end
end

function save_figure_png(fig, output_dir, file_name)
    ensure_dir(output_dir);
    exportgraphics(fig, fullfile(output_dir, file_name), 'Resolution', 300);
end

function ensure_dir(path_name)
    if exist(path_name, 'dir') ~= 7
        mkdir(path_name);
    end
end
