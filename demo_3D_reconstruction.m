clear; clc; close all;

% DEMO_3D_RECONSTRUCTION
% This demo reconstructs a hyperspectral cube x(theta_x, theta_y, lambda).
% The encoder response is loaded from an experimentally calibrated
% transmission matrix. The object and encoded measurements are simulated.

%% Set paths and full-size demo parameters
% The demo uses the same grid as the revised simulation script:
%   50 x 50 angular samples, 601 wavelength samples, 36 measurements.
demo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(demo_root, 'src', 'common'));
addpath(fullfile(demo_root, 'src', 'rec3d'));

data_dir = fullfile(demo_root, 'data');
output_dir = fullfile(demo_root, 'outputs', '3D_demo');
ensure_dir(output_dir);

calibration_path = ensure_demo_calibration_file(data_dir);
wavelength_range = 400:0.5:700;
theta_x_range = linspace(-25, 25, 50);
theta_y_range = linspace(-25, 25, 50);
measurement_num = 36;
rng(2);

x_len = numel(theta_x_range);
y_len = numel(theta_y_range);
w_len = numel(wavelength_range);

%% Load and resample the experimentally calibrated 3D response
% The measured response covers positive theta_y values. The original
% reconstruction code mirrors and pads this response to support the full
% angular field of view used in the 3D simulation.
load(calibration_path, 'response_xyw', 'theta_x_data', 'theta_y_data', 'wavelength_data');

theta_y_aug = cat(2, -50, -40, -30, -20, -18, -14, ...
    -flip(theta_y_data, 2), theta_y_data(:, 2:end), 14, 18, 20, 30, 40, 50);
xyw_response_aug = cat(2, ...
    zeros(numel(theta_x_data), 6, numel(wavelength_data)), ...
    flip(response_xyw, 2), response_xyw(:, 2:end, :), ...
    zeros(numel(theta_x_data), 6, numel(wavelength_data)));

[XX_data, YY_data, WW_data] = ndgrid(theta_x_data, theta_y_aug, wavelength_data);
[XX_samples, YY_samples, WW_samples] = ndgrid(theta_x_range, theta_y_range, wavelength_range);
xyw_sampling_matrix = interpn(XX_data, YY_data, WW_data, xyw_response_aug, ...
    XX_samples, YY_samples, WW_samples, 'spline');
xyw_sampling_matrix = normalize_nonnegative(xyw_sampling_matrix);

%% Build rotated measurement weights
% Each measurement corresponds to one rotation angle of the calibrated
% response. This follows the 3D experimental acquisition model.
xywm_weight = zeros(x_len, y_len, w_len, measurement_num);
rot_angle_seq = linspace(0, 360, measurement_num + 1);
rot_angle_seq = rot_angle_seq(1:end - 1);

fprintf('[3D demo] Building rotated response stack.\n');
for m = 1:measurement_num
    for wl_idx = 1:w_len
        xywm_weight(:, :, wl_idx, m) = imrotate( ...
            xyw_sampling_matrix(:, :, wl_idx), rot_angle_seq(m), 'bilinear', 'crop');
    end
end

%% Construct a simulated hyperspectral ground truth object
% The cube uses a char_A image pattern as the angular-spatial amplitude.
% Every angular position shares the same smooth spectrum generated from
% three random control points in [0.5, 1]. The random points are resampled
% until their contrast is visually clear.
pattern_path = fullfile(demo_root, 'char_A.jpeg');
obj = create_demo_3d_object(theta_x_range, theta_y_range, wavelength_range, pattern_path);

fig = plot_hsi_summary(theta_x_range, theta_y_range, wavelength_range, obj, ...
    '3D demo ground truth');
save_figure_png(fig, output_dir, '3D_ground_truth.png');

%% Simulate the encoded measurements
% The forward model maps the hyperspectral cube to a set of measured spectra,
% one spectrum per rotation angle.
measurements = Forward_model(obj, xywm_weight);
measurements = normalize_nonnegative(measurements);

fig = plot_measurement_stack(wavelength_range, rot_angle_seq, measurements);
save_figure_png(fig, output_dir, '3D_measurements.png');

%% Reconstruct with FISTA and show the built-in live progress window
% FISTA_3D displays spatial/spectral slices and the residual curve every
% opts.display_every iterations. This animation is intentionally not saved.
opts.denoise_method = 'mixed_l1l2';
opts.iter_num = 1200;
opts.record_index = [];
opts.display_every = 50;
opts.lambda_l1 = 1e-6;
opts.lambda_l2 = 1e-0;
opts.lambda_tv_xy = opts.lambda_l1;
opts.lambda_tv_w = opts.lambda_l2;

if can_use_gpu()
    fprintf('[3D demo] Using GPU.\n');
    xywm_weight_solver = gpuArray(xywm_weight);
    measurements_solver = gpuArray(measurements);
else
    fprintf('[3D demo] Using CPU.\n');
    xywm_weight_solver = xywm_weight;
    measurements_solver = measurements;
end

A = @(x) Forward_model(x, xywm_weight_solver);
AT = @(y) Forward_model_transpose(y, xywm_weight_solver);

fprintf('[3D demo] Starting FISTA reconstruction.\n');
tic;
[~, res_record, rec_solver] = FISTA_3D(measurements_solver, A, AT, opts);
runtime_seconds = toc;
fprintf('[3D demo] Finished in %.1f s.\n', runtime_seconds);

rec = normalize_nonnegative(gather(rec_solver));
res_record = gather(res_record);

%% Save final data and visualization
% The final PNG summarizes selected wavelength slices, mean spectra, and the
% residual curve. Full data are saved in the MAT file.
save(fullfile(output_dir, '3D_demo_result.mat'), ...
    'obj', 'rec', 'measurements', 'xywm_weight', ...
    'theta_x_range', 'theta_y_range', 'wavelength_range', 'rot_angle_seq', ...
    'opts', 'res_record', 'runtime_seconds', '-v7.3');

fig = plot_3d_reconstruction_summary(theta_x_range, theta_y_range, wavelength_range, ...
    obj, rec, measurements, Forward_model(rec, xywm_weight), res_record);
save_figure_png(fig, output_dir, '3D_reconstruction_summary.png');

fprintf('[3D demo] Outputs saved to:\n  %s\n', output_dir);

function obj = create_demo_3d_object(theta_x_range, theta_y_range, wavelength_range, pattern_path)
    if exist(pattern_path, 'file') ~= 2
        error('The 3D demo pattern image was not found: %s', pattern_path);
    end

    pattern_image = imread(pattern_path);
    pattern_image = double(pattern_image);
    if ndims(pattern_image) == 3
        pattern_image = mean(pattern_image, 3);
    end
    pattern_image = normalize_nonnegative(pattern_image);

    % The bundled char_A image is usually dark ink on a bright background.
    % Invert it when needed so the character itself becomes the object.
    if mean(pattern_image(:)) > 0.5
        pattern_image = 1 - pattern_image;
    end

    spatial_pattern = imresize(pattern_image, ...
        [numel(theta_x_range), numel(theta_y_range)], 'bilinear');
    spatial_pattern = normalize_nonnegative(spatial_pattern);

    smooth_spectrum = create_random_smooth_spectrum(wavelength_range);
    obj = spatial_pattern .* reshape(smooth_spectrum, 1, 1, []);
    obj = normalize_nonnegative(obj);
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

function fig = plot_hsi_summary(theta_x_range, theta_y_range, wavelength_range, cube, title_text)
    slice_wavelengths = [450, 530, 610, 680];
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 8]);
    for k = 1:numel(slice_wavelengths)
        [~, wl_idx] = min(abs(wavelength_range - slice_wavelengths(k)));
        subplot(1, numel(slice_wavelengths), k);
        imagesc(theta_x_range, theta_y_range, cube(:, :, wl_idx)');
        axis image xy;
        xlabel('\theta_x (deg)');
        if k == 1
            ylabel('\theta_y (deg)');
        end
        title(sprintf('%s nm', num2str(wavelength_range(wl_idx))));
        colormap(parula);
        colorbar;
        box on;
    end
    sgtitle(title_text);
end

function fig = plot_measurement_stack(wavelength_range, rot_angle_seq, measurements)
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 14 9]);

    subplot(1, 2, 1);
    imagesc(rot_angle_seq, wavelength_range, measurements);
    axis xy;
    xlabel('Rotation angle (deg)');
    ylabel('Wavelength (nm)');
    title('Simulated measurements');
    colormap(parula);
    colorbar;
    box on;

    subplot(1, 2, 2);
    plot(wavelength_range, measurements(:, 1:4:end), 'LineWidth', 0.9);
    xlabel('Wavelength (nm)');
    ylabel('Intensity (a.u.)');
    title('Selected measurement spectra');
    grid on;
    box on;
end

function fig = plot_3d_reconstruction_summary(theta_x_range, theta_y_range, wavelength_range, ...
    obj, rec, measurements, measurements_rec, res_record)
    slice_wavelengths = [450, 530, 610, 680];
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 16]);

    for k = 1:numel(slice_wavelengths)
        [~, wl_idx] = min(abs(wavelength_range - slice_wavelengths(k)));

        subplot(4, numel(slice_wavelengths), k);
        imagesc(theta_x_range, theta_y_range, obj(:, :, wl_idx)');
        axis image xy; title(sprintf('GT %g nm', wavelength_range(wl_idx)));
        colormap(parula); colorbar; box on;

        subplot(4, numel(slice_wavelengths), k + numel(slice_wavelengths));
        imagesc(theta_x_range, theta_y_range, rec(:, :, wl_idx)');
        axis image xy; title(sprintf('Rec %g nm', wavelength_range(wl_idx)));
        colormap(parula); colorbar; box on;

        subplot(4, numel(slice_wavelengths), k + 2 * numel(slice_wavelengths));
        imagesc(theta_x_range, theta_y_range, abs(rec(:, :, wl_idx) - obj(:, :, wl_idx))');
        axis image xy; title(sprintf('Abs. err %g nm', wavelength_range(wl_idx)));
        colormap(parula); colorbar; box on;
    end

    subplot(4, 4, 13);
    spatial_gt = mean(obj, 3);
    foreground_mask = spatial_gt(:) > 0.05 * max(spatial_gt(:));
    obj_matrix = reshape(obj, [], numel(wavelength_range));
    rec_matrix = reshape(rec, [], numel(wavelength_range));
    gt_spectrum = normalize_curve(mean(obj_matrix(foreground_mask, :), 1));
    rec_spectrum = normalize_curve(mean(rec_matrix(foreground_mask, :), 1));
    plot(wavelength_range, gt_spectrum, 'LineWidth', 1.2); hold on;
    plot(wavelength_range, rec_spectrum, '--', 'LineWidth', 1.2);
    xlabel('Wavelength (nm)'); ylabel('Normalized intensity');
    title('Foreground mean spectrum'); legend('GT', 'Rec', 'Box', 'off'); grid on; box on;

    subplot(4, 4, 14);
    plot(wavelength_range, measurements(:, 1), 'LineWidth', 1.1); hold on;
    plot(wavelength_range, normalize_nonnegative(measurements_rec(:, 1)), '--', 'LineWidth', 1.1);
    xlabel('Wavelength (nm)'); ylabel('Intensity');
    title('Measurement 1 fit'); legend('Sim', 'Rec', 'Box', 'off'); grid on; box on;

    subplot(4, 4, 15:16);
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
