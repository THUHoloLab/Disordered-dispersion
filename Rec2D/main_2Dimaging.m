clear;clc
addpath('src_2D\')

addpath('data_2D\')

% Define wavelength and spatial angle ranges
wavelength_range = 400:0.5:700; % Wavelength range in nm
theta_x_range = linspace(-25, 25, 50); % X-angle range in degrees, 50 points
theta_y_range = linspace(-25, 25, 50); % Y-angle range in degrees, 50 points

% Get lengths of the ranges
x_len = length(theta_x_range); % Number of points in x-dimension
y_len = length(theta_y_range); % Number of points in y-dimension
w_len = length(wavelength_range); % Number of wavelength points
measurement_num = 36; % Number of measurements (e.g., different rotation angles)


%% LOAD CALIBRATION DATA (System Response)
% Load pre-calibrated system response data (Transmission matrix T)

calibration_data_path = "..\measured_response_data.mat";
load(calibration_data_path, ...
    'response_xyw','theta_x_data','theta_y_data','wavelength_data');

% Data manipulation: Augment theta_y_data and response_xyw
original_theta_y_data = theta_y_data;
theta_y_data_augmented = cat(2, -50, -40, -30, -20, -18, -14, -flip(original_theta_y_data,2), original_theta_y_data(:,2:end), 14, 18, 20, 30, 40, 50);
xyw_response_data_augmented = cat(2, zeros(length(theta_x_data),6,length(wavelength_data)), flip(response_xyw,2), response_xyw(:,2:end,:), zeros(length(theta_x_data),6,length(wavelength_data)));

% Using augmented data for interpolation
theta_y_data_for_interp = theta_y_data_augmented;
xyw_response_data_for_interp = xyw_response_data_augmented;


%% Load Ground Truth (GT) Source Spectrum and Dark Spectrum Data

% Source spectrum
source_spectrum_path = ".\data_2D\Source.csv";
data = readmatrix(source_spectrum_path);
source_spectrum = mean(data(:, 2:end),2); % Average if multiple measurements
source_spectrum = source_spectrum - min(source_spectrum); % Baseline correction
source_spectrum = source_spectrum./(max(source_spectrum)); % Normalize to [0,1]

% Dark spectrum for measurements
dark_spectrum_path = ".\data_2D\300ms_dark.csv";
data_dark = readmatrix(dark_spectrum_path);
spectrum_dark = data_dark(:, 2:end);
spectrum_dark = mean(spectrum_dark,2); % Average if multiple measurements


%% Interpolation of System Response
% Interpolate the loaded system response to the desired sampling grid


[XX_data,YY_data,WW_data] = ndgrid(theta_x_data, theta_y_data_for_interp, wavelength_data);
[XX_samples,YY_samples,WW_samples] = ndgrid(theta_x_range, theta_y_range, wavelength_range);
xyw_sampling_matrix = interpn(XX_data, YY_data, WW_data, xyw_response_data_for_interp, XX_samples, YY_samples, WW_samples, 'spline');

 
% % Dense Interpolation (alternative, currently commented out)
% % This approach first interpolates to a denser grid, then averages.
% dense_sample_factor = 2; % Factor for denser grid
% theta_x_range_dense = linspace(min(theta_x_range), max(theta_x_range), dense_sample_factor*length(theta_x_range)); 
% theta_y_range_dense = linspace(min(theta_y_range), max(theta_y_range), dense_sample_factor*length(theta_y_range)); 
% 
% [XX_dense_samples,YY_dense_samples,WW_dense_samples_constW] = ndgrid(theta_x_range_dense,theta_y_range_dense,wavelength_range);
% 
% xyw_dense_sampling_matrix_interp = interpn(XX_data,YY_data,WW_data,xyw_response_data_for_interp,XX_dense_samples,YY_dense_samples,WW_dense_samples_constW,'spline');
% 
% xyw_sampling_matrix_avg = zeros(x_len,y_len,w_len);
% for i = 1:x_len
%     for j = 1:y_len
%         % Define block for averaging
%         block = xyw_dense_sampling_matrix_interp((i-1)*dense_sample_factor+1 : i*dense_sample_factor, ...
%                                            (j-1)*dense_sample_factor+1 : j*dense_sample_factor, :);
%         xyw_sampling_matrix_avg(i,j,:) = mean(block,[1,2]); % Average over spatial dimensions of the block
%     end
% end
% xyw_sampling_matrix = xyw_sampling_matrix_avg; % Use the averaged matrix


%% Post-process the interpolated sampling matrix
xyw_sampling_matrix(xyw_sampling_matrix < 0) = 0;       % Ensure non-negativity
xyw_sampling_matrix(isnan(xyw_sampling_matrix)) = 0;   % Replace NaNs with 0
xyw_sampling_matrix = xyw_sampling_matrix ./ max(xyw_sampling_matrix(:)); % Normalize to [0,1]

%% Calculate Weight Matrix (System Matrix for different rotations)
% This creates a set of system matrices, each corresponding to a rotation of the base sampling matrix.
xywm_weight = zeros(x_len, y_len, w_len, measurement_num); % Preallocate
rot_angle_seq = linspace(0, 360, measurement_num + 1); % Rotation angles
rot_angle_seq = rot_angle_seq(1:end-1); % Exclude 360 (same as 0)

disp('Generating rotated system matrices...');
for i = 1:measurement_num
    rot_angle = rot_angle_seq(i);
    % Rotate each wavelength slice of the sampling matrix
    current_rotated_matrix = zeros(x_len,y_len,w_len);
    for wl_idx = 1:w_len
         current_rotated_matrix(:,:,wl_idx) = imrotate(xyw_sampling_matrix(:,:,wl_idx), rot_angle, 'bilinear', 'crop');
    end
    xywm_weight(:,:,:,i) = current_rotated_matrix;
end
disp('Rotated system matrices generated.');

% Optional: Visualize a slice of the weight matrix
figure()
imagesc(wavelength_range,theta_x_range,squeeze(xywm_weight(:,round(y_len/2),:,1))); colorbar;
xlabel('Wavelength'); ylabel('X-angle');
title('Slice of Rotated System Matrix (Angle 1, Center Y)');
colormap(parula); % Use a default MATLAB colormap

%% Load Simulation or Experiment data

FLAG='Experiment';
switch FLAG
    case 'Simulation'

        % Create a simple example hyperspectral cube
        spectrum_cube_gt = Create_simple_spectrum_cube(x_len, y_len, w_len); % Assumes this function is defined
        % Interpolate/resize to match dimensions if necessary (already done by Create_simple_spectrum_cube if sizes match)
        obj_sim = spectrum_cube_gt;
        obj_sim(obj_sim < 0) = 0; % Ensure non-negativity
        obj_sim = obj_sim ./ max(obj_sim(:)); % Normalize
        
        measurements = Forward_model(obj_sim, xywm_weight); % Simulate measurements
        Visulize_HS_cube(wavelength_range, obj_sim, 4); % Assumes this function is defined
        title('Simulated Ground Truth Object');
    case 'Experiment'
        % Load Experimental Measurements
        experiment_data_path = ".\windmill_300ms_mea1.csv";
        data = readmatrix(experiment_data_path);
        w_data_experiment = data(:,1); % Wavelengths from spectrometer
        wm_data_raw = data(:, 2:end);  % Raw measurement data (spectra x measurement_angle)
        
        % Process measurements
        wm_data_processed = wm_data_raw - spectrum_dark; % Subtract dark spectrum
        % Resize/interpolate if number of measurements doesn't match measurement_num
        if size(wm_data_processed, 2) ~= measurement_num
            disp(['Resizing measurement data from ', num2str(size(wm_data_processed, 2)), ' to ', num2str(measurement_num), ' columns.']);
            wm_data_processed = imresize(wm_data_processed, [length(w_data_experiment), measurement_num], 'bilinear');
        end
        
        measurements = zeros(w_len, measurement_num); % Preallocate
        for i=1:measurement_num
            measurements(:,i) = interp1(w_data_experiment, wm_data_processed(:,i), wavelength_range, 'linear', 'extrap');
        end
        measurements(measurements < 0) = 0; % Ensure non-negativity
        measurements = measurements ./ max(measurements(:)); % Normalize
end
%% Reconstruction Parameters and Algorithm
opts.denoise_method = 'mixed_l1l2';
opts.iter_num = 4000; % Number of iterations (reduced for faster domo)
opts.record_index = []; % Iteration indices to record (if any)
opts.display_every = 100; % Display progress every N iterations
opts.lambda_l1 = 1e-4; % Example L1 regularization parameter (adjust as needed)
opts.lambda_l2 = 1e-5; % Example L2 regularization parameter (for mixed_l1l2)


% Move data to GPU if available and desired
if gpuDeviceCount > 0
    xywm_weight_gpu = gpuArray(xywm_weight);
    measurements_gpu = gpuArray(measurements);
    disp('Using GPU for reconstruction.');
else
    xywm_weight_gpu = xywm_weight;
    measurements_gpu = measurements;
    disp('No GPU detected, using CPU for reconstruction.');
end

% Define forward and transpose operators for the reconstruction
A = @(x_op) (Forward_model(x_op, xywm_weight_gpu)); % Forward model
AT = @(y_op) (Forward_model_transpose(y_op, xywm_weight_gpu)); % Transpose of forward model (adjoint)

% Run FISTA algorithm
disp('Starting FISTA reconstruction...');
tic; % Start timer
[x_record, res_record, v_reconstructed_gpu] = FISTA_3D(measurements_gpu, A, AT, opts); % Call FISTA
toc % Stop timer

% Retrieve reconstructed data from GPU
rec = gather(v_reconstructed_gpu);
rec = rec ./ max(rec(:)); % Normalize the final reconstruction

%% Save Reconstruction and Measurement Data
% save('Sunflower_rec.mat',"rec","wavelength_range","theta_x_range","theta_y_range");
% save('Sunflower_mea.mat',"measurements","wavelength_range","rot_angle_seq")
disp('Reconstruction and measurement data can be saved (currently commented out).');

%% Visualize Reconstructed Hyperspectral Cube

if exist('Visulize_HS_cube', 'file')
    Visulize_HS_cube(wavelength_range, rec, 4); % Example: display 4 slices
    sgtitle('Reconstructed Hyperspectral Cube'); % Add super title
else
    disp('Visulize_HS_cube function not found. Skipping visualization.');
    % Basic visualization if Visulize_HS_cube is not available:
    figure;
    subplot(1,3,1); imagesc(theta_x_range, theta_y_range, rec(:,:,round(w_len/4))); axis image; title(['Slice at ~', num2str(wavelength_range(round(w_len/4))), 'nm']); colormap(parula); colorbar;
    subplot(1,3,2); imagesc(theta_x_range, theta_y_range, rec(:,:,round(w_len/2))); axis image; title(['Slice at ~', num2str(wavelength_range(round(w_len/2))), 'nm']); colormap(parula); colorbar;
    subplot(1,3,3); imagesc(theta_x_range, theta_y_range, rec(:,:,round(3*w_len/4))); axis image; title(['Slice at ~', num2str(wavelength_range(round(3*w_len/4))), 'nm']); colormap(parula); colorbar;
    sgtitle('Reconstructed Hyperspectral Cube (Basic Slices)');
end


%% Draw PSNR/SSIM Curves (If Ground Truth 'obj_sim' is available from simulation)
% This section is relevant if you run the simulation part to get 'obj_sim'
% if exist('obj_sim', 'var')
%     PSNR_val = zeros(1, w_len);
%     SSIM_val = zeros(1, w_len);
%     for i=1:w_len
%         rec_slice = double(rec(:,:,i));
%         obj_slice = double(obj_sim(:,:,i)); % Ensure obj_sim is the ground truth
%         PSNR_val(i) = psnr(rec_slice, obj_slice);
%         SSIM_val(i) = ssim(rec_slice, obj_slice);
%     end
% 
%     figure()
%     set(gcf, 'Units', 'centimeter', 'Position', [2 2 10 8]); % Adjusted size
%     subplot('Position', [0.15, 0.15, 0.7, 0.7]) % Standard subplot position
%     yyaxis left
%     plot(wavelength_range, PSNR_val,'-.','LineWidth', 1.5)
%     ylabel("PSNR (dB)")
%     ylim([min(PSNR_val)-1 max(PSNR_val)+1]) % Adjust y-axis limits
% 
%     yyaxis right
%     plot(wavelength_range, SSIM_val,'LineWidth', 1.5)
%     ylabel("SSIM")
%     ylim([min(SSIM_val)-0.05 max(SSIM_val)+0.05]) % Adjust y-axis limits
% 
%     legend('PSNR','SSIM','Box','off','Location','best')
%     xlabel("Wavelength (nm)")
%     title('PSNR and SSIM vs. Wavelength')
%     ax = gca;
%     ax.TickDir = 'out';
%     set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1);
%     ax.Box = 'on'; % Keep box for clarity
%     grid on;
% else
%     disp('Ground truth "obj_sim" not available, skipping PSNR/SSIM plot.');
% end


%% Draw Spectrum Slices from Reconstructed Data

figure()
set(gcf, 'Units', 'centimeter', 'Position', [2 2 12 9]); % Adjusted size
subplot('Position', [0.15, 0.15, 0.75, 0.75]) % Standard subplot

% Define points for spectral slices (indices)
pt1_x_idx = round(x_len * 0.38); % Example: 19 for x_len=50
pt1_y_idx = round(y_len * 0.64); % Example: 32 for y_len=50
pt2_x_idx = round(x_len * 0.64); 
pt2_y_idx = round(y_len * 0.64);
pt3_x_idx = round(x_len * 0.64);
pt3_y_idx = round(y_len * 0.38);
pt4_x_idx = round(x_len * 0.38);
pt4_y_idx = round(y_len * 0.38);

spectrum_slice_1_rec = squeeze(rec(pt1_x_idx, pt1_y_idx, :));
spectrum_slice_2_rec = squeeze(rec(pt2_x_idx, pt2_y_idx, :));
spectrum_slice_3_rec = squeeze(rec(pt3_x_idx, pt3_y_idx, :));
spectrum_slice_4_rec = squeeze(rec(pt4_x_idx, pt4_y_idx, :));

hold on
plot(wavelength_range, spectrum_slice_1_rec, 'linestyle','-', 'LineWidth', 1.5, 'Marker','o', 'MarkerIndices',round(linspace(1,w_len,7)), 'MarkerSize',5)
plot(wavelength_range, spectrum_slice_2_rec, 'linestyle',':', 'LineWidth', 1.5, 'Marker','s', 'MarkerIndices',round(linspace(1,w_len,7)), 'MarkerSize',5)
plot(wavelength_range, spectrum_slice_3_rec, 'linestyle','--', 'LineWidth', 1.5, 'Marker','^', 'MarkerIndices',round(linspace(1,w_len,7)), 'MarkerSize',5)
plot(wavelength_range, spectrum_slice_4_rec, 'linestyle','-.', 'LineWidth', 1.5, 'Marker','d', 'MarkerIndices',round(linspace(1,w_len,7)), 'MarkerSize',5)

% Create legend strings dynamically
legend_str = {sprintf('Rec: X=%d, Y=%d', theta_x_range(pt1_x_idx), theta_y_range(pt1_y_idx)), ...
              sprintf('Rec: X=%d, Y=%d', theta_x_range(pt2_x_idx), theta_y_range(pt2_y_idx)), ...
              sprintf('Rec: X=%d, Y=%d', theta_x_range(pt3_x_idx), theta_y_range(pt3_y_idx)), ...
              sprintf('Rec: X=%d, Y=%d', theta_x_range(pt4_x_idx), theta_y_range(pt4_y_idx))};
leg = legend(legend_str);
leg.Box = "off";
leg.Location = 'northwest';
ylabel("Normalized Intensity (a.u.)")
xlabel("Wavelength (nm)")
title('Reconstructed Spectra at Selected Spatial Points')
ax = gca;
ax.TickDir = 'out';
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1);
ax.Box = 'on'; % Keep box
axis xy;
ax.YLim = [0, 1.05*max([spectrum_slice_1_rec; spectrum_slice_2_rec; spectrum_slice_3_rec; spectrum_slice_4_rec])]; % Dynamic YLim
if all([spectrum_slice_1_rec; spectrum_slice_2_rec; spectrum_slice_3_rec; spectrum_slice_4_rec] == 0)
    ax.YLim = [0,1]; % Handle case where all spectra are zero
end
grid on;
hold off

%% Draw Measurements Curves
measurements_norm = measurements ./ max(measurements(:)); % Ensure normalization for plotting
plotted_line_indices = [1, round(measurement_num/2), measurement_num-1]; % Example lines to highlight

figure()

hold on
plotted_handles = []; % Store handles for legend
legend_labels = {};   % Store labels for legend

for i = 1:measurement_num
   h_temp = plot(wavelength_range, measurements_norm(:,i), 'LineWidth', 1);
   plotted_handles(end+1) = h_temp;
   legend_labels{end+1} = ['Meas. Angle: ', num2str(rot_angle_seq(i), '%.1f'), '^{\circ}'];
end
ylabel("Normalized Power (a.u.)")
xlabel("Wavelength (nm)")
title('Measured Spectra at Different Rotation Angles')
ax = gca;
ax.TickDir = 'out';
ax.Box = 'on';
axis xy;
grid on;
hold off

if ~isempty(plotted_handles)
    legend(plotted_handles, legend_labels, 'Box', 'off', 'Location', 'best');
end

disp('Finished processing main_2Dimaging.m');
