clear;clc
addpath('src_1D\')
addpath('data_1D\')
% Define wavelength range and spatial angle range
wavelength_range = 400:0.5:700; % Wavelength range in nm
theta_x_range = -60:1:60;  % Spatial angle range in degrees

% Get lengths of the ranges
x_len = length(theta_x_range); % Number of spatial points
w_len = length(wavelength_range); % Number of wavelength points

%% LOAD Response data
% Load pre-calibrated system response data

response_data_path = "..\measured_response_data.mat";
load(response_data_path, ...
    'response_xyw','theta_x_data','theta_y_data','wavelength_data');

% Extract the 1D response data
xw_response_data = squeeze(response_xyw(:,1,:));


%% Load Ground Truth (GT) and Dark Spectrum Data

% Source spectrum
source_spectrum_path = ".\data_1D\Source.csv";
data = readmatrix(source_spectrum_path);
source_spectrum = mean(data(:, 2:end),2); 
source_spectrum = source_spectrum - min(source_spectrum); 
source_spectrum = source_spectrum./(max(source_spectrum)); 

% Dark spectrum for 400ms exposure
dark_400ms_path = ".\data_1D\400ms_dark.csv";
data_400ms_dark = readmatrix(dark_400ms_path);
spectrum_400ms_dark = data_400ms_dark(:, 2:end);
spectrum_400ms_dark = mean(spectrum_400ms_dark,2); 

%% Interpolation of Response Data

[XX_data,WW_data] = ndgrid(theta_x_data,wavelength_data);
[XX_samples,WW_samples] = ndgrid(theta_x_range,wavelength_range);

% For more accurate averaging, first interpolate to a denser grid
% XX_dense_samples = imresize(XX_samples,[10*x_len,10*w_len],'bilinear'); 
% WW_dense_samples = imresize(WW_samples,[10*x_len,10*w_len],'bilinear');
% Using direct interpolation to target resolution for simplicity here.
% If the dense sampling and averaging block below is preferred, uncomment it and the lines above.

% Perform spline interpolation
% xw_dense_sampling_matrix = interpn(XX_data,WW_data,xw_response_data,XX_dense_samples,WW_dense_samples,'spline');
xw_sampling_matrix = interpn(XX_data,WW_data,xw_response_data,XX_samples,WW_samples,'spline');


% % Averaging from dense grid (if used)
% xw_sampling_matrix=zeros(x_len,w_len);
% for i = 1:x_len
%     for j = 1:w_len
%         block=xw_dense_sampling_matrix((i-1)*10+1:i*10, (j-1)*10+1:j*10);
%         xw_sampling_matrix(i,j)=mean(block,'all');
%     end
% end

% Post-processing for the sampling matrix
xw_sampling_matrix(xw_sampling_matrix<0) = 0; % Ensure non-negativity
xw_sampling_matrix(isnan(xw_sampling_matrix)) = 0; % Replace NaNs with 0

xw_sampling_matrix = xw_sampling_matrix./max(max(xw_sampling_matrix)); % Normalize

% Display the processed sampling matrix
DrawFig(theta_x_range, wavelength_range, xw_sampling_matrix, "x position", "Wavelength (nm)")
title('Interpolated System Response Matrix (A)') % Add title to the figure

%% Load Simulation or Experiment data

FLAG='Experiment';
switch FLAG
    case 'Simulation'
        %Simulate a measurement if a ground truth object (xw_slice) is known.
        load("stuffed_toys.mat"); % Example data
        xw_slice=double(squeeze(spectrum_cube(:,60,:)));
        
        xw_slice=imresize(xw_slice,[x_len,w_len]);
        xw_slice(xw_slice<0)=0;
        xw_slice=xw_slice./(max(max(xw_slice)));
        
        measurement= Forward_model(xw_slice,xw_sampling_matrix); % Simulate measurement
        DrawFig(theta_x_range,wavelength_range,xw_slice,"x position","Wavelength(nm)")
        title('Simulated Ground Truth Object')
    case 'Experiment'
        % Load and process actual experimental data
        experiment_data_path = ".\data_1D\1D_BlueFilter_4mmslit_z=40mm_x=-60-2-60_400ms.csv";
        data = readmatrix(experiment_data_path);
        w_data_experiment = data(:, 1); % Wavelengths from spectrometer
        measurement_seq = data(:, 2:end); % All measurements
        
        % Select a specific measurement (e.g., 35th column)
        measurement = measurement_seq(:,35); 
        measurement = measurement - spectrum_400ms_dark; % Subtract dark spectrum
        measurement = interp1(w_data_experiment, measurement, wavelength_range, 'linear', 'extrap'); % Interpolate to common wavelength range
        measurement(measurement<0) = 0; % Ensure non-negativity
        measurement = measurement./max(max(measurement)); % Normalize
end
% Plot the measurement
figure()
plot(wavelength_range, measurement)
xlabel("Wavelength (nm)")
ylabel("Intensity (a.u.)")
ax.TickDir = 'in';
title('Spectra measurement')
axis xy


%% Reconstruction Parameters and Algorithm
% Set options for the FISTA reconstruction algorithm
opts.denoise_method = 'mixed_l1l2'; 
opts.iter_num = 4000; % Number of iterations
opts.record_index = []; % Iteration indices to record (if any)
opts.display_every = 200; % Display progress every N iterations
opts.lambda_l1 = 1e-3; % Example L1 regularization parameter (adjust as needed)
opts.lambda_l2 = 1e-4; % Example L2 regularization parameter (for mixed_l1l2)
opts.lambda_tv = 1e-2; % Example TV regularization parameter (if using TV)


% Move data to GPU if available and desired
if gpuDeviceCount > 0
    xw_sampling_matrix_gpu = gpuArray(xw_sampling_matrix);
    measurement_gpu = gpuArray(measurement);
else
    xw_sampling_matrix_gpu = xw_sampling_matrix;
    measurement_gpu = measurement;
    disp('No GPU detected, using CPU.');
end


% Define forward and transpose operators for the reconstruction
A = @(x_op) (Forward_model(x_op, xw_sampling_matrix_gpu)); % Forward model operator
AT = @(y_op) (Forward_model_transpose(y_op, xw_sampling_matrix_gpu));  % Transpose of forward model (adjoint)

% Run FISTA algorithm
disp('Starting FISTA reconstruction...');
tic; % Start timer
[x_record, res_record, v_reconstructed_gpu] = FISTA_2D(measurement_gpu, A, AT, opts); % Call FISTA
toc % Stop timer

% Retrieve reconstructed data from GPU
rec = gather(v_reconstructed_gpu);
rec = rec./max(rec(:)); % Normalize reconstructed image

%% Display Reconstructed Image and Slices
DrawFig(theta_x_range, wavelength_range, rec, "x position", "Wavelength (nm)")
title('Reconstructed Image')

% Plot mean spectrum (summed/averaged over x positions)
figure()
plot(wavelength_range, mean(rec,1))
title('Mean Spectrum of Reconstructed Image')
xlabel('Wavelength (nm)')
ylabel('Mean Intensity (a.u.)')
grid on;

% Plot a spatial slice at a specific wavelength (e.g., 10th wavelength index)
if w_len >= 10
    figure()
    plot(theta_x_range, rec(:,10))
    title(['Spatial Slice at Wavelength = ', num2str(wavelength_range(10)), ' nm'])
    xlabel('x position (degrees)')
    ylabel('Intensity (a.u.)')
    grid on;
else
    disp('Not enough wavelength points to plot 10th slice.');
end


% Plot spectral slices at specific x positions 
figure()
hold on
if x_len >=20
    plot(wavelength_range, rec(20,:))
end
title('Spectral Slices at Specific x-positions')
xlabel('Wavelength (nm)')
ylabel('Intensity (a.u.)')
legend_str = {};
if x_len >= 20, legend_str{end+1} = ['x = ', num2str(theta_x_range(20))]; end
if ~isempty(legend_str), legend(legend_str); end
grid on;
hold off

%% Function to Draw Figure
% This function visualizes a 2D matrix (image)
function DrawFig(x_axis_data, y_axis_data, image_data, xaxis_label_str, yaxis_label_str)
    figure()
    imagesc(x_axis_data, y_axis_data, image_data') 
    xlabel(xaxis_label_str)
    ylabel(yaxis_label_str)
    set(gca,'FontName','Arial','FontSize',10,'LineWidth',1); % Increased font size for readability
    box on;    
    axis xy;         
    colormap(parula);
    colorbar; 
end
