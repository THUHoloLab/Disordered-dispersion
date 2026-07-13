function calibration_path = ensure_demo_calibration_file(data_dir)
%ENSURE_DEMO_CALIBRATION_FILE Rebuild the calibration MAT file if needed.
% The full calibration file is larger than the GitHub web upload limit, so
% the demo package stores it as numbered binary parts. This function joins
% those parts on first use and returns the path used by LOAD.

calibration_path = fullfile(data_dir, 'measured_response_data.mat');
if isfile(calibration_path)
    return;
end

part_files = dir(fullfile(data_dir, 'measured_response_data.mat.part*'));
if isempty(part_files)
    error('Calibration file was not found: %s', calibration_path);
end

[~, order] = sort({part_files.name});
part_files = part_files(order);

tmp_path = [calibration_path, '.tmp'];
fid_out = fopen(tmp_path, 'wb');
if fid_out < 0
    error('Could not create calibration file: %s', tmp_path);
end

buffer_size = 8 * 1024 * 1024;
for idx = 1:numel(part_files)
    part_path = fullfile(part_files(idx).folder, part_files(idx).name);
    fid_in = fopen(part_path, 'rb');
    if fid_in < 0
        fclose(fid_out);
        error('Could not read calibration part: %s', part_path);
    end

    while true
        bytes = fread(fid_in, buffer_size, '*uint8');
        if isempty(bytes)
            break;
        end
        fwrite(fid_out, bytes, 'uint8');
    end
    fclose(fid_in);
end

fclose(fid_out);
movefile(tmp_path, calibration_path, 'f');
end
