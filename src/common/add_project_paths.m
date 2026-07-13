function add_project_paths(reconstruction_mode)
%ADD_PROJECT_PATHS Add only the folders needed by one reconstruction mode.
%
% reconstruction_mode:
%   "2D" - spectral-angle reconstruction, output rec(theta_x, wavelength)
%   "3D" - hyperspectral cube reconstruction, output rec(theta_x, theta_y, wavelength)

root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root_dir, 'src', 'common'));

switch string(reconstruction_mode)
    case "2D"
        addpath(fullfile(root_dir, 'src', 'rec2d'));
    case "3D"
        addpath(fullfile(root_dir, 'src', 'rec3d'));
    otherwise
        error('Unknown reconstruction_mode: %s', reconstruction_mode);
end
end
