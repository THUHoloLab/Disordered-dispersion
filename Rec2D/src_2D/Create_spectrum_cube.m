%% Import data cube
folder_path = 'C:\Users\Yuchen Ma\Desktop\Endoscropy\datas\stuffed_toys_ms';
% 获取文件夹中的所有图片文件
image_files = dir(fullfile(folder_path, '*.png')); % 假设你的图片格式为 jpg，如果是其他格式，请相应修改

% 初始化一个空的三维数组来存储图片数据
num_images = length(image_files);
first_image = imread(fullfile(folder_path, image_files(1).name));
[height, width, ~] = size(first_image);

% 读取每张图片并存储到三维数组中
for i = 1:num_images
    image_path = fullfile(folder_path, image_files(i).name);
    current_image = imread(image_path);
    spectrum_cube(:, :, i) = current_image;
end

save_file_name = 'stuffed_toys.mat';

% 使用 save 函数保存 all_images 到 .mat 文件中
save(save_file_name, 'spectrum_cube');
