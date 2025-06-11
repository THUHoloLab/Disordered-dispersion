function Visulize_HS_cube(wavelength_range, hs_cube, fig_number)


%% RGB

% 读取 CIE 1931 数据
cie_data = readmatrix('CIE_xyz_1931_2deg.csv');
cie_wavelength = cie_data(:, 1);
cie_x = cie_data(:, 2);
cie_y = cie_data(:, 3);
cie_z = cie_data(:, 4);

% 插值到高光谱波长范围
x_interp = interp1(cie_wavelength, cie_x, wavelength_range, 'linear', 0);
y_interp = interp1(cie_wavelength, cie_y, wavelength_range, 'linear', 0);
z_interp = interp1(cie_wavelength, cie_z, wavelength_range, 'linear', 0);

% 获取高光谱图像尺寸
[rows, cols, bands] = size(hs_cube);

% 初始化 XYZ 通道
X = zeros(rows, cols);
Y = zeros(rows, cols);
Z = zeros(rows, cols);

% 按波长计算 XYZ
for i = 1:bands
    X = X + hs_cube(:, :, i) * x_interp(i);
    Y = Y + hs_cube(:, :, i) * y_interp(i);
    Z = Z + hs_cube(:, :, i) * z_interp(i);
end
xyz_image=cat(3,X,Y,Z);
rgb_image = xyz2rgb(xyz_image,'WhitePoint','c');
rgb_image=rgb_image./ max(rgb_image, [], "all");

% 显示 RGB 图像
figure;
set(gcf, 'Units', 'centimeter', 'Position', [1 1 3 3] * 2);
subplot('Position', [0.2, 0.2, 0.6, 0.6])
imshow(rgb_image);


%%
indices = round(linspace(1, length(wavelength_range), fig_number+2));
indices=indices(2:end-1);
fig_wavelength = wavelength_range(indices);

for i=1:fig_number
    figure();
    hs_color=wavelength2color(fig_wavelength(i));
    hs_colormap=linspace(0,1,256)'.*hs_color;
    image=squeeze(hs_cube(:,:,indices(i)));
    imagesc(image);
    
    colorbar
    
    set(gcf,'Units','centimeter','Position',[1 1 3 3]*2);
    subplot('Position',[0.2, 0.2, 0.6, 0.6])
    imagesc(image)
    set(gca,'xtick',[],'ytick',[])
    set(gca,'FontName','Arial','FontSize',8,'LineWidth',1);
    
  
    colormap(gca,hs_colormap)
    c = colorbar;
    c.Location="northoutside";
    clim([0, 1]);
    c.Ticks = [0,0.5,1];
    c.FontSize = 8;
    c.LineWidth=1;
    c.Position = [0.4, 0.84, 0.4, 0.02];
    c.YAxisLocation='top';
    c.Label.String = 'Intensity'; % 设置颜色条标签
    c.Label.FontSize = 8; % 设置颜色条标签字体大小
    c.Label.Position = [-0.28, -0.8]; % 设置颜色条标签位置相对于颜色条位置% 在左上角添加indices值的文本标签
    text(0.05, 0.95, sprintf('%d nm', fig_wavelength(i)), 'Units', 'normalized', 'Color', 'w', 'FontSize', 8);
end
end


