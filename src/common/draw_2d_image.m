function draw_2d_image(x_axis, y_axis, image_data, x_label, y_label, fig_title)
%DRAW_2D_IMAGE Display a 2D reconstruction or response matrix.

figure();
imagesc(x_axis, y_axis, image_data');
xlabel(x_label);
ylabel(y_label);
title(fig_title);
set(gca, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1);
box on;
axis xy;
colormap(parula);
colorbar;
end
