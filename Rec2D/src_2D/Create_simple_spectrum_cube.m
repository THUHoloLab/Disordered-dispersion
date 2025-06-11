function spectrum_cube=Create_simple_spectrum_cube(d_x,d_y,d_lambda)
    
    spectrum = 0.5*rand(5, 1)+0.5;
    x = 1:length(spectrum);
    xi = linspace(1, length(spectrum), d_lambda);
    spectrum_interpolated = interp1(x, spectrum, xi, 'spline');
    spectrum_interpolated=spectrum_interpolated./(max(spectrum_interpolated));

    img=imread("char_A.jpeg");
    img=double(img);
    img=img./(max(max(img)));
    img=imresize(img,[d_x,d_y]);
    
    spectrum_cube=zeros(d_x,d_y,d_lambda);
    for i=1:d_lambda
        spectrum_cube(:,:,i)=img.*spectrum_interpolated(i);
    end

    spectrum_cube(spectrum_cube<0)=0;
end

