function w_measurement = Forward_model(xw_cube,xw_sampling_matrix)
    w_measurement=mean(xw_cube.*xw_sampling_matrix,1);
end
