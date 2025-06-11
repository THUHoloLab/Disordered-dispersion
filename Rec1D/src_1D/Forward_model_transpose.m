function xw_output = Forward_model_transpose(w_measurement,xw_sampling_matrix)
    [x,w]=size(xw_sampling_matrix);
    w_measurement=reshape(w_measurement,[1,w]);
    xw_output=ones(x,1)*w_measurement.*xw_sampling_matrix;
end