function xyw_cube = Forward_model_transpose(measurement_seq,xywm_weight_transpose)

[x,y,w,m]=size(xywm_weight_transpose);

temp=reshape(measurement_seq,[1,1,w,m]);
temp=temp.*ones(x,y,1,1);
temp=sum(temp.*xywm_weight_transpose,4);
xyw_cube=squeeze(temp);

xyw_cube=xyw_cube/x/y;

end





% function xyw_cube = Forward_model_transpose(measurement_seq,xw_sampling_matrix,y_response)
% 
% [x,w]=size(xw_sampling_matrix);
% [w,measurement_num]=size(measurement_seq);
% y=length(y_response);
% xyw_cube=zeros(x,y,w);
% 
% rot_angle_seq=linspace(0,360,measurement_num+1);
% rot_angle_seq=rot_angle_seq(1:end-1);
% 
% xw_sampling_matrix=reshape(xw_sampling_matrix,[x,1,w]);
% y_response=reshape(y_response,[1,y,1]);
% weight=xw_sampling_matrix.*y_response;
% 
% 
% for i=1:measurement_num
%     rot_angle=rot_angle_seq(i);
%     measurement=reshape(measurement_seq(:,i),[1,1,w]);
%     temp=ones(x,y,1).*measurement/x/y;
%     temp=temp.*weight;
%     temp=imrotate(temp,rot_angle,'crop');
%     xyw_cube=xyw_cube+temp;
% end
% xyw_cube=xyw_cube./measurement_num;
