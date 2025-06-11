function measurement_seq = Forward_model(xyw_cube,xywm_weight)

[x,y,w,m]=size(xywm_weight);

xyw_cube=reshape(xyw_cube,[x,y,w,1]);
xywm_cube=xyw_cube.*ones(1,1,1,m);
temp=sum(xywm_cube.*xywm_weight,[1,2]);
measurement_seq=squeeze(temp);

measurement_seq=measurement_seq/m;

end




% function measurement_seq = Forward_model(xyw_cube,xw_sampling_matrix,y_response,measurement_num)
% [x,y,w]=size(xyw_cube);
% measurement_seq=zeros(w,measurement_num);
% rot_angle_seq=linspace(0,360,measurement_num+1);
% rot_angle_seq=rot_angle_seq(1:end-1);
% 
% xw_sampling_matrix=reshape(xw_sampling_matrix,[x,1,w]);
% y_response=reshape(y_response,[1,y,1]);
% weight=xw_sampling_matrix.*y_response;
% 
% for i=1:measurement_num
%     rot_angle=rot_angle_seq(i);
%     image_cube_rot=imrotate(xyw_cube,rot_angle,'crop');
%     measurement_seq(:,i)=sum(image_cube_rot.*weight,[1,2]);
% end
% 
% 
% end