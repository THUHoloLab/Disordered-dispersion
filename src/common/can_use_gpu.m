function tf = can_use_gpu()
%CAN_USE_GPU True when MATLAB can create gpuArray data in this session.

tf = false;
if exist('gpuDeviceCount', 'file') ~= 2
    return;
end

try
    tf = gpuDeviceCount > 0;
catch
    tf = false;
end
end
