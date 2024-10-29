function RTI_func = RTI_cal_fast(var_tmp01,var_tmp02,var_tmp04,...
        ToF_TRx01,ToF_TRx02,ToF_TRx04,AnchorPos,re_SampTime,x_grid,y_grid)
% function RTI_func = RTI_cal_fast(var_tmp01,var_tmp02,var_tmp04,var_tmp12,var_tmp14,var_tmp24,...
%         ToF_TRx01,ToF_TRx02,ToF_TRx04,ToF_TRx12,ToF_TRx14,ToF_TRx24,AnchorPos,re_SampTime,x_grid,y_grid)
% Radio tomographic imaging
% x->(-5,5)  y->(-5,1)
%% system setting
light_speed = 299792458;

var_tmp01 = mat2gray(smooth(var_tmp01,10)');
var_tmp02 = mat2gray(smooth(var_tmp02,10)');
var_tmp04 = mat2gray(smooth(var_tmp04,10)');
% var_tmp12 = mat2gray(smooth(var_tmp12,10)');
% var_tmp14 = mat2gray(smooth(var_tmp14,10)');
% var_tmp24 = mat2gray(smooth(var_tmp24,10)');

dis01 = sqrt((AnchorPos(1,1)-x_grid).^2+(AnchorPos(1,2)-y_grid).^2)+sqrt((AnchorPos(2,1)-x_grid).^2+(AnchorPos(2,2)-y_grid).^2);
dis02 = sqrt((AnchorPos(1,1)-x_grid).^2+(AnchorPos(1,2)-y_grid).^2)+sqrt((AnchorPos(3,1)-x_grid).^2+(AnchorPos(3,2)-y_grid).^2);
dis04 = sqrt((AnchorPos(1,1)-x_grid).^2+(AnchorPos(1,2)-y_grid).^2)+sqrt((AnchorPos(4,1)-x_grid).^2+(AnchorPos(4,2)-y_grid).^2);
% dis12 = sqrt((AnchorPos(2,1)-x_grid).^2+(AnchorPos(2,2)-y_grid).^2)+sqrt((AnchorPos(3,1)-x_grid).^2+(AnchorPos(3,2)-y_grid).^2);
% dis14 = sqrt((AnchorPos(2,1)-x_grid).^2+(AnchorPos(2,2)-y_grid).^2)+sqrt((AnchorPos(4,1)-x_grid).^2+(AnchorPos(4,2)-y_grid).^2);
% dis24 = sqrt((AnchorPos(3,1)-x_grid).^2+(AnchorPos(3,2)-y_grid).^2)+sqrt((AnchorPos(4,1)-x_grid).^2+(AnchorPos(4,2)-y_grid).^2);

delay_tmp01 = dis01(:)/light_speed*1e9-ToF_TRx01+1;
delay_tmp02 = dis02(:)/light_speed*1e9-ToF_TRx02+1;
delay_tmp04 = dis04(:)/light_speed*1e9-ToF_TRx04+1;
% delay_tmp12 = dis12(:)/light_speed*1e9-ToF_TRx12+1;
% delay_tmp14 = dis14(:)/light_speed*1e9-ToF_TRx14+1;
% delay_tmp24 = dis24(:)/light_speed*1e9-ToF_TRx24+1;

[~,idx_tmp01] = min(abs(repmat(delay_tmp01,[1 numel(re_SampTime)])-repmat(re_SampTime,[numel(delay_tmp01) 1])),[],2);
[~,idx_tmp02] = min(abs(repmat(delay_tmp02,[1 numel(re_SampTime)])-repmat(re_SampTime,[numel(delay_tmp02) 1])),[],2);
[~,idx_tmp04] = min(abs(repmat(delay_tmp04,[1 numel(re_SampTime)])-repmat(re_SampTime,[numel(delay_tmp04) 1])),[],2);
% [~,idx_tmp12] = min(abs(repmat(delay_tmp12,[1 numel(re_SampTime)])-repmat(re_SampTime,[numel(delay_tmp12) 1])),[],2);
% [~,idx_tmp14] = min(abs(repmat(delay_tmp14,[1 numel(re_SampTime)])-repmat(re_SampTime,[numel(delay_tmp14) 1])),[],2);
% [~,idx_tmp24] = min(abs(repmat(delay_tmp24,[1 numel(re_SampTime)])-repmat(re_SampTime,[numel(delay_tmp24) 1])),[],2);

RTI_func01 = var_tmp01(idx_tmp01);
RTI_func02 = var_tmp02(idx_tmp02);
RTI_func04 = var_tmp04(idx_tmp04);
% RTI_func12 = var_tmp12(idx_tmp12);
% RTI_func14 = var_tmp14(idx_tmp14);
% RTI_func24 = var_tmp24(idx_tmp24);

RTI_func = RTI_func01+RTI_func02+RTI_func04;
% RTI_func = RTI_func01+RTI_func02+RTI_func04+RTI_func12+RTI_func14+RTI_func24;
RTI_func = mat2gray(RTI_func);

end
