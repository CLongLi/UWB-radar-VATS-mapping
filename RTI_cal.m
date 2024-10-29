function [x_pos,y_pos,RTI_func] = RTI_cal(VAR,ToF_TRx,AnchorPos,re_SampTime)
AnchorPos = AnchorPos(:,1:2);
% Radio tomographic imaging
% x->(-5,5)  y->(-5,1)
%% system setting
light_speed = 299792458;
numSize = 100;
RTI_func = zeros(numSize,numSize);
VAR = mat2gray(smooth(VAR,10)');
offset_cali = 1;
for ii = 1:numSize
    for jj = 1:numSize
        x_pos(ii,jj) = -5+9.5/numSize*ii;
        y_pos(ii,jj) = -5.5+7.5/numSize*jj;
        xy_pos = [x_pos(ii,jj) y_pos(ii,jj)];
        delay_tmp = (norm(AnchorPos(1,:)-xy_pos)+norm(AnchorPos(2,:)-xy_pos))/light_speed*1e9-ToF_TRx+offset_cali;
        [~,idx_tmp] = min(abs(delay_tmp-re_SampTime));
        RTI_func(ii,jj) = VAR(idx_tmp);

    end
end
% figure;hold on;xlabel('$x$ (m)');ylabel('$y$ (m)')
% pcolor(x_pos,y_pos,RTI_func);
% plot(AnchorPos(1,1),AnchorPos(1,2),'r^');
% plot(AnchorPos(2,1),AnchorPos(2,2),'r^')

end
