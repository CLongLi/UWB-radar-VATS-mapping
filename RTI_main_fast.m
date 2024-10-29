function RTI_main_fast(tau_idx,AnchorPos,tUWB,MU,...
                    tUWB01,tUWB02,tUWB12,...
                    var_Dyn01,var_Dyn02,var_Dyn12,...
                    var_Bg01,var_Bg02,var_Bg12,...
                    ToF_TRx01,ToF_TRx02,ToF_TRx12)
% function RTI_main_fast(tau_idx,AnchorPos,tUWB,MU,...
%                     tUWB01,tUWB02,tUWB04,tUWB12,tUWB14,tUWB24,...
%                     var_Dyn01,var_Dyn02,var_Dyn04,var_Dyn12,var_Dyn14,var_Dyn24,...
%                     var_Bg01,var_Bg02,var_Bg04,var_Bg12,var_Bg14,var_Bg24,...
%                     ToF_TRx01,ToF_TRx02,ToF_TRx04,ToF_TRx12,ToF_TRx14,ToF_TRx24)
% (21.4<tUWB & tUWB<100); % dynamic
tempIDX01 = find(21.4<tUWB01 & tUWB01<100);
tUWB01 = tUWB01(tempIDX01);
var_Dyn01 = var_Dyn01(tempIDX01,:);

tempIDX02 = find(21.4<tUWB02 & tUWB02<100);
tUWB02 = tUWB02(tempIDX02);
var_Dyn02 = var_Dyn02(tempIDX02,:);

% tempIDX04 = find(21.4<tUWB04 & tUWB04<100);
% tUWB04 = tUWB12(tempIDX04);
% var_Dyn04 = var_Dyn04(tempIDX04,:);

tempIDX12 = find(21.4<tUWB12 & tUWB12<100);
tUWB12 = tUWB12(tempIDX12);
var_Dyn12 = var_Dyn12(tempIDX12,:);
% 
% tempIDX14 = find(21.4<tUWB14 & tUWB14<100);
% tUWB14 = tUWB14(tempIDX14);
% var_Dyn14 = var_Dyn14(tempIDX14,:);
% 
% tempIDX24 = find(21.4<tUWB24 & tUWB24<100);
% tUWB24 = tUWB24(tempIDX24);
% var_Dyn24 = var_Dyn24(tempIDX24,:);


%
tUWB_start = max([tUWB01(1,1) tUWB02(1,1) tUWB12(1,1)]);% tUWB12(1,1) tUWB14(1,1) tUWB24(1,1)]);
tUWB_end = min([tUWB01(end) tUWB02(end) tUWB12(end)]);% tUWB12(end) tUWB14(end) tUWB24(end)]);
duration = tUWB_end-tUWB_start;
time_step = 0.1;
numWin = floor(duration/time_step);
t_unify = [tUWB_start:time_step:tUWB_start+numWin*time_step];
GT_MU = zeros(numWin,3);
PF_est = zeros(numWin,2);
factor = 1;
factor2 = 1;

%% Parameters
std_RTI = 0.1;
Np = 2000;
% x_std_noise = 0.4;
v_std_noise = 0.3;

x_begin = -5; x_end = 9.5;
y_begin = -5.5; y_end = 7.5;
resol = 0.05;
[x_grid,y_grid] = ndgrid(x_begin:resol:x_end,y_begin:resol:y_end);

TranMatrix = [1 0 time_step 0;0 1 0 time_step;0 0 1 0;0 0 0 1];

for kk = 1:numWin
    [~,idx01_temp] = min(abs(tUWB01-t_unify(kk)));
    [~,idx02_temp] = min(abs(tUWB02-t_unify(kk)));
%     [~,idx04_temp] = min(abs(tUWB04-t_unify(kk)));
    [~,idx12_temp] = min(abs(tUWB12-t_unify(kk)));
%     [~,idx14_temp] = min(abs(tUWB14-t_unify(kk)));
%     [~,idx24_temp] = min(abs(tUWB24-t_unify(kk)));
    % ground truth
    [~,idx_tmp] = min(abs(tUWB-t_unify(kk)));
    GT_MU(kk,:) = MU(idx_tmp,:);
    %%
    var_tmp01 = factor2*var_Dyn01(idx01_temp,:)-factor*var_Bg01;
    var_tmp02 = factor2*var_Dyn02(idx02_temp,:)-factor*var_Bg02;
%     var_tmp04 = factor2*var_Dyn04(idx04_temp,:)-factor*var_Bg04;
    var_tmp12 = factor2*var_Dyn12(idx12_temp,:)-factor*var_Bg12;
%     var_tmp14 = factor2*var_Dyn14(idx14_temp,:)-factor*var_Bg14;
%     var_tmp24 = factor2*var_Dyn24(idx24_temp,:)-factor*var_Bg24;
    
    if kk == 1
        RTI_func = RTI_cal_fast(var_tmp01,var_tmp02,var_tmp12,...
            ToF_TRx01,ToF_TRx02,ToF_TRx12,AnchorPos,tau_idx,x_grid,y_grid);
%         RTI_func = RTI_cal_fast(var_tmp01,var_tmp02,var_tmp04,var_tmp12,var_tmp14,var_tmp24,...
%             ToF_TRx01,ToF_TRx02,ToF_TRx04,ToF_TRx12,ToF_TRx14,ToF_TRx24,AnchorPos,tau_idx,x_grid,y_grid);
        [~,idx_temp] = max(RTI_func,[],'all');
        [idx,idy] = ind2sub(size(x_grid),idx_temp);
        x_pre = x_grid(idx,idy);
        y_pre = y_grid(idx,idy);

        x_P(1,:) = x_pre+0.1*rand(1,Np);
        x_P(2,:) = y_pre+0.1*rand(1,Np);
        x_P(3,:) = zeros(1,Np);
        x_P(4,:) = zeros(1,Np);
    end
    
    x_P_update = TranMatrix*x_P+[time_step*v_std_noise*randn(1,Np);time_step*v_std_noise*randn(1,Np);v_std_noise*randn(1,Np);v_std_noise*randn(1,Np)];
    x_P_update(1,x_P_update(1,:)<-4.7) = -4.7+abs(randn);
    x_P_update(1,x_P_update(1,:)>4.2) = 4.2-abs(randn);
    x_P_update(2,x_P_update(2,:)<-5.25) = -5.25+abs(randn);
    x_P_update(2,x_P_update(2,:)>1.9) = 1.9-abs(randn);
    
    RTI_func = RTI_cal_fast(var_tmp01,var_tmp02,var_tmp12,...
                ToF_TRx01,ToF_TRx02,ToF_TRx12,AnchorPos,tau_idx,x_P_update(1,:),x_P_update(2,:));
%     RTI_func = RTI_cal_fast(var_tmp01,var_tmp02,var_tmp04,var_tmp12,var_tmp14,var_tmp24,...
%                 ToF_TRx01,ToF_TRx02,ToF_TRx04,ToF_TRx12,ToF_TRx14,ToF_TRx24,AnchorPos,tau_idx,x_P_update(1,:),x_P_update(2,:));
    % ---PF weight
    P_w = log(1/(sqrt(2*pi)*std_RTI))-0.5*(1-RTI_func).^2/std_RTI^2;
    P_w = exp(P_w-max(P_w));
    % ---normalization 
    P_w = P_w./sum(P_w);  
    % ---Resampling    
    [~, idx_tmp] = histc(rand(1,Np), cumsum(P_w,2)); 
    x_P = x_P_update(:,idx_tmp+1);
    %% ---output
    PF_est(kk,:) = mean(x_P(1:2,:)',1); 
    
end
% timeCost = toc;
% avgTimeCost = toc/numWin;
smooth_idx = 0;
if smooth_idx == 1
    PF_est(:,1) = smooth(PF_est(:,1)',floor(0.5/time_step));
    PF_est(:,2) = smooth(PF_est(:,2)',floor(0.5/time_step));
end

disERR_PF = sqrt((GT_MU(:,1)-PF_est(:,1)).^2+(GT_MU(:,2)-PF_est(:,2)).^2);
figure;hold on;
cdfdraw(disERR_PF,'color','red','LineStyle','-.','Marker','none')
xlabel('Tracking error (m)')

figure;hold on;
plot(GT_MU(:,1),GT_MU(:,2),'k.-');
plot(PF_est(:,1) ,PF_est(:,2),'rx-');
xlabel('$x$ (m)');
ylabel('$y$ (m)');
legend('Ground truth','Tracking results')
L = legend;L.ItemTokenSize(1) = 15;

figure;scatter(GT_MU(:,1),GT_MU(:,2),7,disERR_PF);
xlabel('$x$ (m)');
ylabel('$y$ (m)');
end

