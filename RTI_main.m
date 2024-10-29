function RTI_main(tau_idx,AnchorPos,tUWB,MU,...
                    tUWB01,tUWB02,tUWB04,tUWB12,tUWB14,tUWB24,...
                    var_Dyn01,var_Dyn02,var_Dyn04,var_Dyn12,var_Dyn14,var_Dyn24,...
                    var_Bg01,var_Bg02,var_Bg04,var_Bg12,var_Bg14,var_Bg24,...
                    ToF_TRx01,ToF_TRx02,ToF_TRx04,ToF_TRx12,ToF_TRx14,ToF_TRx24)
% (21.4<tUWB & tUWB<100); % dynamic
tempIDX01 = find(21.4<tUWB01 & tUWB01<100);
tUWB01 = tUWB01(tempIDX01);
var_Dyn01 = var_Dyn01(tempIDX01,:);

tempIDX02 = find(21.4<tUWB02 & tUWB02<100);
tUWB02 = tUWB02(tempIDX02);
var_Dyn02 = var_Dyn02(tempIDX02,:);

tempIDX04 = find(21.4<tUWB04 & tUWB04<100);
tUWB04 = tUWB04(tempIDX04);
var_Dyn04 = var_Dyn04(tempIDX04,:);

tempIDX12 = find(21.4<tUWB12 & tUWB12<100);
tUWB12 = tUWB12(tempIDX12);
var_Dyn12 = var_Dyn12(tempIDX12,:);

tempIDX14 = find(21.4<tUWB14 & tUWB14<100);
tUWB14 = tUWB14(tempIDX14);
var_Dyn14 = var_Dyn14(tempIDX14,:);

tempIDX24 = find(21.4<tUWB24 & tUWB24<100);
tUWB24 = tUWB24(tempIDX24);
var_Dyn24 = var_Dyn24(tempIDX24,:);


%
tUWB_start = max([tUWB01(1,1) tUWB02(1,1) tUWB04(1,1) tUWB12(1,1) tUWB14(1,1) tUWB24(1,1)]);
tUWB_end = min([tUWB01(end) tUWB02(end) tUWB04(end) tUWB12(end) tUWB14(end) tUWB24(end)]);
duration = tUWB_end-tUWB_start;
time_step = 0.1;
numWin = floor(duration/time_step);
t_unify = [tUWB_start:time_step:tUWB_start+numWin*time_step];
x_est = zeros(1,numWin);
y_est = zeros(1,numWin);
GT_MU = zeros(numWin,3);
PF_est = zeros(numWin,2);
factor = 1;
factor2 = 1;

%% tracking-fig
figure;hold on;
filename = 'Tracking RTI_4anchors.gif';
[ttt,sss] = title('Ground Truth (black) vs VATS Mapping (white)','4 anchors');
ttt.FontSize = 12;
sss.FontSize = 12;
set(gcf,'color','w');
xlabel('$x$ (m)','interpreter','latex');ylabel('$y$ (m)','interpreter','latex')
plot(AnchorPos(1,1),AnchorPos(1,2),'w^','LineWidth',2,'MarkerSize',10);
plot(AnchorPos(2,1),AnchorPos(2,2),'w^','LineWidth',2,'MarkerSize',10);
plot(AnchorPos(3,1),AnchorPos(3,2),'w^','LineWidth',2,'MarkerSize',10);
plot(AnchorPos(4,1),AnchorPos(4,2),'w^','LineWidth',2,'MarkerSize',10);
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
xlim([-5, 4.5])
ylim([-5.5, 2])
%% Parameters
std_RTI = 0.1;
Np = 2000;
% x_std_noise = 0.4;
v_std_noise = 0.3;

x_P_update = zeros(2,Np);
x_P(1,:) = -4.5+8.5*rand(1,Np);
x_P(2,:) = -4.5+6*rand(1,Np);
x_P(3,:) = zeros(1,Np);
x_P(4,:) = zeros(1,Np);

P_w = zeros(1,Np);
tic 
for kk = 1:numWin
    [~,idx01_temp] = min(abs(tUWB01-t_unify(kk)));
    [~,idx02_temp] = min(abs(tUWB02-t_unify(kk)));
    [~,idx04_temp] = min(abs(tUWB04-t_unify(kk)));
    [~,idx12_temp] = min(abs(tUWB12-t_unify(kk)));
    [~,idx14_temp] = min(abs(tUWB14-t_unify(kk)));
    [~,idx24_temp] = min(abs(tUWB24-t_unify(kk)));
    % ground truth
    [~,idx_tmp] = min(abs(tUWB-t_unify(kk)));
    GT_MU(kk,:) = MU(idx_tmp,:);
    %%
    [x_pos,y_pos,RTI_func01_Dyn] = RTI_cal(factor2*var_Dyn01(idx01_temp,:)-factor*var_Bg01,ToF_TRx01,AnchorPos([1 2],:),tau_idx);
    [~,~,RTI_func02_Dyn] = RTI_cal(factor2*var_Dyn02(idx02_temp,:)-factor*var_Bg02,ToF_TRx02,AnchorPos([1 3],:),tau_idx);
    [~,~,RTI_func04_Dyn] = RTI_cal(factor2*var_Dyn04(idx04_temp,:)-factor*var_Bg04,ToF_TRx04,AnchorPos([1 4],:),tau_idx);
    [~,~,RTI_func12_Dyn] = RTI_cal(factor2*var_Dyn12(idx12_temp,:)-factor*var_Bg12,ToF_TRx12,AnchorPos([2 3],:),tau_idx);
    [~,~,RTI_func14_Dyn] = RTI_cal(factor2*var_Dyn14(idx14_temp,:)-factor*var_Bg14,ToF_TRx14,AnchorPos([2 4],:),tau_idx);
    [~,~,RTI_func24_Dyn] = RTI_cal(factor2*var_Dyn24(idx24_temp,:)-factor*var_Bg24,ToF_TRx24,AnchorPos([3 4],:),tau_idx);
    % 3 anchors
%     RTI_func = RTI_func01_Dyn+RTI_func02_Dyn+RTI_func12_Dyn;
%     RTI_func = RTI_func01_Dyn+RTI_func04_Dyn+RTI_func14_Dyn;
%     RTI_func = RTI_func02_Dyn+RTI_func04_Dyn+RTI_func24_Dyn;
%     RTI_func = RTI_func12_Dyn+RTI_func14_Dyn+RTI_func24_Dyn;
    % 4 anchors
    RTI_func = RTI_func01_Dyn+RTI_func02_Dyn+RTI_func04_Dyn+RTI_func12_Dyn+RTI_func14_Dyn+RTI_func24_Dyn;
    RTI_func = mat2gray(RTI_func);
    
    for iter = 1:Np
        x_P_update(1,iter) = x_P(1,iter) + time_step*x_P(3,iter) + time_step*v_std_noise*randn;  
        x_P_update(2,iter) = x_P(2,iter) + time_step*x_P(4,iter) + time_step*v_std_noise*randn;
        % particles limits
        if x_P_update(1,iter)>4.2
            x_P_update(1,iter) = 4.2;
        elseif x_P_update(1,iter)<-4.7
            x_P_update(1,iter) = -4.7;
        end
        if x_P_update(2,iter)>1.9
            x_P_update(2,iter) = 1.9;
        elseif x_P_update(2,iter)<-5.25
            x_P_update(2,iter) = -5.25;
        end

        x_P_update(3,iter) = x_P(3,iter) + v_std_noise*randn;  
        x_P_update(4,iter) = x_P(4,iter) + v_std_noise*randn;
        
        [~,idx_x] = min(abs(x_P_update(1,iter)-x_pos(:,1)));
        [~,idx_y] = min(abs(x_P_update(2,iter)-y_pos(1,:)));
        P_w(iter) = log(1/(sqrt(2*pi)*std_RTI))-0.5*(1-RTI_func(idx_x,idx_y))^2/std_RTI^2;
    end
    P_w = exp(P_w - max(P_w));
    % normalization 
    P_w = P_w./sum(P_w);  
    %% Resampling    
    [~, idx_tmp] = histc(rand(1,Np), cumsum(P_w,2)); 
    x_P = x_P_update(:,idx_tmp+1);
    PF_est(kk,:) = mean(x_P(1:2,:)',1); 
    
    hold on;
    xlim([-5, 4.5])
    ylim([-5.5, 2])
    xlabel('$x$ (m)','interpreter','latex');ylabel('$y$ (m)','interpreter','latex')
    h = pcolor(x_pos,y_pos,abs(RTI_func));
    set(h, 'EdgeColor', 'none');
    plot(x_P(1,:),x_P(2,:),'rx','MarkerSize',8);  
    plot(AnchorPos(1,1),AnchorPos(1,2),'w^','LineWidth',2,'MarkerSize',10);
    plot(AnchorPos(2,1),AnchorPos(2,2),'w^','LineWidth',2,'MarkerSize',10);
    plot(AnchorPos(3,1),AnchorPos(3,2),'w^','LineWidth',2,'MarkerSize',10);
    plot(AnchorPos(4,1),AnchorPos(4,2),'w^','LineWidth',2,'MarkerSize',10);
    plot(GT_MU(1:kk,1),GT_MU(1:kk,2),'k.-','LineWidth',2);
    plot(PF_est(1:kk,1),PF_est(1:kk,2),'w.-','LineWidth',2);

    set(gca,'FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if kk == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.01,'WriteMode','append');
    end
    
    
    [idx_x,idx_y] = find(RTI_func == max(RTI_func,[],'all'));
    x_est(kk) = x_pos(idx_x(1,1),idx_y(1,1));
    y_est(kk) = y_pos(idx_x(1,1),idx_y(1,1));
    
    kk/numWin
    kk
    disERR = sqrt((GT_MU(kk,1)'-x_est(kk)).^2+(GT_MU(kk,2)'-y_est(kk)).^2)
    
end
timeCost = toc;
avgTimeCost = toc/numWin;
smooth_idx = 0;
if smooth_idx == 1
    x_est = smooth(x_est,floor(0.5/time_step))';
    y_est = smooth(y_est,floor(0.5/time_step))';
    PF_est(:,1) = smooth(PF_est(:,1)',floor(0.5/time_step));
    PF_est(:,2) = smooth(PF_est(:,2)',floor(0.5/time_step));
end

disERR = sqrt((GT_MU(:,1)'-x_est).^2+(GT_MU(:,2)'-y_est).^2);
disERR_PF = sqrt((GT_MU(:,1)-PF_est(:,1)).^2+(GT_MU(:,2)-PF_est(:,2)).^2);
figure;hold on;
cdfdraw(disERR,'color','black','LineStyle','-.','Marker','none')
cdfdraw(disERR_PF,'color','red','LineStyle','-.','Marker','none')
xlabel('Tracking error (m)')

figure;hold on
plot(GT_MU(:,1),GT_MU(:,2),'k.-');
plot(x_est,y_est),'rx-';
xlabel('$x$ (m)');
ylabel('$y$ (m)');
legend('Ground truth','Tracking results')
L = legend;L.ItemTokenSize(1) = 15;

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

