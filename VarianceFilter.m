function [tau_idx,var_h_Bg,var_h_Dyn] = VarianceFilter(Z_mag,tUWB,tOffset,SampleTimes)
% beta_c -> variance scaling parameter
[numSample,numBin] = size(Z_mag);
delta_taus = 1;
m = 4; % spacing
delta_tau_knot = delta_taus/m;
%% CIR ini
meas_Time = ones(numSample,1)*SampleTimes-tOffset'*ones(1,numBin);
tau_start = -4; % min(meas_Time,[],'all');
tau_end = 27.25; % max(meas_Time,[],'all');
tau_start_end = tau_end-tau_start;
N_knot = ceil(tau_start_end/delta_tau_knot)+1;

tau_idx = zeros(1,N_knot);
h_idx = zeros(1,N_knot);
for i = 1:N_knot
    idx = i-1;
    tau_idx(i) = tau_start+delta_tau_knot*idx;
    if tau_idx(i)<meas_Time(1,1)
        h_idx(i) = Z_mag(1,1);
    elseif tau_idx(i)>meas_Time(1,end) || tau_idx(i)==meas_Time(1,end)
        h_idx(i) = Z_mag(1,end);
    else % 
        idx_tmp = discretize(tau_idx(i), meas_Time(1,:));
        h_idx(i) = Z_mag(1,idx_tmp);
    end
end
%%
alpha = 0.05;
y = zeros(numSample,numBin);
h_tau = zeros(numSample,numBin);
var_h = 4*ones(1,N_knot-1);
var_h_Dyn = zeros(numSample,N_knot-1);
for k = 1:numSample
    if  tUWB(k)<5.5 || tUWB(k)==5.5 % k~1000
        var_alpha = 0.1;
    elseif tUWB(k)>5.5 && tUWB(k)<17.5
        var_alpha = 0.1;
    else % dynamic
        var_alpha = 0.1;
    end
    for j = 1:numBin
        idx_tmp = discretize(meas_Time(k,j), tau_idx);
        h_tau(k,j) = (1-(meas_Time(k,j)-tau_idx(idx_tmp))/delta_tau_knot)*h_idx(idx_tmp)+...
            ((meas_Time(k,j)-tau_idx(idx_tmp))/delta_tau_knot)*h_idx(idx_tmp+1); 
        y(k,j) = Z_mag(k,j)-h_tau(k,j);
        h_idx(idx_tmp) = h_idx(idx_tmp)+alpha*(1-(meas_Time(k,j)-tau_idx(idx_tmp))/delta_tau_knot)*2*y(k,j);
        h_idx(idx_tmp+1) = h_idx(idx_tmp+1)+alpha*((meas_Time(k,j)-tau_idx(idx_tmp))/delta_tau_knot)*2*y(k,j);
        var_h(idx_tmp) = var_h(idx_tmp)+var_alpha*(abs(y(k,j))-var_h(idx_tmp));
    end
    % background variance
    if tUWB(k)>10 && tUWB(k)<15
        var_h_Bg = var_h;
    else
%         var_h_Bg = [];
    end
    % dynamic variance
    var_h_Dyn(k,:) = var_h;
    
    if tUWB(k)>100
        break
    end
end
tau_idx = tau_idx(1:N_knot-1);
end
