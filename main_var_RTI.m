function main_var_RTI
clear all;clc;close all
%% --------------------Data loading----------------------------------------
% CIR
load('CIRReal.mat')
load('CIRImag.mat')
CIR_meas = CIRReal+1i*CIRImag;
% anchor locations
load('AnchorPos.mat')
% MU ground truth
load('posDynObjInterp.mat')
% Tx/Rx ID pairs
load('rxId.mat') % #0-1-2-4
load('txId.mat') % #0-1-2-4
% UWB recorded time (minus the first recorded time stamp)
load('tUWB.mat')
% (ns) FP offset with respect to FP sample 
load('tOffset.mat') 
%% system setting
light_speed = 299792458;
Ts = 1/(2*(4*124.8e6)); % (s) sampling period for CIR measurement
FPOffset = -3; % offset of logged samples from first path sample as determined by Decawave's leading edge detection algorithm.
numCIRSamples = 31; % number of logged CIR samples
sampT_Len = FPOffset+(numCIRSamples-0.5)*Ts*1e9;
SampleTimes = [FPOffset:Ts*1e9:sampT_Len];
%% --------------------Data processing-------------------------------------
%% background: tUWB<17.5s OR tUWB>103.5s
%  Moving Person: 21.4<tUWB & tUWB<100
% tempIDX = find(1<tUWB & tUWB<15); % background
% tempIDX = find(tUWB<21.4); % background
% tempIDX = find(21.4<tUWB & tUWB<100); % dynamic
tempIDX = find(tUWB<100);
CIR_meas = CIR_meas(tempIDX,:);
tUWB = tUWB(tempIDX);
tOffset = tOffset(tempIDX);
MU = posDynObjInterp(tempIDX,:);
txId = txId(tempIDX);
rxId = rxId(tempIDX);
%% select Tx/RX pair 0-1-2-4
TxID = [0 1 2 4];
RxID = [0 1 2 4];
tic
%% Tx/Rx-pair-0/1
ToF_TRx01 = norm(AnchorPos(1,:)-AnchorPos(2,:))/light_speed*1e9;
[CIR_meas01,MU01,tUWB01,tOffset01] = TRxPairSelection(TxID(2),RxID(1),txId,rxId,CIR_meas,MU,tUWB,tOffset);
[tau_idx,var_Bg01,var_Dyn01] = VarianceFilter(abs(CIR_meas01),tUWB01,tOffset01,SampleTimes);
%% Tx/Rx-pair-0/2
ToF_TRx02 = norm(AnchorPos(1,:)-AnchorPos(3,:))/light_speed*1e9;
[CIR_meas02,MU02,tUWB02,tOffset02] = TRxPairSelection(TxID(3),RxID(1),txId,rxId,CIR_meas,MU,tUWB,tOffset);
[~,var_Bg02,var_Dyn02] = VarianceFilter(abs(CIR_meas02),tUWB02,tOffset02,SampleTimes);
% Tx/Rx-pair-0/4
ToF_TRx04 = norm(AnchorPos(1,:)-AnchorPos(4,:))/light_speed*1e9;
[CIR_meas04,MU04,tUWB04,tOffset04] = TRxPairSelection(TxID(4),RxID(1),txId,rxId,CIR_meas,MU,tUWB,tOffset);
[~,var_Bg04,var_Dyn04] = VarianceFilter(abs(CIR_meas04),tUWB04,tOffset04,SampleTimes);
%% Tx/Rx-pair-1/2
ToF_TRx12 = norm(AnchorPos(2,:)-AnchorPos(3,:))/light_speed*1e9;
[CIR_meas12,MU12,tUWB12,tOffset12] = TRxPairSelection(TxID(3),RxID(2),txId,rxId,CIR_meas,MU,tUWB,tOffset);
[~,var_Bg12,var_Dyn12] = VarianceFilter(abs(CIR_meas12),tUWB12,tOffset12,SampleTimes);
%% Tx/Rx-pair-1/4
ToF_TRx14 = norm(AnchorPos(2,:)-AnchorPos(4,:))/light_speed*1e9;
[CIR_meas14,MU14,tUWB14,tOffset14] = TRxPairSelection(TxID(4),RxID(2),txId,rxId,CIR_meas,MU,tUWB,tOffset);
[~,var_Bg14,var_Dyn14] = VarianceFilter(abs(CIR_meas14),tUWB14,tOffset14,SampleTimes);
%% Tx/Rx-pair-2/4
ToF_TRx24 = norm(AnchorPos(3,:)-AnchorPos(4,:))/light_speed*1e9;
[CIR_meas24,MU24,tUWB24,tOffset24] = TRxPairSelection(TxID(4),RxID(3),txId,rxId,CIR_meas,MU,tUWB,tOffset);
[~,var_Bg24,var_Dyn24] = VarianceFilter(abs(CIR_meas24),tUWB24,tOffset24,SampleTimes);
%% input all CIR series, align the recording time

RTI_main(tau_idx,AnchorPos,tUWB,MU,...
                    tUWB01,tUWB02,tUWB04,tUWB12,tUWB14,tUWB24,...
                    var_Dyn01,var_Dyn02,var_Dyn04,var_Dyn12,var_Dyn14,var_Dyn24,...
                    var_Bg01,var_Bg02,var_Bg04,var_Bg12,var_Bg14,var_Bg24,...
                    ToF_TRx01,ToF_TRx02,ToF_TRx04,ToF_TRx12,ToF_TRx14,ToF_TRx24);
% RTI_main_fast(tau_idx,AnchorPos,tUWB,MU,...
%                     tUWB01,tUWB02,tUWB12,...
%                     var_Dyn01,var_Dyn02,var_Dyn12,...
%                     var_Bg01,var_Bg02,var_Bg12,...
%                     ToF_TRx01,ToF_TRx02,ToF_TRx12);

avgTime = toc/785;


end

function [Sel_CIR_meas,Sel_MU,Sel_tUWB,Sel_tOffset] = TRxPairSelection(TxID,RxID,txId,rxId,CIR_meas,MU,tUWB,tOffset)
% pinpoint the specific Tx/Rx pair
Sel_ID = find((txId==TxID & rxId==RxID) | (txId==RxID & rxId==TxID));
% Sel_ID = find((txId==TxID & rxId==RxID));
Sel_CIR_meas = CIR_meas(Sel_ID,:);
Sel_CIR_magn = sqrt(real(Sel_CIR_meas).^2+imag(Sel_CIR_meas).^2);
Sel_MU = MU(Sel_ID,:);
Sel_tUWB = tUWB(Sel_ID);
Sel_tOffset = tOffset(Sel_ID);
% outlier rejection
r_idx_tmp = find(sum(Sel_CIR_magn(:,1:3),2)<3 & sum(Sel_CIR_magn(:,1:end),2)>260 & max(Sel_CIR_magn,[],2)<70);
% reduce sampling rate
cut_ratio = 1;
r_idx = r_idx_tmp(1:cut_ratio:end,1);
Sel_CIR_meas = Sel_CIR_meas(r_idx,:);
Sel_MU = Sel_MU(r_idx,:);
Sel_tUWB = Sel_tUWB(r_idx);
Sel_tOffset = Sel_tOffset(r_idx);

end
