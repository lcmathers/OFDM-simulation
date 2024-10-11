clear all;
%% 参数设定

T_sim=0.1;  % 信道响应仿真时间
N=55;  % 针对起飞降落场景和巡航场景的谐波数
N_i=20; % 针对滑行和停止场景的最小谐波数
fd=3000; % 最大多普勒频率 
fmin=1000; % 巡航场景的fmin
fmax=2000; % 巡航场景的fmax

K_takeoff=15;   %起飞，降落场景的莱斯因子
K_takeoff_1= 10^(K_takeoff/10);
K_enroute=15;   % 巡航场景的莱斯因子
K_enroute_1= 10^(K_enroute/10);
K_taxi=6.9;   % 滑行场景的莱斯因子
K_taxi_1= 10^(K_taxi/10);
K_park_1=0;  % 停止场景没有LOS

f_los_takeoff=fd;  % 起飞降落场景LOS径的频偏
f_los_enroute=-fd; % 巡航场景LOS径的频偏
f_los_taxi=0.7*fd; % 滑行场景LOS径的频偏
f_los_flytower=0.5*fd;  % 飞过塔台LOS径的频偏

fs=20*1e6;    % 信号带宽 20MHz 采样率
%fs=100*max(max(f_result));  % 最大频率的100倍的采样率
Ts=1/fs; %采样周期



%% 指数衰减型功率谱的延时抽头生成
max_delay=10e-6; % 最大时延
edge_time=1e-6;  % 边沿时间
N_delay=20;  % 时延抽头数目 
[delay_actual,delay_actual_dB]=delay_pdp1(N_delay,Ts,max_delay,edge_time);  % 注意：与采样周期匹配后，由于抽头的重复，抽头的数目可能变小


%% 巡航场景是两径模型
delay_route=33*1e-6; %经典延时
delay_tap_1= [0,round(delay_route/Ts)];
delay_dB_1=[0 -10];

K=2; % 两径模型生成两组参数


%% 起飞降落场景

% 返回h(tau,t),对应的时间区间和对应的采样率（100倍的最大频率）
[mu_k,tdomain,fs1]=take_off(T_sim,fd,K_takeoff_1,f_los_takeoff,fs,delay_actual,delay_actual_dB,N);


%% 巡航场景
% [mu_k,tdomain,fs1]=en_route(T_sim,fd,fmin,fmax,K_enroute_1,f_los_enroute,fs,delay_tap_1,delay_dB_1,N);


%% 滑行场景
%[mu_k,tdomain,fs1,delay_max]=taxi(T_sim,fd,K_taxi_1,f_los_taxi,fs,delay_actual,delay_actual_dB,N_i);

%% 停泊场景
%[mu_k,tdomain,fs1,delay_max]=park(T_sim,fd,K_park_1,0,fs,delay_actual,delay_actual_dB,N_i);

%% 飞过塔台场景
% [mu_k,tdomain,fs1,delay_max]=flytower(T_sim,fd,f_los_flytower,fs);
