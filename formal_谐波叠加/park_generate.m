clear all;

%% 参数设定

T_sim=0.1;  % 信道响应仿真时间
N_i=20;  % 最小的谐波数 每一径的谐波叠加的数目不同
fd=3000; % 最大多普勒频率 
K_park=0;   %起飞，降落场景的莱斯因子 (dB)

f_los_park=0;  % LOS径的频偏 0.7倍的最大频偏

fs=20*1e6;    % 信号带宽 20MHz 采样率
%fs=100*max(max(f_result));  % 最大频率的100倍的采样率
Ts=1/fs; %采样周期


%% 指数衰减型功率谱的延时抽头生成
max_delay=10e-6; % 最大时延
edge_time=1e-6;  % 边沿时间
N_delay=20;
[delay_actual,delay_actual_dB]=delay_pdp1(N_delay,Ts,max_delay,edge_time);  % 注意：与采样周期匹配后，由于抽头的重复，抽头的数目可能变小
delay_linear=10.^(delay_actual_dB/10);
num_of_taps=length(delay_actual);


%% MEDS 参数生成
max_Ni=N_i+2*num_of_taps-1;
num_of_Ni=max_Ni-1:-2:N_i;

[F1,F2,C1,C2,Th1,Th2]=MEDS_parameter(num_of_Ni,fd,delay_actual_dB);


% 直射分量的影响

%% SOS 谐波叠加

%fs=20*1e6;
fs1=100*max(max(max(F1)),max(max(F2)));  %100倍最大频率的采样率
%fs1=100*fd;
[mu_k,t_domain]=SOS_MEDS(F1,F2,C1,C2,Th1,Th2,K_park,f_los_park,T_sim,fs1);

% 时延抽头匹配
[delay_tap,mu_k]=delay_match(delay_actual,mu_k,Ts);

%% 画出h(tau,t)
[X,Y]=meshgrid(t_domain,delay_tap);
figure;
surf(X,Y*Ts,abs(mu_k));
xlabel('time');
ylabel('delay');
zlabel('amplitude(dB)');

%% 画出功率谱 和第一径比较

S_theory=Doppler_spectrum(fd,1024);
S_theory=fftshift(S_theory);
fdomain=linspace(-fd,fd,1024);

figure;
plot(fdomain,S_theory(1:1024)*10);
hold on;
stem(F1(1,[1:end-1]),C1(1,[1:end-1]).^2/4);
stem(-F1(1,[1:end-1]),C1(1,[1:end-1]).^2/4);
stem(-F2(1,:),C2(1,:).^2/4);
stem(F2(1,:),C2(1,:).^2/4);
xlabel('frequency (Hz)');
ylabel('Amplitude');


%% 自相关函数

delay_max=max_Ni/(2*fd);
Ns=1e3;
Ts_delay=delay_max/Ns;
t=linspace(0,delay_max,Ns);


h_theory=besselj(0,2*pi*fd*t); % Jakes 谱理论自相关函数


h_simulation_1=zeros(1,Ns);  
h_simulation_2=zeros(1,Ns);


% 理论遍历时间自相关函数
for l=0:Ns-1
    for i=1:max_Ni
    h_simulation_1(l+1)=h_simulation_1(l+1)+C1(2,i)^2/2.*cos(2*pi*F1(2,i)*l*Ts_delay);
    h_simulation_2(l+1)=h_simulation_2(l+1)+C2(2,i)^2/2.*cos(2*pi*F2(2,i)*l*Ts_delay);
    end
end

% 实际样本自相关函数 以第一径为列

% 计算数据的均值和标准差
mu = mean(mu_k(2,:));
sigma = std(mu_k(2,:));

% 对数据进行归一化处理
normalized_data = (mu_k(2,:) - mu) / sigma;

[h_simulation,lags]=xcorr(normalized_data,'biased')  % biased estimator
[h_simulation1,lags1]=xcorr(normalized_data,'unbiased')  % biased estimator
lags=lags/fs1;
lag_range = lags >= 0 & lags <= delay_max;
acf_selected = h_simulation(lag_range);
acf_selected1=h_simulation1(lag_range);
lags_selected = lags(lag_range);


figure;
plot(t,abs(h_theory),'-r');
hold on;
plot(t,abs(h_simulation_1+h_simulation_2),'-g');
plot(lags_selected,abs(acf_selected),'m--','LineWidth',1);
plot(lags_selected,abs(acf_selected1),'r--','LineWidth',1);
xlabel('delay');
ylabel('Amplitude');
legend('theory','ergodic','time sampled BE','time sampled UE');


% 概率分布
figure;
h=histogram(abs(mu_k(1,:)),'Normalization','pdf');
xlabel('Amplitude');
ylabel('Probability Density');
title('PDF of the Amplitude of mu_k(1,:)');

figure;
h1=histogram(abs(mu_k(2,:)),'Normalization','pdf');
xlabel('Amplitude');
ylabel('Probability Density');
title('PDF of the Amplitude of mu_k(2,:)');