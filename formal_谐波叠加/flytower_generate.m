clear all;

%% 飞过塔台场景只有一个 有频偏的直射径
%% 参数设定

T_sim=0.01;  % 信道响应仿真时间
fd=3000; % 最大多普勒频率 
f_los_flytower=0.5*fd;  % LOS径的频偏 0.5倍的最大频偏

%fs=20*1e6;    % 信号带宽 20MHz 采样率
fs=100*fd;  % 最大频率的100倍的采样率
Ts=1/fs; %采样周期
Ns=round(T_sim/Ts);
mu_k=zeros(1,Ns);
tdomain=[0:Ts:(Ns-1)*Ts];

for l=0:Ns-1
mu_k(l+1)=exp(1j*2*pi*f_los_flytower*l*Ts);
end

figure
plot(tdomain,abs(mu_k));
xlabel('time');
ylabel('Amplitude');

% 时限自相关函数和多普勒功率谱

delay_max=1/(2*f_los_flytower);

% 计算数据的均值和标准差
mu = mean(mu_k);
sigma = std(mu_k);

% 对数据进行归一化处理
normalized_data = (mu_k - mu) / sigma;

[h_simulation,lags]=xcorr(normalized_data,'biased')  % biased estimator
[h_simulation1,lags1]=xcorr(normalized_data,'unbiased')  % biased estimator
lags=lags/fs;
lag_range = lags >= 0 & lags <= delay_max;
acf_selected = h_simulation(lag_range);
acf_selected1=h_simulation1(lag_range);
lags_selected = lags(lag_range);
h_theory=ones(1,length(lags_selected));

figure;
plot(lags_selected,abs(h_theory),'r-','LineWidth',1);
hold on;
plot(lags_selected,abs(acf_selected),'m--','LineWidth',1);
ylim([0,1])
xlabel('delay');
ylabel('Amplitude');
legend('theory','time sampled BE')


figure;
Doppler=fftshift(fft(h_simulation));
fdomain=linspace(0,fd,length(Doppler));
plot(fdomain,abs(Doppler));
xlabel('frequency');
ylabel('Amplitude');

