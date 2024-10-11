

function [mu_k,t_domain,fs1]=en_route(T_sim,fd,fmin,fmax,K_route_1,f_los_route,fs,delay_tap,delay_dB_1,N)
%% 参数设定

% T_sim=0.5;  % 信道响应仿真时间
% N=55;  % 谐波数
% fd=3000; % 最大多普勒频率 
% fmin=1000;
% fmax=2000;
% K_route=15;     %巡航场景的莱斯因子
% K_route_1= 10^(K_route/10);
% f_los_route=-fd;


% fs=20*1e6;    % 信号带宽 20MHz 采样率
%fs=100*max(max(f_result));  % 最大频率的100倍的采样率
Ts=1/fs; %采样周期


%% 巡航场景是两径模型
% delay_route=33*1e-6; %经典延时
% delay_tap= [0,round(delay_route/Ts)];
% delay_dB_1=[0 -10];
% delay_actual=10.^(delay_dB_1/10);
K=2; % 两径模型生成两组参数
[F,C,Th]=function_MMEA_SP_enroute(N,K,fd,fmin,fmax,delay_dB_1);

%% SOC

%fs1=fs; % 信道响应的抽样率
fs1=100*max(max(F)); % 100倍的最大频率
%num_ponits=1e3;
%T_sim=num_ponits*Ts;

[mu_k,t_domain]=SOC_MMEA_SP(F,C,Th,fs1,T_sim,K_route_1,f_los_route);


%% 画出h(tau,t)
[X,Y]=meshgrid(t_domain,delay_tap);
figure;
plot3(X(1,:),Y(1,:)*Ts,abs(mu_k(1,:)));
hold on;
plot3(X(2,:),Y(2,:)*Ts,abs(mu_k(2,:)));
xlabel('time');
ylabel('delay');
zlabel('amplitude');

%% 画出多普勒功率谱 和第一径相比较

S_theory=Doppler_spectrum(fd,1024);
fdomain=linspace(0,fd,512);

% 只画出关心的部分Jakes谱
match1=abs(round(fmin/fd*512));
match2=abs(round(fmax/fd*512));
match1_1=max(match1,match2);
match2_2=min(match1,match2);
figure;
S_theory([1:match2_2-1])=0;
S_theory([match1_1+1:512])=0;
plot(fdomain,S_theory(1:512)*1e2);  %理论

hold on;
stem(F(1,:),C(1,:).^2);  % 实际的仿真
xlabel('frequency (Hz)');
ylabel('Amplitude');



%% 自相关函数的比较

delay_max=N/(2*fd); % 关心的最大时延
%delay_max=0.1;
Ns1=1e3;  % 生成理论自相关函数的点数
T1=delay_max/Ns1;
delay_tau=[0:T1:(Ns1-1)*T1];
t1=0;
h_theory=zeros(1,Ns1);  % 理论自相关函数
h_model=zeros(K,Ns1);   % 时间遍历自相关函数
alpha1=acos(fmin/fd);
alpha2=acos(fmax/fd);
alpha=abs(alpha2-alpha1);



% 理论自相关函数

for i=1:Ns1
    AngleFunction=@(angle) exp(1j*(2*pi*fd*t1*cos(angle)));
    h_theory(i) = 1/alpha*integral(AngleFunction,min(alpha1,alpha2),max(alpha1,alpha2));
    t1=t1+T1;
end


% 时间遍历的自相关函数

t1=0;
for i=1:Ns1
    for k=1:K
        h_model(k,i)=sum(exp(1j*2*pi*F(k,:)*t1)./N);
    end
    t1=t1+T1;
end



% 时限（样本函数）的自相关函数 以第二径为列


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



% 画出对比图
figure;
plot(delay_tau,abs(h_theory),'b-','LineWidth',1);
hold on;
plot(delay_tau,abs(h_model(1,:)),'k--','LineWidth',1);
plot(lags_selected,abs(acf_selected),'m--','LineWidth',1);
plot(lags_selected,abs(acf_selected1),'r--','LineWidth',1);
ylim([0,1]);
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

end
