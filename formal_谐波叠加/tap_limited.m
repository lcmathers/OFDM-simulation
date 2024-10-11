clear all;

% 
% delay_tap=[0 0.75];
% 
% band_matrix=zeros(10,length(delay_tap));
% 
% for n=1:1:11
%     for i=1:length(delay_tap)
%         n1=n-6;
%         band_matrix(n,i)=sinc(delay_tap(i)-n1);
%     end
% end
% 
% at1=[1 1 1 1 1 1 1];
% at2=[0.1 0.1 0.1 0.1 0.1 0.1 0.1];
% 
% at=[at1;
%     at2];
% 
% 
% tap=band_matrix*at;



fs=20e6; % 采样频率
Ts=1/fs;
fd=1e3;

% fmin=250; % 巡航场景的fmin
% fmax=750; % 巡航场景的fmax

max_delay=10e-6; % 最大时延
edge_time=1e-6;  % 边沿时间
N_delay=20;  % 时延抽头数目 

[delay_actual,delay_actual_dB]=delay_pdp1(N_delay,Ts,max_delay,edge_time);

K_takeoff=15;
K_takeoff_linear=10^(K_takeoff/10);

N=55;
t_initial=0;
num_points=2024;


[mu_k_1,t_domain1,fs1,delay_tap1]= take_off(t_initial,num_points,fd,K_takeoff_linear,fd,fs,delay_actual,delay_actual_dB,N)

t_iter=t_initial+num_points*Ts;

[mu_k_2,t_domain2,fs2,delay_tap2]= take_off(t_iter,num_points,fd,K_takeoff_linear,fd,fs,delay_actual,delay_actual_dB,N)

[mu_limited,n1]=bandlimited_channel(mu_k_1,delay_tap1);


% % 
% for i=1:length(n1)
%     if n1(i)<0
%         x_data=[x(-n1(i)+1:end) zeros(1,-n1(i))]
%     else
%         x_data=[zeros(1,n1(i)) x(1:end-n1(i))];
%     end
% 
%     x_con=x_data.*mu_limited(n1(i)+101,:)
% 
% end



x=[1,2,3,4,5,6];
n1i=1;
%x_data=[x(-n1i+1:end) zeros(1,-n1i)]
 x_data=[zeros(1,n1i) x(1:end-n1i)];
x_data1=filter([0 1],1,x);


% for i=1:size(mu_limited,1)
%     if sum(abs(mu_limited(i,:)))<1e-10
%         n1(i)=0
%     end
% end

% s_takeoff=doppler('Asymmetric Jakes',[0,1]); % 半边Jakes谱
% 
% channel_takeoff=comm.RicianChannel('SampleRate',fs,'PathDelays',delay_tap*Ts,'AveragePathGains',delay_actual_dB, ...
%     'DopplerSpectrum',s_takeoff,'MaximumDopplerShift',fd,'DirectPathDopplerShift',fd,'KFactor',K_takeoff_linear,'PathGainsOutputPort',true);
% 
% x=randi([0 3],1,2048*5);
% x_tr=qammod(x,4,"gray","UnitAveragePower",true);
% 
% [y pathgains]=channel_takeoff(x_tr')

subplot(2,1,1)
plot(abs(mu_k_1(1,:)));
subplot(2,1,2)
plot(abs(mu_k_2(1,:)));
