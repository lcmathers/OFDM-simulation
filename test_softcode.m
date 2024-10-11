clc;clear all; close all;
rng default;

%% Parameters setting

Nfft=2048;      
Ng=Nfft/4;    % CP length
% Ng=3;
Nvc=400;        % Vitural carrier
Nframe=5;

M_mod=16;  % Modulation order

channeltype=3; % 1 滑行 2 停泊 3 起飞/降落 4 巡航
fd=0; % Doppler shift
fmin=100;  % enroute min Doppler
fmax=500; % enroute max Doppler

R_code=1/2; % channel code rate

EbN0=[0:5:30];

% SNR_dB=[0:5:30];
ExtraNoise=0; % Extra noise sample
EndNoise=0;  % End noise sample
CFO=0;

trellis = poly2trellis(7,[171 133]);

max_iter=10; % number of iter

Nps=2; % the space of pilot symbol


[STS,LTS]=generate_train(Nfft,Ng);
Train_sym=[STS LTS];
% Train_sym=[STS];
%Train_sym=[];

total_biterrors=zeros(1,length(EbN0));
total_biterrors_hard=zeros(1,length(EbN0));

BER_ofdm=zeros(1,length(EbN0));
BER_ofdm_hard=zeros(1,length(EbN0));



%% Parameters CalCulation

Nsym=Nfft+Ng;
Ndata=Nfft-Nvc;

[X_pilot,pilot_loc]=generate_pilot(Nfft,Nps,Nvc,Ndata);

M_bit=log2(M_mod);
R1=(Nfft/Nps)/(Nfft); % the ratio of pilot symbol and Nfft

sym_perframe=Ndata*Nframe*R1;
bits_perframe=M_bit*Ndata*Nframe*R1;
rawbit_perframe=bits_perframe*R_code;

SNR_dB=EbN0+10*log10(R_code*M_bit*(Ndata/Nfft));

eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2)); % 调制信号平均功率
A=1/eng_sqrt;   % QAM归一化因子

delay_result=zeros(length(SNR_dB),max_iter);
CFO_est=zeros(length(SNR_dB),max_iter);
FFO_est=zeros(length(SNR_dB),max_iter);
STO=zeros(length(SNR_dB),max_iter);


%% Main Loop

for i=1:length(SNR_dB)
    
    for iter=1:max_iter
        %% Tx

        raw_data=randi([0,1],rawbit_perframe,1);

        code_data=convenc(raw_data,trellis);

        % inter_data=tx_interleaver(code_data,Ndata/2,M_bit,Nframe); % interleaver
        inter_data=matintrlv(code_data,length(code_data)/32,32);

        X_data_1=qammod(inter_data,M_mod,'gray','InputType','bit','UnitAveragePower',true);

        X_data=X_data_1';

        x_tr=frame_generate(X_data,Nfft,Nsym,Ndata,Nvc,Nframe,Train_sym,X_pilot,Nps);  % generate a frame

        

        % x_tr= [Train_sym x_tr];

        %% Channel

        [y,h]=channel_project(x_tr,channeltype,fd,fmin,fmax);
        
%         delay_tap=[0 2 6 10];
% 
%         channel_coeff=sqrt(1/2)*(randn(1,length(delay_tap))+1j*randn(1,length(delay_tap)));
%         h1=zeros(1,delay_tap(end)+1);
%         h1(delay_tap+1)=channel_coeff;
% 
%         y=conv(x_tr,h1);
%         h=1;
% 
%         y_CFO=add_CFO(y,CFO,Nfft);

        [y_re,noise_var]=add_noise(y,SNR_dB(i),ExtraNoise,EndNoise);
        


        %% Rx


        delay_result(i,iter)=ExtraNoise;
        % delay_result=ExtraNoise;
        y_data=y_re(delay_result(i,iter)+1:end);


        % 

        h=1;
        % Y_est_data=frame_decompose(y_data,Nfft,Nsym,Ndata,Nvc,Nframe,channeltype,X_pilot,Nps,pilot_loc);

        [Y_est_data,H_DFT,noise_mean]=frame_decompose(y_data,Nfft,Nsym,Ndata,Nvc,Nframe,channeltype,X_pilot,Nps,pilot_loc,noise_var);

        X_est_data=qamdemod(Y_est_data',M_mod,'gray','OutputType','bit','UnitAveragePower',true);

        % X_est_deinter=matdeintrlv(X_est_data,length(X_est_data)/32,32);

        noise_var_1=10.^(-SNR_dB(i)/10)*mean(Y_est_data.*conj(Y_est_data));

        % noise_var=10.^(-SNR_dB(i)/10);

        X_est_data_llr=qamdemod(Y_est_data',M_mod,'OutputType','approxllr','UnitAveragePower',true);

        X_est_data_llr_deinter=matdeintrlv(X_est_data_llr,length(X_est_data_llr)/32,32);

        % X_est_decode_soft_1=channel_decode_soft(X_est_data_llr,noise_var,trellis);
        X_est_decode_soft_1=vitdec(X_est_data_llr_deinter,trellis,32,'cont','unquant');
       

 
        % X_est_reshape=reshape(X_est_data,bits_perframe,1);

        %X_est_deinter=rx_deinterleave(X_est_reshape,Ndata/2,M_bit,Nframe); % deinterleaver
        % X_est_deinter=matdeintrlv(X_est_data,length(X_est_data)/32,32);
        % X_est_deinter=X_est_reshape;

        % X_est_decode=channel_decode(X_est_data); % channel decode (2,1,7)

        X_est_decode_1=vitdec(X_est_data,trellis,32,'cont','hard');

      

        % X_est_decode=matdeintrlv(X_est_decode_1,length(X_est_decode_1)/32,32);

        % BER calculation
        total_biterrors(i)=total_biterrors(i)+sum(sum(X_est_decode_soft_1(33:end)~=raw_data(1:end-32)));
        % total_biterrors(i)=total_biterrors(i)+sum(sum(X_est_decode_soft_1~=raw_data));
        total_biterrors_hard(i)=total_biterrors_hard(i)+sum(sum(X_est_decode_1(33:end)~=raw_data(1:end-32)));
        % total_biterrors_hard(i)=total_biterrors_hard(i)+sum(sum(X_est_decode_1~=raw_data));

%         BER_ofdm(i)=total_biterrors(i)./(iter*rawbit_perframe);
%         BER_ofdm_hard(i)=total_biterrors_hard(i)./(iter*rawbit_perframe);


    end

end

%% plot the SNR-BER curve

BER_ofdm=total_biterrors./(max_iter*rawbit_perframe);
BER_ofdm_hard=total_biterrors_hard./(max_iter*rawbit_perframe);
%BER_ofdm_undecode=total_biterrors_undecode./(max_iter*bits_perframe);

% EbN0=SNR_dB-10*log10(M_bit);

% BER_ofdm(5)=398/(1e4*rawbit_perframe);
% BER_ofdm(6)=98/(1e4*rawbit_perframe);
% BER_ofdm(7)=37/(1e4*rawbit_perframe);
% SNR_dB=[0:5:30];max_iter*rawbit_perframe

figure(2);
semilogy(EbN0,BER_ofdm);
hold on;
semilogy(EbN0,BER_ofdm_hard);
xlabel('EbN0'),ylabel('BER')
legend('Soft','Hard');