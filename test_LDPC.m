clear all;


[cfgLDPCEnc,decodercfg] = generateConfigLDPC(1/2);

k_bits = cfgLDPCEnc.NumInformationBits;

maxnumiter = 50; % Number of iterations for LDPC decoder

raw_data=randi([0 1],10240,1);

[code_data,last_len]=LDPC_code(raw_data);

M=4;
modSignal = qammod(code_data,M,'InputType','bit','UnitAveragePower',true);

snrdB = 10; % SNR in dB

receivedSignal = awgn(modSignal,snrdB);

% Computing the effective noise variance for given snr
noiseVar = 1/10^(snrdB/10);

% Computing soft estimate i.e., llr
code_llr = qamdemod(receivedSignal,M,'OutputType','llr','UnitAveragePower',true,'NoiseVariance',noiseVar);


decode_data = LDPC_decode(code_llr,last_len);

final=sum(decode_data==raw_data);
