function decode_data=LDPC_decode(inllr,maxnumiter)

rate = 1/2;
len = 1944;

% [cfgLDPCEnc,decodercfg] = generateConfigLDPC(rate,len);

P = [57 -1 -1 -1 50 -1 11 -1 50 -1 79 -1  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    3 -1 28 -1  0 -1 -1 -1 55  7 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1
    30 -1 -1 -1 24 37 -1 -1 56 14 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1
    62 53 -1 -1 53 -1 -1  3 35 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1
    40 -1 -1 20 66 -1 -1 22 28 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1
    0 -1 -1 -1  8 -1 42 -1 50 -1 -1  8 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1
    69 79 79 -1 -1 -1 56 -1 52 -1 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1
    65 -1 -1 -1 38 57 -1 -1 72 -1 27 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1
    64 -1 -1 -1 14 52 -1 -1 30 -1 -1 32 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1
    -1 45 -1 70  0 -1 -1 -1 77  9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1
    2 56 -1 57 35 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0
    24 -1 61 -1 60 -1 -1 27 51 -1 -1 16  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0
    ];

blockSize  = 81;

pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);

decodercfg = ldpcDecoderConfig(pcmatrix,'offset-min-sum');

loop_num = length(inllr)/len;

decode_data = zeros((len/2)*loop_num,1);

for k=1:loop_num

decode_data((k-1)*(len/2)+1:k*(len/2)) = ldpcDecode(inllr((k-1)*len+1:k*len),decodercfg,maxnumiter);

end