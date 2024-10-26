function [code_data,last_len] = LDPC_code(raw_data)

%% 1/2 码率 LDPC 编码  324/648

[cfgLDPCEnc,decodercfg] = generateConfigLDPC(1/2);

k_bits = cfgLDPCEnc.NumInformationBits;

% encData=0;

raw_len = length(raw_data);  

last_len = mod(raw_len,k_bits);  % 计算是否需要补0

if last_len>0

    full_bits = [raw_data.' zeros(1,k_bits-last_len)];

    full_bits = full_bits.';

    encData = zeros(length(full_bits)*2,1);
    
    for i=1:length(full_bits)/k_bits

    encData((i-1)*(2*k_bits)+1:i*(2*k_bits)) = ldpcEncode(full_bits((i-1)*k_bits+1:i*k_bits),cfgLDPCEnc);

    end

    code_data = encData(1:end-(k_bits-last_len)*2);
    
elseif last_len==0
    
    encData = zeros(length(raw_data)*2,1);

    for i=1:length(raw_data)/k_bits

    encData((i-1)*(2*k_bits)+1:i*(2*k_bits)) = ldpcEncode(raw_data((i-1)*k_bits+1:i*k_bits),cfgLDPCEnc);

    end

    code_data=encData;

end


end