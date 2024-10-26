function [decode_data] = LDPC_decode(code_llr,last_len)

%% code_lrr: 编码数据的LLR
%% last_len: 最后一个324 LDPC编码中数据的长度

[cfgLDPCEnc,decodercfg] = generateConfigLDPC(1/2);

k_bits = cfgLDPCEnc.NumInformationBits;

maxnumiter = 50; % Number of iterations for LDPC decoder

if last_len==0
    
    decode_data = zeros(length(code_llr)/2,1);

    for i=1:length(code_llr)/(k_bits*2)
 
    decode_data((i-1)*k_bits+1:k_bits*i) = ldpcDecode(code_llr((i-1)*(2*k_bits)+1:i*(2*k_bits)),decodercfg,maxnumiter);

    end

elseif last_len>0
    
    decode_data = zeros(length(code_llr)/2,1);

    % enc_decode = zeros(length(code_llr)/2+(k_bits-last_len),1);

    % 先解码前面完整的
    for i=1:(length(code_llr)-2*last_len)/(2*k_bits)

    decode_data((i-1)*k_bits+1:k_bits*i) = ldpcDecode(code_llr((i-1)*(2*k_bits)+1:i*(2*k_bits)),decodercfg,maxnumiter);
    
    end

    rest_llr = code_llr(end-2*last_len+1:end);
    rest_llr_decode = [rest_llr.' zeros(1,2*(k_bits-last_len))].';  % 补0
    
    rest_decode = ldpcDecode(rest_llr_decode,decodercfg,maxnumiter);  % 补0后的解码
    
    decode_data (end-last_len+1:end) = rest_decode (1:last_len);  % 最后的比特

end

end

