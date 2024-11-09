function [codedata]=LDPC_code(indata)

codedata = zeros(length(indata)*2,1);

rate = 1/2;
len = 1944;

[cfgLDPCEnc,decodercfg] = generateConfigLDPC(rate,len);

loop_num = length(indata)/(len/2);

for k=1:loop_num

    codedata((k-1)*len+1:k*len) = ldpcEncode(indata((k-1)*(len/2)+1:k*(len/2)),cfgLDPCEnc);
 
end

end