function [mu_k,tdomain] = SOS_MEDS(F1,F2,C1,C2,Th1,Th2,K_rice,f_los,T_sim,fs)

% 莱斯因子输入为线性值

K_taxi=K_rice;
Ts=1/fs;
Ns=round(T_sim/Ts); % 生成点数
tdomain=[0:Ts:(Ns-1)*Ts];

sizeF=size(F1);
mu_k=zeros(sizeF(1),Ns)

for i=0:Ns-1
mu_k(:,i+1)=sum(C1.*cos(2*pi*F1*(i-1)*Ts+Th1),2)+1j*sum(C2.*cos(2*pi*F2*(i-1)*Ts+Th2),2);
end


for l=0:Ns-1
mu_k(1,l+1)=sqrt(1/(K_taxi+1))*mu_k(1,l+1)+sqrt(K_taxi/(K_taxi+1))*exp(1j*(2*pi*f_los*l*Ts));
end

mu_k([2:end],:)=sqrt(1/(K_taxi+1)).*mu_k([2:end],:);

end