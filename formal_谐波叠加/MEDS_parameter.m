function [F1,F2,C1,C2,Th1,Th2] = MEDS_parameter(N_i,fd,delay_dB)

% N_i 生成谐波数目的一列数据
% fd  最大多普勒频偏

delay_linear=10.^(delay_dB/10);
max_num=max(N_i)+1;
F1=zeros(length(N_i),max_num);
F2=zeros(length(N_i),max_num);
C1=zeros(length(N_i),max_num);
C2=zeros(length(N_i),max_num);
Th1=zeros(length(N_i),max_num);
Th2=zeros(length(N_i),max_num);

for i=1:length(N_i)
n=(1:N_i(i));
F1(i,[1:N_i(i)])=sin(pi/(2*N_i(i))*(n-1/2));
C1(i,[1:N_i(i)])=sqrt(delay_linear(i))*sqrt(1/N_i(i))*ones(1,N_i(i));
Th1(i,[1:N_i(i)])=rand(1,N_i(i))*2*pi;

n1=(1:N_i(i)+1);
F2(i,[1:N_i(i)+1])=sin(pi/(2*(N_i(i)+1))*(n1-1/2));
C2(i,[1:N_i(i)+1])=sqrt(delay_linear(i))*sqrt(1/(N_i(i)+1))*ones(1,N_i(i)+1);
Th2(i,[1:N_i(i)+1])=rand(1,N_i(i)+1)*2*pi;

end

F1=F1*fd;
F2=F2*fd;

end
