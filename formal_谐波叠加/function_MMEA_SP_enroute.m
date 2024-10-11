
function [f_result,C,Th] = function_MMEA_SP_enroute(N,K,fd,fmin,fmax,delay_dB)

%Define the Jakes power spectral density function
S = @(f, fd, scale) (scale./ (pi * fd * sqrt(1 - (f / fd).^2)));

% Define the integral function with the negative maximum Doppler frequency
f_integral = @(f, fd, C, scale, fmin) integral(@(f_prime) S(f_prime, fd, scale), fmin, f) - C;

alpha1=acos(fmin/fd);
alpha2=acos(fmax/fd);
alpha=abs(alpha2-alpha1);
alpha_re=alpha/pi;
scale=1/alpha_re;
f_result=zeros(K,N);

for n=1:N
    for k=1:K

        C=(n-1/2)/N;
        rotation_k=(k-((K+1)/2))/(K*N);
        C=C+rotation_k;  % The desired integral result

        % Define the function to find the root 
        root_func = @(f) f_integral(f, fd, C,scale,fmin);

        % Use fzero to find the root
        f_initial_guess = [fmin fmax]; % Initial guess for the root
        root = fzero(root_func, f_initial_guess);
        f_result(k,n)=root;

    end
end

% 参数c
C=zeros(K,N);
delay_linear=zeros(length(delay_dB));
delay_linear=10.^(delay_dB/10);
C(:,:)=1/sqrt(N);
for i=1:K
    C(i,:)=C(i,:)*sqrt(delay_linear(i))
end

% 参数theta
Th=rand(K,N)*2*pi;

end
