function [X_data,H_DFT,H_frame] = frame_decompose(y_data,Nfft,Nsym,Ndata,Nvc,Nframe,channeltype,X_pilot,Nps)

%% y_data  / rsceving data after syn

Ng=Nsym-Nfft; % CP length


if Nvc==0

    H_frame=zeros(Nframe,Nfft);

    noise_eq=1;

    X_mod=[];
    Xmod_r=0;

    kk1= 1:Nsym;
    kk2= 1:Nfft;
    kk3= 1:Ndata;
    kk4= Ndata/2+Nvc+1:Nfft;
    kk5= (Nvc~=0)+[1:Ndata/2];


    % outcp,FFT,Equalization

    for k=1:Nframe

        y_handle=y_data(kk1);
        y_outcp=outcp(y_handle,Nfft,Ng); % Remove CP

        Y(kk2)=fft(y_outcp,Nfft)/sqrt(Nfft);

        Y_shift=[Y(kk4),Y(kk5)];


        if channeltype==0
            Xmod_r=Y_shift;
            H_DFT=1;

        elseif (channeltype==1)||(channeltype==2)||(channeltype==3)||(channeltype==4)
            % Xmod_r=Y_shift./H_shift;

            %% LS/DFT channel estimation

            H_est=LS_test(Y_shift,pilot_loc,X_pilot,Nfft,Nvc);

            channel_length=201;
            H_DFT=LS_DFT(H_est,channel_length);

            H_frame(k,:)=H_DFT;


            % Xmod_r=Y_shift./H_DFT([[(Nvc~=0)+[1:Ndata/2]],[Ndata/2+Nvc+1:Nfft]]);
            Xmod_r=Y_shift./H_DFT;

            % conj(H_tf)./(H_tf.*conj(H_tf)+noisepower);
            % Xmod_r=Y_shift./H_shift;
        end
% 
        % figure(1);
        % subplot(2,1,1)
        % plot(Y_shift,'.','MarkerSize',5);
        % axis([-1.5 1.5 -1.5 1.5]);
        % subplot(2,1,2)
        % plot(Xmod_r,'.','MarkerSize',5);
        % axis([-1.5 1.5 -1.5 1.5]);



        X_mod=[X_mod Xmod_r];

        kk1=kk1+Nsym;
        kk2=kk2+Nfft;
        kk3=kk3+Ndata;
        kk4=kk4+Nfft;
        kk5=kk5+Nfft;

    end

    X_data=extract_data(X_mod,Ndata,Nframe,Nps); % extract the data from the pilot


elseif Nvc>0

    X_mod=[];

    H_frame = zeros(Nframe/2,Nfft);

 
    kk1=1:Nsym
    kk2=1:Nfft;
    kk3=1:Ndata;
    kk4= Ndata/2+Nvc+1:Nfft;
    kk5= (Nvc~=0)+[1:Ndata/2];

    for k=1:Nframe


        y_handle = y_data((k-1)*Nsym+1:k*Nsym);
        y_outcp = outcp(y_handle,Nfft,Ng); % Remove CP

        Y = fft(y_outcp,Nfft)/sqrt(Nfft);

        % Y_shift=fftshift(Y(kk2));

        if (mod(k,2)==1)

            channel_length=201;

            [H_final,H_DFT]=LS_vc(Y,X_pilot,Nfft,Nvc,channel_length);


        elseif (mod(k,2)==0)

            H_frame(k/2,:) = H_DFT;

            % 用前一帧的信道估计做均衡

            H_shift(kk3) = [H_DFT(kk4) H_DFT(kk5)];

            Y_shift = [Y(kk4) Y(kk5)];

            Xmod_r = Y_shift./H_shift(kk3);



            figure(1);
            subplot(2,1,1)
            plot(Y,'.','MarkerSize',5);
            axis([-1.5 1.5 -1.5 1.5]);
            subplot(2,1,2)
            plot(Xmod_r,'.','MarkerSize',5);
            axis([-1.5 1.5 -1.5 1.5]);


            %H_DFT_dB=10*log10(abs(H_DFT.*conj(H_DFT)));
        
            X_mod=[X_mod Xmod_r]
            % 
            % kk1=kk1+Nsym;
            % kk2=kk2+Nfft;
            % kk3=kk3+Ndata;
            % kk4=kk4+Nfft;
            % kk5=kk5+Nfft;

        end


    end

    % noise_mean=mean(noise_eq);

    X_data=X_mod;



end


end