clc
clear 
close all

s = rng(211);            % Set RNG state for repeatability
CFO = 0:0.1:1;     
%% OFDM
Nfft = 1024;
Ng = Nfft/4;
%K = 64000*2;
M = Nfft*1000;
%
signal_power = (Nfft+Ng)*2;
SNRdB = 10; %%%%%%%%%%%%%%%%%%%%Ng
SNR = 10.^(SNRdB/10);
No = signal_power./SNR*Nfft/(Nfft+Ng);
%
count1 = 10;
for n = 1:length(CFO)%%%%%%%%%%%%%%%%%%%
    for k = 1:count1
        t = randi([0 1],1,2*M);
        for p = 1:M
            if (t(2*p-1)==1)&&(t(2*p)==1)
                o(p) = 1+j;
            elseif (t(2*p-1)==0)&&(t(2*p)==1)
                o(p) = -1+j; 
            elseif (t(2*p-1)==1)&&(t(2*p)==0)
                o(p) = 1-j; 
            else
                o(p) = -1-j;
            end
        end
        x = reshape(o,Nfft,M/Nfft).';%%S/P
        for m = 1:M/Nfft
            x_ifft(m,:) = Nfft*ifft(x(m,:));%%Nfft   ifft
        end
        x_cp = [x_ifft(:,Nfft-Ng+1:end) x_ifft];
        tx = reshape(x_cp.',1,M*(Nfft+Ng)/Nfft);%%P/S
        noise = sqrt(No/2)*(randn(1,M*(Nfft+Ng)/Nfft) + j*randn(1,M*(Nfft+Ng)/Nfft));
        rx = tx + noise;
        cfo=fftshift(exp(1i*2*pi*CFO(n)*(0:length(rx)-1)/length(rx)));
        yr =(cfo).*rx;%%%%%%%%%%%%%%%%%%%%%%
        X = reshape(yr,Nfft+Ng,M/Nfft).';%%%%%%%%%%%%%%%%%%%%%%%
        x_de_cp = X(:,Ng+1:end);
        for m = 1:M/Nfft
            x_fft(m,:) = fft(x_de_cp(m,:))/Nfft;%%%%%%%%%%%%%%%%%%%fft
        end
        xx = reshape(x_fft.',1,M);
        s_h_1 = sign(real(xx));
        s_h_2 = sign(imag(xx));
        temp = zeros(1,2*M);
        temp(1:2:end) = s_h_1;
        temp(2:2:end) = s_h_2;
        s_h = temp;
        detect = (s_h ~= 2*t-1);
        error(k) = mean(detect);
    end 
    ERROR_OFDM(n) = mean(error)
end
% Restore RNG state
rng(s);
figure(1)
semilogy(CFO, ERROR_OFDM,'b-*');
xlabel('normalized CFO');
ylabel('Bit Error Rate');

title('OFDM CFO in QPSK Modulation'); 