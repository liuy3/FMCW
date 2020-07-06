clc
clear 
close all

%s = rng(211);            % Set RNG state for repeatability
CFO = 0:0.05:0.2;     
%% OFDM
Nfft = 1024;
Ng = Nfft/4;
%K = 64000*2;
M = Nfft*1000;
%
signal_power =(Nfft+Ng)*2;
SNRdB = 10; %%%%%%%%%%%%%%%%%%%%Ng
SNR = 10.^(SNRdB/10);
No = signal_power./SNR*Nfft/(Nfft+Ng);
%
count1 =5;
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
        x = reshape(o,Nfft,M/Nfft).';
        for m = 1:M/Nfft
            x_ifft(m,:) = Nfft*ifft(x(m,:));%%Nfft   ifft
        end
        x_cp = [x_ifft(:,Nfft-Ng+1:end) x_ifft];
        tx = reshape(x_cp.',1,M*(Nfft+Ng)/Nfft);
        noise = sqrt(No/2)*(randn(1,M*(Nfft+Ng)/Nfft) + j*randn(1,M*(Nfft+Ng)/Nfft));
        rx = tx + noise;
        yr =(exp(1i*2*pi*CFO(n)*(0:length(rx)-1)/length(rx))).*rx;%%%%%%%%%%%%%%%%%%%%%%
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
%rng(s);

%% FBMC vs. OFDM Modulation


%s = rng(211);            % Set RNG state for repeatability

%% System Parameters

numFFT = 1024;           % Number of FFT points
numGuards = 1;%numFFT/8;         % Guard bands on both sides
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
numSymbols = 1000;        % Simulation length in symbols
bitsPerSubCarrier = 2;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 10;              % SNR in dB

signal_power = numFFT*K*K;
SNRdB = 0:1:10;
SNR = 10.^(SNRdB/10);
No = signal_power./SNR;

%CFO = 0:0.2:8;     
%% Filter Bank Multi-Carrier Modulation
% Prototype filter
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end
% Build symmetric filter
Hk = [fliplr(HkOneSided) 1 HkOneSided];

% QAM symbol mapper
qamMapper = comm.RectangularQAMModulator(...
    'ModulationOrder', 2^bitsPerSubCarrier, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%   Initialize arrays
L = numFFT-2*numGuards;  % Number of complex symbols per OFDM symbol
KF = K*numFFT;
KL = K*L;
dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);

sumFBMCSpec = zeros(KF*2, 1);
sumOFDMSpec = zeros(numFFT*2, 1);

numBits = bitsPerSubCarrier*L/2;    % account for oversampling by 2
inpData = zeros(numBits, numSymbols); 
rxBits = zeros(numBits, numSymbols);
txSigAll = complex(zeros(KF, numSymbols));
symBuf = complex(zeros(2*KF, 1));

% Loop over symbols
for symIdx = 1:numSymbols
    
    % Generate mapped symbol data
    inpData(:, symIdx) = randi([0 1], numBits, 1);
    modData = qamMapper(inpData(:, symIdx));
    
    % OQAM Modulator: alternate real and imaginary parts
    if rem(symIdx,2)==1     % Odd symbols
        dataSubCar(1:2:L) = real(modData);%%%%%%%%%%%%%%%mismatch??
        dataSubCar(2:2:L) = 1i*imag(modData);
    else                    % Even symbols
        dataSubCar(1:2:L) = 1i*imag(modData);
        dataSubCar(2:2:L) = real(modData);
    end

    % Upsample by K, pad with guards, and filter with the prototype filter
    dataSubCarUp(1:K:end) = dataSubCar;
    dataBitsUpPad = [zeros(numGuards*K,1); dataSubCarUp; zeros(numGuards*K,1)];
    X1 = filter(Hk, 1, dataBitsUpPad);
    % Remove 1/2 filter length delay
    X = [X1(K:end); zeros(K-1,1)];%%%%%%%%%%%%%??????????

    % Compute IFFT of length KF for the transmitted symbol
    txSymb = K*numFFT*(fftshift(ifft(X)));

    % Transmitted signal is a sum of the delayed real, imag symbols
    symBuf = [symBuf(numFFT/2+1:end); complex(zeros(numFFT/2,1))];%%%%%%??????512+1:8192;512
    symBuf(KF+(1:KF)) = symBuf(KF+(1:KF)) + txSymb;

    currSym = complex(symBuf(1:KF));

    % Store transmitted signals for all symbols
    txSigAll(:,symIdx) = currSym;
end

for symIdx = numSymbols+1:numSymbols+2*K
    symBuf = [symBuf(numFFT/2+1:end); complex(zeros(numFFT/2,1))];
    symBuf(KF+(1:KF)) = symBuf(KF+(1:KF)) + complex(zeros(KF,1));
    txSigAll(:,symIdx) = complex(symBuf(1:KF));
end

for symIdx = 1:numSymbols
    txSigAll(:,symIdx) = txSigAll(:,symIdx+2*K);
end

%% FBMC Receiver with No Channel


% QAM demodulator
qamDemod = comm.RectangularQAMDemodulator(...
    'ModulationOrder', 2^bitsPerSubCarrier, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

ERROR_FBMC=zeros(1, length(CFO));
count1 = 5;
% Process symbol-wise
for t= 1:length(CFO)%%%%%%%%%%%%%%%%%%%
    for k=1:count1
        BER = comm.ErrorRate;    
        for symIdx = 1:numSymbols

            rxSig = txSigAll(:, symIdx);

            % Add WGN
           % rxNsig = awgn(rxSig, snrdB, 'measured');%%%%%%%%%%%%%%%%%%why its wrong??
            noise = sqrt(No(11)/2)*complex(randn(KF,1),randn(KF,1));
                    rxNsig = rxSig + noise;
            %Add CFO         
            %rxNsig=fftshift(rxNsig);
            cfo=fftshift((exp(1i*2*pi*CFO(t)*(0:length(rxNsig)-1)/length(rxNsig))));
            yr =cfo.*rxNsig.';


            % Perform FFT
            rxf = fft((fftshift(yr.')))/numFFT/K;

            % Matched filtering with prototype filter
            rxfmf = filter(Hk, 1, rxf);
            % Remove K-1 delay elements%%%%%%%%%%%%%%%%%%%%necessary??
            rxfmf = [rxfmf(K:end); zeros(K-1,1)];
            % Remove guards
            rxfmfg = rxfmf(numGuards*K+1:end-numGuards*K);

            %if(CFO>0.5)rxfmfg=-rxfmfg;end
            % OQAM post-processing
            %  Downsample by 2K, extract real and imaginary parts
            if rem(symIdx, 2)
                % Imaginary part is K samples after real one
                r1 = real(rxfmfg(1:2*K:end));
                r2 = imag(rxfmfg(K+1:2*K:end));
                rcomb = complex(r1, r2);
            else
                % Real part is K samples after imaginary one
                r1 = imag(rxfmfg(1:2*K:end));
                r2 = real(rxfmfg(K+1:2*K:end));
                rcomb = complex(r2, r1);
            end
            %  Normalize by the upsampling factor
            rcomb = (1/K)*rcomb;

            % Demapper: Perform hard decision
            rxBits(:, symIdx) = qamDemod(rcomb);    
            %if(CFO(t)>0.5)
           %     rxBits(:, symIdx)=1-rxBits(:, symIdx);
           % end
        end

    % Measure BER with appropriate delay
        %BER.ReceiveDelay = bitsPerSubCarrier*KL;
        ber = BER(inpData(:), rxBits(:));
        ERROR_FBMC(t)=ERROR_FBMC(t)+ber(1);
    end
     ERROR_FBMC(t)= ERROR_FBMC(t)/count1
end    

%rng(s);


figure(1)
semilogy(CFO, ERROR_OFDM,'b-*',CFO, ERROR_FBMC,'r-*');
xlabel('normalized CFO');
ylabel('Bit Error Rate');
legend('OFDM','FBMC')
title('OFDM vs. FBMC CFO QPSK Modulation'); 