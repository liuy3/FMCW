%% FBMC Modulation
clear all;
s = rng(211);            % Set RNG state for repeatability

%% System Parameters

numFFT = 1024;           % Number of FFT points
numGuards =1;%numFFT/8;         % Guard bands on both sides
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
numSymbols = 100;      % Simulation length in symbols
bitsPerSubCarrier = 2;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM

signal_power =0%numFFT*K*K;
SNRdB = 0:1:10;
SNR = 10.^(SNRdB/10);
No = signal_power./SNR;

CFO = 0:0.02:0.5;     
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
        dataSubCar(1:2:L) = real(modData);
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
    txSymb = numFFT*K*fftshift(ifft(X));
%txSymb = numFFT*K*ifft(X);
    % Transmitted signal is a sum of the delayed real, imag symbols
    symBuf = [symBuf(numFFT/2+1:end); complex(zeros(numFFT/2,1))];%%%%%%??????512+1:8192;512
    symBuf(KF+(1:KF)) = symBuf(KF+(1:KF)) + txSymb;

    % Compute power spectral density (PSD)
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

% Process symbol-wise
for t= 1:length(CFO)%%%%%%%%%%%%%%%%%%%
    BER = comm.ErrorRate;
    
    for symIdx = 1:numSymbols
        
        rxSig = txSigAll(:, symIdx);
       %rxSig=[complex(zeros(numFFT/2,1));rxSig(1:KF-numFFT/2)];
        % Add WGN
       % rxNsig = awgn(rxSig, snrdB, 'measured');%%%%%%%%%%%%%%%%%%why its wrong??
        noise = sqrt(No(11)/2)*complex(randn(KF,1),randn(KF,1));
                rxNsig = rxSig + noise;
        %Add CFO         
        yr =(exp(1i*2*pi*CFO(t)*(0:length(rxNsig)-1)/length(rxNsig))).*rxNsig.';
        
        %rxNsig=fftshift(rxNsig);
        %yr =(exp(1i*2*pi*CFO(t)*(0:length(rxNsig)-1)/length(rxNsig))).*rxNsig.';

%rxf = fft(yr.')/numFFT/K;
            % Perform FFT
        rxf = fft(fftshift(yr.'))/numFFT/K;
        % Perform FFT
       % rxf = fft(fftshift(yr.'))/K/numFFT;
        %rxf = fft(yr.')/K/numFFT;
        rxfAll(:, t)=rxf;
        
        % Matched filtering with prototype filter
        rxfmf = filter(Hk, 1, rxf);
        % Remove K-1 delay elements%%%%%%%%%%%%%%%%%%%%necessary??
        rxfmf = [rxfmf(K:end); zeros(K-1,1)];
        rxfmfAll(:, t)=rxfmf;
        % Remove guards
        rxfmfg = rxfmf(numGuards*K+1:end-numGuards*K);

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
    end
    
% Measure BER with appropriate delay
    %BER.ReceiveDelay = bitsPerSubCarrier*KL;
    ber = BER(inpData(:), rxBits(:));
    ERROR_FBMC(t)=ber(1)
    
end    

% Restore RNG state
rng(s);
figure(1)
semilogy(CFO, ERROR_FBMC,'b-*');
xlabel('CFO');
ylabel('BER');
title('FBMC CFO in QPSK Modulation'); 
%% Conclusion and Further Exploration
