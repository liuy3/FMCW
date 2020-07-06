%fmcw for RD-map CFAR
%last edited time :6/11/2018
clear all; close all; clc;
%parameter
rng(112);

S=8e12;  %8MHz/us
Tc=51e-6;  %51us
ns=2^8;
ts=Tc/ns;
t=0:ts:Tc-ts;
fs=1/ts;
f0=76e9;%76GHz
A=1;
c=3e8;
N=2^7;%128chirps
numRX=12;
f=linspace(-fs/2,fs/2,length(t));
r_max=fs/2*c/2/S;
v_max=pi*c/f0/Tc/4/pi;
SNRindB=-12;
SNR=10^(SNRindB/10);
E=1/2;

%object parameter
d=c/f0/2;
r=[20 10 30 10.5 10];%(r_max-0.3)*rand(1, 5)+0.3;%r max=40m% r res=1m ??r res=0.36m
v=[5 10 -15 10.5 v_max-0.5];%(v_max-0.15)*2*rand(1,100)-v_max+0.3;%v max=25m/s away:'+'%  v res=0.52m/s2*
theta=[pi/8 pi/4 pi/3 -pi/4 -pi/8 -pi/4 pi/4 pi/4 pi/8 pi/4];%ones(1, 100);%
r_hat=zeros(1, length(r));
v_hat=zeros(1, length(r));
theta_hat=zeros(1, length(r));

%tx signal
frame_tx=frame(0, 0, N, f0, S, t);
%rx signal
if_sig=zeros( N, length(t),numRX);
IF_sig=zeros( N, length(t),numRX);
R_FFT_SUM=zeros(N, length(t));
for j=1:numRX
    rx=zeros(1,length(t));
    for m=1:length(r)
        noise=sqrt(E/SNR)*randn(N, length(t));
        frame_rx=frame(2*v(m)*Tc/c, 2*r(m)/c+d*sin(theta(m))*(j-1)/c, N, f0, S, t);
        rx=rx+frame_rx+noise;
    end
    %IF signal
    if_sig( :, :,j)=frame_tx.*rx;
    %IF_sig(j)=zeros(N, length(t));
    for m=1:N
        IF_sig( m, :,j)=fftshift(fft(if_sig( m, :, j)))/fs;
    end
    R_FFT_SUM=R_FFT_SUM+IF_sig(:, :, j);
end

range=linspace(0, r_max, length(if_sig( 1, :, 1))/2);
%doppler FFT
RD_FFT=ones(N, length(if_sig( 1, :,1))/2, numRX);
RD_FFT_SUM=zeros(N, length(if_sig( 1, :,1))/2);
tic
for j=1:numRX
    for m=1:Tc/ts
        doppler_FFT(:,m, j)=fftshift(fft(IF_sig(:, m, j)));
    end
    RD_FFT(:,:,j)=IF_sig( :, length(if_sig( 1, :,1))/2+1:end,j)+doppler_FFT( :, length(if_sig( 1, :,1))/2+1:end, j);
    RD_FFT_SUM=RD_FFT_SUM+ abs(RD_FFT(:,:,j));
end
w=linspace(-pi, pi, N);
velocity=w*c/f0/Tc/4/pi;
RD_FFT_SUM1=R_FFT_SUM+doppler_FFT(:, :, 1);
RD_FFT_SUM_r=RD_FFT_SUM1(:, length(if_sig( 1, :,1))/2+1:end);
figure;imagesc(abs(RD_FFT(:, :,1)));title('R-D map from one RX');
%figure;imagesc(abs(RD_FFT_SUM_r));
figure;imagesc(abs(RD_FFT_SUM));title('R-D map from all RX');
figure;imagesc(abs(RD_FFT(:, :,1)));title('R-D map from one RX');
xticklabels = 0:5:r_max;
xticks = linspace(1,128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = linspace(pi, -pi, 10)*c/f0/Tc/4/pi;
yticks = linspace(1, 128, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Range');ylabel('velocity');
%{
%CFAR
refLength=8;
guardLength=3;
cfarWin=ones((refLength+guardLength)*2+1,1);
cfarWin(refLength+1:refLength+1+2*guardLength)=0;
cfarWin=cfarWin/sum(cfarWin);

Pfa=1e-2;
K=refLength*2*(Pfa^(-1/refLength/2)-1);
temp=abs(RD_FFT(:, :,1));%abs(RD_FFT_SUM);%abs(RD_FFT(:, :,1));%
temp=[temp(1:(refLength+guardLength), :);temp;temp(end-(refLength+guardLength)+1:end, :)];
noiseLevel=conv(temp(:),cfarWin,'same');
cfarThreshold=K*noiseLevel;
map1=(temp(:)>cfarThreshold);
map1=reshape(map1, [128+2*(refLength+guardLength) ,128]);
map1=map1((refLength+guardLength)+1:end-(refLength+guardLength), :);

Pfa=5e-2;
K=refLength*2*(Pfa^(-1/refLength/2)-1);
temp=abs(RD_FFT(:, :,1));%abs(RD_FFT_SUM);%abs(RD_FFT(:, :,1));%
temp=temp';
temp=[temp(1:(refLength+guardLength), :);temp;temp(end-(refLength+guardLength)+1:end, :)];
noiseLevel=conv(temp(:),cfarWin,'same');
cfarThreshold=K*noiseLevel;
map2=(temp(:)>cfarThreshold);
map2=reshape(map2, [128+2*(refLength+guardLength) ,128]);
map2=map2((refLength+guardLength)+1:end-(refLength+guardLength), :);
map2=map2';
cfarRDmap=map1.*map2;
%figure;imagesc(map1);
%figure;imagesc(map2);
figure;
imagesc(cfarRDmap);

%group peaking
temp1=cfarRDmap;
for n=1:128
    for j=1:128
        if(temp1(n, j)==1)indexvec= getvec(temp1, n, j);
            [I, I]=max(abs(RD_FFT(n, indexvec, 1)));
            for k=1:length(indexvec)
                temp1(n, j+k-1)=0;
            end
            temp1(n, j+I-1)=1;
        end
    end
end
temp2=cfarRDmap';
for n=1:128
    for j=1:128
        if(temp2(n, j)==1)indexvec= getvec(temp2, n, j);
            [I, I]=max(abs(RD_FFT( indexvec, j, 1)));
            for k=1:length(indexvec)
                temp2(n, j+k-1)=0;
            end
            temp2(n, j+I-1)=1;
        end
    end
end
PG_RDmap=temp1.*temp2';
%figure;imagesc(temp1);
%figure;imagesc(temp2');
figure;
imagesc(PG_RDmap);

%DOA
stepsize=pi/80;
[mm,ii]=sort(PG_RDmap(:));
for m=1:sum(PG_RDmap(:))
    for j=1:numRX
        temp= RD_FFT(:, :,j);
        RD_vec(:, j)=temp(:);
    end
    for j=1:numRX
        a(j)= RD_vec(ii(end+1-m), j);%steering vector
    end
    for k=1:pi/stepsize+1
        for j=1:numRX
            ww(j)=exp(-i*2*pi/c*f0*d*sin(-pi/2+(k-1)*stepsize)*(j-1));
        end
        p(k)=ww*a.';
    end
    r_hat(m)=floor(ii(end+1-m)/128)/128*r_max;
    v_hat(m)=mod(ii(end+1-m)-1, 128)/128*2*v_max-v_max;
    [I, I]=max(abs(p));
    theta_hat(m)=((I-1)*stepsize-pi/2)*180/pi;
    angle=linspace(-90, 90, pi/stepsize+1);
    %figure;plot(angle,abs(p));
end
sum(PG_RDmap(:))
%r=sort(r);r_hat=sort(r_hat);
toc
%}
%CFAR
refLength=8;
guardLength=3;
cfarWin=ones((refLength+guardLength)*2+1,1);
cfarWin(refLength+1:refLength+1+2*guardLength)=0;
cfarWin=cfarWin/sum(cfarWin);

Pfa=1.5e-1;
K=refLength*2*(Pfa^(-1/refLength/2)-1);
temp=abs(RD_FFT(:, :,1));%abs(RD_FFT_SUM);%
temp=[temp(1:(refLength+guardLength), :);temp;temp(end-(refLength+guardLength)+1:end, :)];
noiseLevel=conv(temp(:),cfarWin,'same');
cfarThreshold=K*noiseLevel;
map1=(temp(:)>cfarThreshold);
map1=reshape(map1, [128+2*(refLength+guardLength) ,128]);
map1=map1((refLength+guardLength)+1:end-(refLength+guardLength), :);

Pfa=2e-1;
K=refLength*2*(Pfa^(-1/refLength/2)-1);
temp=abs(RD_FFT(:, :,1));%abs(RD_FFT_SUM);%
temp=temp';
temp=[temp(1:(refLength+guardLength), :);temp;temp(end-(refLength+guardLength)+1:end, :)];
noiseLevel=conv(temp(:),cfarWin,'same');
cfarThreshold=K*noiseLevel;
map2=(temp(:)>cfarThreshold);
map2=reshape(map2, [128+2*(refLength+guardLength) ,128]);
map2=map2((refLength+guardLength)+1:end-(refLength+guardLength), :);
map2=map2';
figure;imagesc(map1);
figure;imagesc(map2);
cfarRDmap=map1.*map2;
figure;imagesc(cfarRDmap);


%group peaking
temp1=cfarRDmap;
for n=1:128
    for j=1:128
        if(temp1(n, j)==1)indexvec= getvec(temp1, n, j);
            [I, I]=max(abs(RD_FFT(n, indexvec, 1)));
            for k=1:length(indexvec)
                temp1(n, j+k-1)=0;
            end
            temp1(n, j+I-1)=1;
        end
    end
end
temp2=cfarRDmap';
for n=1:128
    for j=1:128
        if(temp2(n, j)==1)indexvec= getvec(temp2, n, j);
            [I, I]=max(abs(RD_FFT( indexvec, j, 1)));
            for k=1:length(indexvec)
                temp2(n, j+k-1)=0;
            end
            temp2(n, j+I-1)=1;
        end
    end
end
PG_RDmap=temp1.*temp2';
%figure;imagesc(temp1);
%figure;imagesc(temp2');
figure;imagesc(PG_RDmap);
xticklabels = 0:5:r_max;
xticks = linspace(1,128, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = linspace(pi, -pi, 10)*c/f0/Tc/4/pi;
yticks = linspace(1, 128, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Range');ylabel('velocity');
%DOA
stepsize=pi/80;
[mm,ii]=sort(PG_RDmap(:));
for m=1:sum(PG_RDmap(:))
    for j=1:numRX
        temp= RD_FFT(:, :,j);
        RD_vec(:, j)=temp(:);
    end
    for j=1:numRX
        a(j)= RD_vec(ii(end+1-m), j);%steering vector
    end
    for k=1:pi/stepsize+1
        for j=1:numRX
            ww(j)=exp(-i*2*pi/c*f0*d*sin(-pi/2+(k-1)*stepsize)*(j-1));
        end
        p(k)=ww*a.';
    end
    r_hat(m)=floor(ii(end+1-m)/128)/128*r_max;
    v_hat(m)=mod(ii(end+1-m)-1, 128)/128*2*v_max-v_max;
    [I, I]=max(abs(p));
    theta_hat(m)=((I-1)*stepsize-pi/2)*180/pi;
    angle=linspace(-90, 90, pi/stepsize+1);
    figure;plot(angle,abs(p));xlabel('angle');title(['Target ',num2str(m)]);
end
sum(PG_RDmap(:))