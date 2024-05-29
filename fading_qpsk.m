close all
clear all
clc

j = sqrt(-1);                           % define complex number
N = 1e+6;                               % define number of symbols
s = [-1+j -1-j +1+j +1-j];              % define BPSK symbols
s = s/norm(s)*sqrt(length(s));          % normalize the energy
txd = randsrc(1,N,s);                   % generate QPSK symbols
SNR_dB = linspace(0,30,11);             % range of SNR dB values
SNR = db2pow(SNR_dB);                   % convert dB to power
sigma = 1./sqrt(2*SNR);                 % define noise standard deviation
noise = randn(1,N) + j*randn(1,N);      % generate noise
ch = (randn(1,N)+j*randn(1,N))/sqrt(2); % generate complex-gaussian channel with unit power
error = zeros(1,length(SNR));           % initialize the error

for i = 1:length(SNR)    

    % rayleigh fading and AWGN channel porcess
    rxd = ch.*txd + noise*sigma(i);
    
    % ML decoding (find minimum distance)
    ML =  [   (rxd-(s(1)*ones(1,length(rxd))).*ch).^2;
              (rxd-(s(2)*ones(1,length(rxd))).*ch).^2;
              (rxd-(s(3)*ones(1,length(rxd))).*ch).^2;
              (rxd-(s(4)*ones(1,length(rxd))).*ch).^2];           
    [ ~, index] = min(ML);               
    index(index == 1) = s(1);
    index(index == 2) = s(2);  
    index(index == 3) = s(3);
    index(index == 4) = s(4);     
    
    % count the error
    error(i) = size(find(index-txd),2);
    
end
error = error/N;

M = 4;
b = sqrt((2*(sin(pi/M))^2)^2.*SNR/2)./sqrt(1 + (2*(sin(pi/M))^2)^2.*SNR/2);
Pe = ((M-1)/M).*(1-b.*((M-1)/M)^-1/pi.*(pi/2+atan(b*cot(pi/M))));    

figure(1)
semilogy(SNR_dB,Pe,'k-o');
hold on
semilogy(SNR_dB,error,'r-x');
grid on
legend('QPSK Theory','QPSK Simulation');
axis([min(SNR_dB) max(SNR_dB) 10^-6 10^0])
xlabel('SNR [dB]');
ylabel('SER');
