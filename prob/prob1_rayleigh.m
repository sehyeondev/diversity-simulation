close all
clear all
clc

j = sqrt(-1);                           % define complex number
N = 1e6;                               % define number of symbols
s = [-1+j -1-j +1+j +1-j];              % define BPSK symbols
s = s/norm(s)*sqrt(length(s));          % normalize the energy
SNR_dB = 0:2:20;             % range of SNR dB values
SNR = db2pow(SNR_dB);                   % convert dB to power
sigma = 1./sqrt(2*SNR);                 % define noise standard deviation
error = zeros(1,length(SNR));           % initialize the error
error_awgn = zeros(1,length(SNR));           % initialize the error

% Ps_theory = 2.*(qfunc(sqrt(SNR))) + (qfunc(sqrt(SNR))).^2;
Ps_rayleigh = zeros(length(SNR), 1);
Ps_gaussian = zeros(length(SNR), 1);

a = 1;
for i = 1:length(SNR)    
    txd = randsrc(1,N,s);                   % generate QPSK symbols
    txd = repmat(txd,[a,1]);

    noise = randn(a,N) + j*randn(a,N);      % generate noise
    ch = (randn(a,N)+j*randn(a,N))/sqrt(2); % generate complex-gaussian channel with unit power

    % rayleigh fading and AWGN channel porcess
    rxd = ch.*txd + noise.*sigma(i);
    rxd_awgn = txd + noise.*sigma(i);
    
    % ML decoding (find minimum distance)
    ML =  [   (rxd-(s(1)*ones(a,length(rxd))).*ch).^2;
              (rxd-(s(2)*ones(a,length(rxd))).*ch).^2;
              (rxd-(s(3)*ones(a,length(rxd))).*ch).^2;
              (rxd-(s(4)*ones(a,length(rxd))).*ch).^2];           
    
    [ ~, index] = min(ML,[],1);               
    index(index == 1) = s(1);
    index(index == 2) = s(2);  
    index(index == 3) = s(3);
    index(index == 4) = s(4);     
    
    ML_awgn =  [   (rxd_awgn-(s(1)*ones(a,length(rxd)))).^2;
              (rxd_awgn-(s(2)*ones(a,length(rxd_awgn)))).^2;
              (rxd_awgn-(s(3)*ones(a,length(rxd_awgn)))).^2;
              (rxd_awgn-(s(4)*ones(a,length(rxd_awgn)))).^2];           
    [ ~, index_awgn] = min(ML_awgn,[],1);               
    index_awgn(index_awgn == 1) = s(1);
    index_awgn(index_awgn == 2) = s(2);  
    index_awgn(index_awgn == 3) = s(3);
    index_awgn(index_awgn == 4) = s(4);     
    
    % count the error
    error(i) = size(find(index-txd),2);

    error_awgn(i) = size(find(index_awgn-txd),2);

    % theory
    gamma_i = SNR(i);
    Ps_gaussian(i) = 2.*(qfunc(sqrt(gamma_i))) + (qfunc(sqrt(gamma_i))).^2;
    Ps_rayleigh(i) = 1 / (1+gamma_i) - (1/4)*(1 / (2*gamma_i+1));
    % Ps_rayleigh(i) = 0.5 * (1 - sqrt(gamma_i / (2 + gamma_i)));
    
end
error = error/N;
error_awgn = error_awgn/N;

M = 4;
b = sqrt((2*(sin(pi/M))^2)^2.*SNR/2)./sqrt(1 + (2*(sin(pi/M))^2)^2.*SNR/2);
Pe = ((M-1)/M).*(1-b.*((M-1)/M)^-1/pi.*(pi/2+atan(b*cot(pi/M))));    

figure(1)
% semilogy(SNR_dB,Ps_theory,'k-o');
semilogy(SNR_dB, Ps_rayleigh, 'ko-', 'LineWidth', 0.5); % Rayleigh
hold on;
semilogy(SNR_dB,error,'r-x');
hold on;
semilogy(SNR_dB, Ps_gaussian, 'bo-', 'LineWidth', 0.5); % AWGN
hold on;
semilogy(SNR_dB,error_awgn,'r-x');
hold on;

grid on;

axis([min(SNR_dB) max(SNR_dB) 10^-8 10^0])
legend('Rayleigh Theory','Rayleigh Simulation', "AWGN theory", "AWGN simulation");
xlabel('SNR [dB]');
ylabel('SER');

SISO_SER=error;
save("SISO_SER","SISO_SER");


