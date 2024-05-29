clear
load SC_SER.mat
load MRC_SER.mat

SNR_dB = 0:2:20;            % range of SNR dB values
SNR = db2pow(SNR_dB);                   % dB to power

figure;%selection combining

% semilogy(SNR_dB,SC_SER(1,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,SC_SER(2,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,SC_SER(3,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,SC_SER(4,:),"s-","LineWidth",2);hold on;

semilogy(SNR_dB,MRC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,MRC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,MRC_SER(3,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,MRC_SER(4,:),"s-","LineWidth",2);hold on;
legend("M=1","M=2","M=3","M=4","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("Maximal Ratio Combining");
ylim([10^-5 1]);grid on;axis square;
