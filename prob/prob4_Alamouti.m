clear
load SC_SER.mat
load MRC_SER.mat
load A21_SER.mat

SNR_dB = 0:2:20;            % range of SNR dB values
SNR = db2pow(SNR_dB);                   % dB to power


%a21 vs. sc vs. mrc
figure;
semilogy(SNR_dB,SC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,MRC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,A21_SER(1,:),"s-","LineWidth",2);hold on;
legend("SC | L=2","MRC | L=2","2x1 Alamouti");
xlabel("SNR(dB)");ylabel("SER");title("MRC with 2 Antennas vs. 2x1 Alamouti");
ylim([10^-5 1]);grid on;axis square;