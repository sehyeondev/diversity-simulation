clear
load SC_SER.mat
load MRC_SER.mat
load MRC_SC_SER.mat
load SISO_SER.mat
load A21_SER.mat

SNR_dB = 0:2:20;            % range of SNR dB values
SNR = db2pow(SNR_dB);                   % dB to power

figure;%selection combining

semilogy(SNR_dB,SISO_SER(:,:),"s-","LineWidth",2);hold on;

% semilogy(SNR_dB,SC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,SC_SER(2,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,SC_SER(3,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,SC_SER(4,:),"s-","LineWidth",2);hold on;

% semilogy(SNR_dB,MRC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_dB,MRC_SER(2,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,MRC_SER(3,:),"s-","LineWidth",2);hold on;
% semilogy(SNR_dB,MRC_SER(4,:),"s-","LineWidth",2);hold on;

semilogy(SNR_dB,MRC_SC_SER(1,:),"s-","LineWidth",2);hold on;

semilogy(SNR_dB,A21_SER(1,:),"s-","LineWidth",2);hold on;

legend("SISO", "SC | L=2","MRC | L=2","MRC+SC | L=2","2x1 Alamouti");
