clear;close all;
snr_db=0:2:20;%snr values in db
N02=1./(10.^(snr_db./10));%snr values
sigmas=sqrt(N02);%std dev of awgn
fading=sqrt(1);%fading std dev

antennas=[1 2 3 4];%antenna count
S=length(snr_db);%constants
A=length(antennas);

snr_ber=zeros(1,length(snr_db));%snr vs ber vector
snr_ser=zeros(1,length(snr_db));%snr vs fer vector

M=4;%modulation alphabet size
m=1:M;%symbol index

spf=1e6;%symbols per frame
bps=log2(M);%bits per symbol
bpf=bps*spf;%bits per frame

bb=de2bi((m-1)')';%bit book
sym=(exp(1j*2*pi*(m-1)/M))';%complex symbol book

for a=1:A
    for s=1:S
        %variance of awgn noise
        sigma=sigmas(s);
            
        %transmitter
        bits=randi([0 1],[1 bpf]);%generate information sequence
        bits=reshape(bits,[bps spf]);%divide sequence into log2(M)-bits
        sym_ind=(bi2de(bits')+1)';%compute symbol indexes
        x=zeros(1,spf);%transmit signal
        for i=1:spf
            %map to symbols
            x(1,i)=sym(sym_ind(:,i),:);
        end
        x=repmat(x,[a,1]);%repeat for all antennas
        
        %channel
        %complex gaussian noise
        n_ip=normrnd(0,sigma,[a,spf]);
        n_q=normrnd(0,sigma,[a,spf]);
        n=(n_ip+1j*n_q)./sqrt(2);
        %rayleigh fading coeffs
        h_ip=normrnd(0,fading,[a,spf]);
        h_q=normrnd(0,fading,[a,spf]);
        h=(h_ip+1j*h_q)./sqrt(2);
        %impose channel
        r=h.*x+n;
        
        %receiver
        %combining
        r_tmp=zeros(1,spf);
        h_tmp=zeros(1,spf);
        [~, max_ind]=max(h.*conj(h),[],1);
        for i=1:spf
            r_tmp(1,i)=r(max_ind(1,i),i);
            h_tmp(1,i)=h(max_ind(1,i),i);
        end
        r_comb=r_tmp.*exp(-1j*angle(h_tmp));
        %detection
        r_comb=repmat(r_comb,[M 1]);%stack received signal, divide by channel impulse
        distance=abs(sym-r_comb);%calculate distance between received symbols and mod. symbols
        [~,det_sym_ind]=min(distance,[],1);%minimum distance for mapping symbols back
        %errors
        sym_err=nnz(sym_ind-det_sym_ind);%symbol errors
        bit_err=nnz(bits-bb(:,det_sym_ind));%calculate bit error
            
        [a snr_db(s) sym_err]

        snr_ser(a,s) = sym_err/spf;
    end%end snr loop
end%end antenna loop
SC_SER=snr_ser;
save("SC_SER","SC_SER");