
%% OFDM-IM simulation with ML, GD and LLR detectors and imperfect CSI
% Performance metics: SEP symbol error probability [1], BER in [2]
% Matlab version 2015b, also working well on 2019a.
% Note that MCIK-OFDM is another name of OFDM-IM see [1], [2].

%% Author information
% Thien Van Luong, Queen's University Belfast, UK, now with
% University of Southampton, UK.
% Email: tluong01@qub.ac.uk or thien.luong@soton.ac.uk.
% Personal page: https://tvluong.wordpress.com

%% References
% [1] T. V. Luong and Y. Ko, “A tight bound on BER of MCIK-OFDM with
% greedy detection and imperfect CSI,” IEEE Commun. Lett., vol. 21,
% no. 12, pp. 2594 – 2597, Dec. 2017.
% [2] T. V. Luong and Y. Ko, “Impact of CSI uncertainty on MCIK-OFDM:
% tight, closed-form symbol error probability analysis,” IEEE Trans. Veh.
% Technol., vol. 67, no. 2, pp. 1272 – 1279, Feb. 2018.

%% ==============================OFDM-IM/MCIK-OFDM=================================
clc

%% System parameters
M=4; % M-ary modulation size
N=4; % number of sub-carriers
K=1; % number of active sub-carriers

% imperfect CSI setting
var = 0.05; % fixed imperfect CSI variance, see [1], [2]
mmse = 1; % variable CSI
CSI=1; % 1 perfect CSI, 2 fixed CSI error variance, 3 MMSE variable CSI error variance

Detect_method =1; % to select 1 ML, 2 LLR, 3 Greedy GD detector
LLR = 1; % select one of two types of LLR detector

ro=0;
Mary=1; % 1 PSK, 2 QAM
if(M==8)
    QAM = (5*M-4)./6; % QAM power scale factor
else
    QAM = (2/3)*(M-1);
end

tic
%% ======================= Misc Parameters ================================
iter = 4;  % # Iterations
nSymPerFrame = 1e4; % Number of symbol per frame(1 OFDM symbol)
EbN0dB = 0:5:40;
EsN0 = 10.^(EbN0dB/10);
sigma = sqrt(1./EsN0); % additive noise variance

PwrSC = N/K; % Average Tx power per active sub-carrier
bps = log2(M); % bits per M-ary symbol
c = 2^floor(log2(nchoosek(N,K))); % Effective Carrier Combinations
p1 = floor(log2(nchoosek(N,K)));  % index bit length per cluster
p2 = K*bps; % information bit length per cluster
p=p1+p2; % total number of bits

%% Reference M-ary and index symbols used for detection
if(K==2&&N==4)
    index_all = [1 0;2 0;3 1;3 2]; % optimal combination for this case
    %index_all = Combin_Md(N,K);
else
    index_all = Combin_Md(N,K);
end
index_allz=index_all+1;

sym_test=zeros(M,1);
for qq=1:M
    if(Mary==1)
        sym_test(qq)=pskmod(qq-1,M,ro*pi./M,'gray');
    else
        sym_test(qq)=qammod(qq-1,M,0,'gray');
    end
end

ref_sym = sym_test;
if(Mary==1)
    ref_symmm = ref_sym.*(1./abs(ref_sym)); % PSK
else
    ref_symmm = ref_sym.*(1/sqrt(QAM)); % QAM
end

%% ==================== Loop for SNR =========================
PEP = zeros(1,size(sigma,2)); % index symbol error IEP
OFDM_SER = zeros(1,size(sigma,2)); % M-ary symbol error
Total_SER = zeros(1,size(sigma,2)); % SEP overall
BER=zeros(1,size(sigma,2));
BER1=zeros(1,size(sigma,2)); % index bit error rate
BER2=zeros(1,size(sigma,2)); % M-ary bit error rate

for s1 = 1:size(sigma,2)
    fprintf('== EbN0(dB) is %g == \n',EbN0dB(s1))
    %% ==================== Loop for iteration =======================
    symerr_mcik = zeros(1,iter);
    symerr_ofdm = zeros(1,iter);
    symerr_iter= zeros(1,iter);
    BER_iter= zeros(1,iter);
    BER_iter_1= zeros(1,iter);
    BER_iter_2= zeros(1,iter);
    for s2 = 1:iter
        fprintf('== EbN0(dB) is %g and iteration is %g == \n',EbN0dB(s1),s2)
        %% ===================== Bit generator =========================
        % bit = (index bit + M-ary bps) * symbols in OFDM frame
        bit = randi([0 1],1,(p1+p2)*nSymPerFrame);
        % bit split - reshape bit stream (p1+p2)
        bit2 = reshape(bit.',p1+p2,nSymPerFrame).';
        %% ================= Index selector =========================
        % information bits (p2)
        info_bit = bit2(:,p1+1:end);      
        % M-ary data symbol
        sym=[];
        x=1;
        for i=1:K
            y=bps*i;
            info_bit_i= info_bit(:,x:y);
            x=y+1;
            info_dec_i = bi2de(info_bit_i);
            % sym_i = sym_test(info_dec_i+1);
            if(Mary==1)
                sym_i = pskmod(info_dec_i,M,ro*pi./M,'gray');
            else
                sym_i = qammod(info_dec_i,M,0,'gray');
            end
            sym(:,i)=sym_i;
        end
        % index bits (p1)
        index_bit = bit2(:,1:p1);
        % index symbol ( bit to decimal ), select indices from combinatorial method
        index_sym = BitoDe(index_bit);
        % Set average symbol power to 1
        sym_norm = sym.*(1./abs(sym));
        % Power reallocation
        sym_tx = sym_norm.*sqrt(PwrSC);
        % transmitted OFDM symbols
        tx_sym = zeros(N,nSymPerFrame);
        for kk = 1:nSymPerFrame
            kk_index = index_sym(kk)+1;
            indices = index_all(kk_index,:)+1;
            tx_sym(indices,kk) = sym_tx(kk,:);
        end
        
        % CSI error variance
        if(CSI==1)
            eps=0; % perfect CSI
        elseif(CSI==2)
            eps=var; % fixed CSI
        else
            eps=1./(1+mmse*EsN0(s1)); % variable CSI
        end
        
        %% Transmission over Rayleigh fading channel and AWGN noise
        noise = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)));
        h = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)))*sqrt(1-eps);
        e=sqrt(eps)./sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)));
        h1=h+e;
        y = sqrt(EsN0(s1))*h1.*tx_sym+noise;
        avSNR=sqrt(EsN0(s1));
        
        %% ================== ML / LLR / Greedy detectors ====================
        index_sym_de = zeros(1,nSymPerFrame);
        indices_de = zeros(nSymPerFrame,K);
        re_sym = zeros(nSymPerFrame,K);
        OutputData = zeros(N,nSymPerFrame);
        for jj=1:nSymPerFrame
            %% ML detector
            if Detect_method == 1  % ML or low complexity ML detectors (all have same performance)
                [BB,MM] = ML_Detector_LowC(avSNR,M,K,p1,PwrSC,index_all,y,h,jj,ref_sym,ref_symmm,N,Mary,QAM);
                %                 [BB,MM] = ML_Detector_NearML(avSNR,M,K,p1,PwrSC,index_allz,y,h,jj,ref_sym,ref_symmm,N);
                index_sym_de(jj) = BB-1;
                re_sym(jj,:) = MM;
                %% LLR detector
            elseif Detect_method == 2
                s=sqrt(PwrSC)*ref_symmm;
                lamda = zeros(N,1);
                for i = 1:N
                    if (LLR==1)
                        lamda(i) = LLR_Detector(y(i,jj),h(i,jj),avSNR,s);
                    else
                        lamda(i) = LLR_Detector2(y(i,jj),h(i,jj),avSNR,s);
                    end
                end
                
                [Adist,Bindex] = sort(lamda);
                AcIndex = sort(Bindex(N-K+1:N),'descend');
                AcIndex = AcIndex';
                indices_de(jj,:)=AcIndex;
                index_sym_de(jj)=-1;
                for ii=1:2.^p1
                    if(sum(index_allz(ii,:)==AcIndex)==K)
                        index_sym_de(jj) = ii-1;
                    end
                end
                [NN] = ML_Detector_Mary(avSNR,K,M,PwrSC,y,h,jj,ref_sym,ref_symmm,indices_de);
                re_sym(jj,:) = NN;
                %% Greedy Detector
            else
                Y = abs(y(:,jj));
                [Adist,Bindex] = sort(Y);
                AcIndex = sort(Bindex(N-K+1:N),'descend');
                AcIndex = AcIndex';
                indices_de(jj,:)=AcIndex;
                index_sym_de(jj)=-1;
                for ii=1:2.^p1
                    if(sum(index_allz(ii,:)==AcIndex)==K)
                        index_sym_de(jj) = ii-1;
                    end
                end
                [NN] = ML_Detector_Mary(avSNR,K,M,PwrSC,y,h,jj,ref_sym,ref_symmm,indices_de);
                re_sym(jj,:) = NN;
            end
        end
        
        %% =================error rate computation====================
        % ofdm symbol error
        ofdm_symerr = sum(sum(sym~=re_sym));
        % index symbol error
        ind_symerr = sum(index_sym~=index_sym_de);
        % index symbol to bit, index bit error
        index_bit_de = DetoBit(index_sym_de,p1);
        index_bit_err=sum(sum(index_bit~=index_bit_de));
        
        % QAM symbol to bit
        if(Mary==1)
            info_de_re=pskdemod(re_sym,M,ro*pi./M,'gray');
        else
            info_de_re=qamdemod(re_sym,M,0,'gray');
        end
        info_bit_re= zeros(nSymPerFrame,K*bps);
        for kk=1:K
            info_bit_re(:,(kk-1)*bps+1:kk*bps)=de2bi(info_de_re(:,kk),bps);
        end
        info_bit_err=sum(sum(info_bit~=info_bit_re));
        
        %% ===========symbol & bit error rate  1 iteration==========
        % MCIK sym error/ index
        symerr_mcik(s2) = ind_symerr/nSymPerFrame;
        % OFDM sym error / M-ary symbols
        symerr_ofdm(s2) = ofdm_symerr/(K*nSymPerFrame);
        % symbol error rate
        symerr_iter(s2) = (ind_symerr+ofdm_symerr)/(nSymPerFrame+K*nSymPerFrame);
        
        %%% Bit error rate BER
        BER_iter(s2)=(info_bit_err+index_bit_err)./((p1+p2)*nSymPerFrame);
        BER_iter_1(s2) = index_bit_err./p1./nSymPerFrame;
        BER_iter_2(s2) = info_bit_err./p2./nSymPerFrame;
    end
    
    %% =============average bit error rate================
    PEP(s1) = sum(symerr_mcik)/iter;
    OFDM_SER(s1) = sum(symerr_ofdm)/iter;
    Total_SER(s1) = sum(symerr_iter)/iter;
    BER(s1)= sum(BER_iter)./iter;
    BER1(s1)= sum(BER_iter_1)./iter;
    BER2(s1)= sum(BER_iter_2)./iter;
end
fprintf('¬¬ N = %g / K = %g / M = %g / ¬¬ \n',N,K,M)

semilogy(EbN0dB,BER,'r >-','LineWidth',2.5,'MarkerSize',12)
hold on
semilogy(EbN0dB,Total_SER,'b +-','LineWidth',2.5,'MarkerSize',12)
hold on

axis([0 40 10^-5 10^0])
grid on
hold on
title('')
xlabel('Es/No (dB)')
ylabel('BER/SEP')
toc