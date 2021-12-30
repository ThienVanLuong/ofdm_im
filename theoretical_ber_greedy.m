%% Main objective of this file
% To calcuate theoretical bound for BER of Greedy Detection of MCIK-OFDM/OFDM-IM in ref. [1]

%% OFDM-IM simulation with ML, GD and LLR detectors and imperfect CSI (see file OFDM_IM.m) 
% Performance metics: SEP symbol error probability [1], BER in [2]
% Matlab version 2015b, also working well on 2019a.
% Note that MCIK-OFDM is another name of OFDM-IM see [1], [2].

%% Author information
% Thien Van Luong, Queen's University Belfast, UK, now with
% University of Southampton, UK.
% Email: thienctnb@gmail.com.
% Personal page: https://tvluong.wordpress.com
% Note that I am now with Phenikaa University, Vietnam. 

%% References
% [1] T. V. Luong and Y. Ko, “A tight bound on BER of MCIK-OFDM with
% greedy detection and imperfect CSI,” IEEE Commun. Lett., vol. 21,
% no. 12, pp. 2594 – 2597, Dec. 2017.
% [2] T. V. Luong and Y. Ko, “Impact of CSI uncertainty on MCIK-OFDM:
% tight, closed-form symbol error probability analysis,” IEEE Trans. Veh.
% Technol., vol. 67, no. 2, pp. 1272 – 1279, Feb. 2018.


%% ==============================OFDM-IM/MCIK-OFDM=================================
clear

M=4;
N=4;
K=1;

var = 0.01;
mmse = 1;
CSI=1;

Detect_method = 3;
LLR = 1;


Plot_type = 1;
ro=0;

Mary=1; % 1 PSK, 2 QAM
if(M==8)
    QAM = (5*M-4)./6;
else
    QAM = (2/3)*(M-1);
end

if(M==2)
    xi=1;
else
    xi=2;
end

L=1;

%% ======================= Misc Parameters ================================
iter = 2;  % Iterations
nSymPerFrame = 1e4; % Number of symbol per frame(1 OFDM symbol)
EbN0dB = 0:5:40;
EbN0 = 10.^(EbN0dB/10);

% Es/N0 parameter
PwrSC = N/K; % Average Tx power per active sub-carrier
QAM_Scaling_Factor = (2/3)*(M-1); 
bps = log2(M); % bits per symbol 
EsN0dB = EbN0dB; % + 10*log10(bps);%+10*log10(1/PwrSC);
EsN0 = 10.^(EsN0dB/10);
c = 2^floor(log2(nchoosek(N,K))); % Effective Carrier Combinations
p1 = floor(log2(nchoosek(N,K)));  % index bit length per cluster
p2 = K*bps; % information bit length per cluster
p2_re=bps;
p=p1+p2;
sigma = sqrt(1./EsN0);

%% ==================== Loop for SNR =========================
PEP = zeros(1,size(sigma,2)); % index symbol error
OFDM_SER = zeros(1,size(sigma,2)); % ofdm symbol error
Total_SER = zeros(1,size(sigma,2));
BER=zeros(1,size(sigma,2));
BER1=zeros(1,size(sigma,2));
BER2=zeros(1,size(sigma,2));

A=(1:M)';
for x=1:K-1
    E=cell(1,M);
    for i=1:M
        C=ones(length(A),1)*i;
        E{i}=[C A];
    end
    combs=cell2mat(E');
    A=combs;
end
%       sym_test=zeros(1,M);
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
        
        % data symbol              
        info_bit_=[];
        info_dec_=[];
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
        
        if(K==2&&N==4)
            index_all = [1 0;2 0;3 1;3 2];
            %index_all = Combin_Md(N,K);
        else
            index_all = Combin_Md(N,K);
        end
        % Set average symbol power to 1
        if(Mary==1)
            sym_norm = sym.*(1./abs(sym));
        else
            sym_norm = sym./sqrt(QAM);
        end
       
        % Power reallocation
        sym_tx = sym_norm.*sqrt(PwrSC);
        % transmitted OFDM symbols
        tx_sym = zeros(N,nSymPerFrame);
        for kk = 1:nSymPerFrame            
            kk_index = index_sym(kk)+1;            
            indices = index_all(kk_index,:)+1;
            tx_sym(indices,kk) = sym_tx(kk,:);
        end        
        index_allz=index_all+1;
        MCIK_LUT=zeros(c,N);        
        for row=1:c
            for col=1:K
                MCIK_LUT(row,index_allz(row,col))=1;
            end
        end
               
        if(CSI==1)
            eps=0;
        elseif(CSI==2)
            eps=var;
        else
            eps=1./(1+mmse*EsN0(s1));
        end
        avSNR=sqrt(EsN0(s1));
        noise = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)));        
        h = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)))*sqrt(1-eps);                     
        e=sqrt(eps)./sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym))); 
        h1=h+e; 
        y = avSNR*h1.*tx_sym+noise;         
         
        
        %% ================== ML / LLR / Greedy detect ====================
        index_sym_de = zeros(1,nSymPerFrame);
        indices_de = zeros(nSymPerFrame,K);
        re_sym = zeros(nSymPerFrame,K);
        OutputData = zeros(N,nSymPerFrame);   
        tic
        for jj=1:nSymPerFrame
            %% ML detector
            if Detect_method == 1               
                %[BB,MM] = ML_Detector_MCIK_PSK(avSNR,M,K,p1,PwrSC,index_all,y,h,N,jj,A,sym_test);
                [BB,MM] = ML_Detector_MCIK_LowC(avSNR,M,K,p1,PwrSC,index_all,y,h,jj,ref_sym,ref_symmm,N,Mary,QAM);
                %[BB,MM] = ML_Detector_MCIK_NearML(avSNR,M,K,p1,PwrSC,index_allz,y,h,jj,ref_sym,ref_symmm,N,Mary,QAM);
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
%                 [MLsymLUT,minInd] = ML_Detector_PSK(avSNR,M,K,PwrSC,y,h,jj,Diversity_Technique,L,indices_de,ro);
%                 re_sym(jj,:) = MLsymLUT(minInd);   
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
        % MCIK sym error
        symerr_mcik(s2) = ind_symerr/nSymPerFrame;
        % OFDM sym error
        symerr_ofdm(s2) = ofdm_symerr/(K*nSymPerFrame);
        % symbol error rate
        symerr_iter(s2) = (ind_symerr+ofdm_symerr)/(nSymPerFrame+K*nSymPerFrame);
        
        %%% Bit error rate BER
        BER_iter(s2)=(info_bit_err+index_bit_err)./((p1+p2)*nSymPerFrame);
        BER_iter_1(s2) = index_bit_err./p1./nSymPerFrame;
        BER_iter_2(s2) = info_bit_err./p2./nSymPerFrame;
%         toc
%         toc/nSymPerFrame
    end    
    
    %% =============average bit error rate================
    PEP(s1) = sum(symerr_mcik)/iter;
    OFDM_SER(s1) = sum(symerr_ofdm)/iter;
    Total_SER(s1) = sum(symerr_iter)/iter; 
    BER(s1)= sum(BER_iter)./iter;
    BER1(s1)= sum(BER_iter_1)./iter;
    BER2(s1)= sum(BER_iter_2)./iter;
end

fprintf('¬¬ N = %g / K = %g / M = %g / L = %g ¬¬ \n',N,K,M,L)

%% display plot figure
if Plot_type == 1
    plot = PEP;
else
    plot = Total_SER;
end

%% Theoretical bound caculation of MCIK-OFDM/OFDM-IM, for bit error rate (BER)
SNR_av=EsN0*PwrSC;
PEP_theo=zeros(size(EbN0dB));
SEP_theo=zeros(size(EbN0dB));
SEP_theo_1=zeros(size(EbN0dB));

PEP_GD=zeros(size(EbN0dB));
PEP_GD_ex=zeros(size(EbN0dB));
SEP_GD=zeros(size(EbN0dB));

PM=zeros(size(EbN0dB));

SEP_theo_asym=zeros(size(EbN0dB));
SEP_GD_asym=zeros(size(EbN0dB));

BER_theo=zeros(size(EbN0dB));
BER_theo_1=zeros(size(EbN0dB));
BER_theo_2=zeros(size(EbN0dB));

BER_asym=zeros(size(EbN0dB));


MM=(sin(pi./M)).^2;
B1=13*xi./48./MM;
B2=B1*2;
w1=(2-1./M)*K*(N-K)./(K+1);
w2 = K./(K+1);
v1=(2-1./M);

%v2=1+bps*(1-1./M);
if(p1==1)
    v2=(1+bps./2);
else
    v2=(p1./2+bps./2);
end

v2=(1+bps./2);

%a = 1-eps;
GDindexgain = 0;
for i=1:N-K
    GDindexgain = GDindexgain + (-1).^(i+1).*nchoosek(N-K,i)./i;
end
omega = GDindexgain;
GDindexgain=v1*GDindexgain;
if(p1==1)
    nu = 2;
else
    nu = 1;
end
pm=nu*p1+bps;
pc=K./(2*p*N);
pd=K*xi./12./p;

for i=1:length(EbN0dB)
     snr=SNR_av(i);
     Es = EsN0(i);
        if(CSI==1)
            eps=0;
        elseif(CSI==2)
            eps=var;
        else
            eps=1./(1+mmse*EsN0(i));
        end
        a=1-eps;
    PM(i)=xi*(1./(1+(1-eps)*MM*snr./(1+eps*snr./1))./12+1./(1+4./3.*MM*(1-eps)*snr./(1+eps*snr./1))./4);
    
    if(Detect_method==1) % ML detection
      % PEP_theo(i)=(K*(N-K)./2)*(1./(1+(1-eps)*snr./(4*(1+eps*snr./2))).^2);
       PEP_theo(i)=(K*(N-K)./12)*(1./(1+(1-eps)*snr./(4*(1+eps*snr./2))).^2+3./(1+(1-eps)*snr./(3*(1+eps*snr./2))).^2);  
       SEP_theo(i)=(v1*PEP_theo(i)+K*PM(i))./(K+1);
       if(CSI==1)
           SEP_theo_asym(i)=B1*w2./snr;
       elseif(CSI==2)          
           SEP_theo_asym(i)=w1./12.*((1+a./2./eps).^(-2)+3*(1+2*a./3./eps).^(-2))+xi*w2./12.*(1./(1+a*MM./eps)+3./(1+4*a*MM./3./eps));
       else
           SEP_theo_asym(i)=B1*(N+K)./(K+1)./snr;
       end  
       BER_theo(i) = (v2*PEP_theo(i)+K*PM(i))./(p1+p2);
       BER_theo_1(i) = PEP_theo(i)./2;
       BER_theo_2(i) = PEP_theo(i)./K./2+PM(i)./bps;
    else   %Greedy detection        
        PEP_GD(i)=K*(N-K)./(2+SNR_av(i));
        pe=zeros(1,N-K);
        for j=1:N-K
            pe(j)= K*nchoosek(N-K,j)*(-1).^(j+1)./((j+1)*(1+j*SNR_av(i)./(j+1)));          
        end
        PEP_GD_ex(i)=sum(pe);
        join = xi*K*(N-K)./24.*(1./(1+(0.5+MM)*snr)+3./(1+(0.5+4*MM./3)*snr));       
        SEP_GD(i)=v1*PEP_GD_ex(i)./(K+1)+K.*PM(i)./(K+1);%-join./(K+1);  %-join./(K+1)
       if(CSI==1)
           SEP_GD_asym(i)=(GDindexgain+B1)*w2./snr;
       elseif(CSI==2)          
           SEP_GD_asym(i)=xi*w2./12.*(1./(1+a*MM./eps)+3./(1+4*a*MM./3./eps));
       else
           SEP_GD_asym(i)=(B1*(N+K)+K*GDindexgain)./(K+1)./snr;
       end
       BER_theo(i) = (v2*PEP_GD_ex(i)+K*PM(i))./(p1+p2);
       BER_theo_1(i) = p1*PEP_GD_ex(i)./2;
       BER_theo_2(i) = PEP_GD_ex(i)./K./2+PM(i)./bps;
       if(CSI==1)
           BER_asym(i)=K*pc*(omega*pm+B2)./Es;
       elseif(CSI==2)
           BER_asym(i)=pd*(1./(1+a*MM./eps)+3./(1+4*a*MM./3./eps));
       else
           BER_asym(i)=pc*(K*omega*pm+B2*(N+K))./Es;
       end
       
    end    
    
end


figure (23)
semilogy(EbN0dB,BER,'r >-','LineWidth',1.5,'MarkerSize',10)
hold on
semilogy(EbN0dB,BER_theo,'b >:','LineWidth',1.5,'MarkerSize',10)
hold on
semilogy(EbN0dB,BER_asym,'b --','LineWidth',1.5,'MarkerSize',10)
hold on
legend('Simulation, GD','Theo, GD','Asym, GD')

axis([0 40 10^-5 10^0])
grid on
hold on
title('')
xlabel('Es/No (dB)')
if Plot_type == 1
    ylabel('Average IEP')
elseif(Plot_type == 2)
    ylabel('Average SEP')
else
    ylabel('BER')
end
toc
