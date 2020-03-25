function [NN] = ML_Detector_Mary(avSNR,K,M,PwrSC,y,h,jj,ref_sym,ref_symmm,indices_de)

sym_norm = ref_symmm;
idx = indices_de(jj,:);

sym_m = zeros(1,K);
for i=1:K
    n=idx(i);
    dis = zeros(1,M);
    for k=1:M
        dis(k)=abs(y(n,jj)-avSNR*sqrt(PwrSC).*h(n,jj).*sym_norm(k)).^2;
    end
    [~,I] = min(dis);
    sym_m(i)=ref_sym(I);
end
NN=sym_m;

