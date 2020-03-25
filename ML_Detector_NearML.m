% ============Near ML Maximum Likelihood Detector
function [BB,MM] = ML_Detector_NearML(avSNR,M,K,p1,PwrSC,index_allz,y,h,jj,ref_sym,ref_symmm,N)
 
  sym_k = zeros(1,N);
  dis = zeros(N,1);
  for k1=1:N        
        disA = zeros(1,M);
        r=y(k1,jj)./avSNR./h(k1,jj)./sqrt(PwrSC);
        for k2=1:M
            sym_m = ref_symmm(k2);
            disA(k2)=(abs(r-sym_m)).^2;
        end       
        [~,iA] = min(disA);
        sym_k(k1)=ref_sym(iA);
        s=ref_symmm(iA);
        
        hh=avSNR*sqrt(PwrSC)*h(k1,jj);
        dis(k1)=(abs(s*hh)).^2-2*real(conj(y(k1,jj))*hh*s);
       %  dis(k1)=abs(y(k1,jj)-s*hh).^2;
  end  

    [~,Bindex] = sort(dis);
    AcIndex = sort(Bindex(1:K),'descend'); 
    AcIndex = AcIndex';
    BB=-1;
    for ii=1:2^p1
        if(sum(index_allz(ii,:)==AcIndex)==K)
            BB = ii; 
        end
    end  
                
MM = sym_k(AcIndex);
