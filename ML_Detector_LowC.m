% ============LowC ML Maximum Likelihood Detector
function [BB,MM] = ML_Detector_LowC(avSNR,M,K,p1,PwrSC,index_all,y,h,jj,ref_sym,ref_symmm,N,Mary,QAM)

dis_id = zeros(1,2^p1);
dis_sym = zeros(K,2^p1);

for bb=1:2^p1
  id = index_all(bb,:)+1; 
  sym_k = zeros(K,1);
  sum = 0;
  for k1=1:K        
        n1 = id(k1);
        disA = zeros(1,M);
        for k2=1:M
            sym_m = sqrt(PwrSC).*ref_symmm(k2);
            disA(k2)=norm(y(n1,jj)-avSNR*h(n1,jj)*sym_m).^2;
        end       
        [minA,iA] = min(disA);
        sum = sum + minA;   % abs(y-sym_hat*H)     
        sym_k(k1)=ref_sym(iA);     
  end  
   % dis_id(bb)=sum;
    dis_sym(:,bb)=sym_k; %%
%%  New Low-Com ML detector  
      sym_b = zeros(N,1);       
      sym_b(id) = sym_k; %
      if(Mary==1)
            tmp = norm(y(:,jj)-avSNR.*sqrt(PwrSC).*h(:,jj).*sym_b).^2;     
      else
            tmp = norm(y(:,jj)-avSNR.*sqrt(PwrSC).*h(:,jj).*sym_b./sqrt(QAM)).^2;  
      end
      dis_id(bb)=tmp;

end

[~, I] = min(dis_id);
BB = I; 
MM = dis_sym(:,I).';