% it can be if n=k/2
% using Jacobian logarithm
% from OFDM-IM paper


function lamda = LLR_Detector2(y,h,avSNR,ref_sym)
    s = ref_sym;    
    a = abs(y-h.*s.*avSNR).^2;  
    lamda=-min(a)+abs(y).^2;
   