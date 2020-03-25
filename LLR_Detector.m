% it can be if n=k/2
% using Jacobian logarithm
% from OFDM-IM paper


function lamda = LLR_Detector(y,h,avSNR,ref_sym)
    M = length(ref_sym);
    s = ref_sym;    
    a = -abs(y-h.*s.*avSNR).^2;           
    output = max(a(1),a(2)) + log(1+exp(-abs(a(1)-a(2))));

    if M>3

        for i = 3:M

            input = output;

            ai = a(i);

            output = max(input,ai);


        end
    end

lamda = output+abs(y).^2;