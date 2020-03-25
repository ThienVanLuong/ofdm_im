
% Combinatorial Method

% nn : number of subcarriers
% kk : number of active subcarriers
% output : all possible realization of active subcarrier indices


function [output] = Combin_Md(nn,kk)

n = nn;  
k = kk;

% all possible realization, binomial
N = binomial(n,k); 

% For loop for all J sequence
JJ = zeros(N,k);
for ii = 1:N

% Initialization
n = nn;
k = kk;
N = binomial(n,k);

% Bit length of indices
L = floor(log2(N));

% J sequence
Z = binomial(n,k)-ii;
tmp = zeros(1,k); 
    
         % For loop to find combination sum of J(ii)
         for jj = 1:k
                 
            ck = 0;
            C = 0; 
            
            % While loop to find maximum of ck while C <=Z
            while C <= Z
                
            C = binomial(ck,k);
             
            ck = ck+1;
             
            end
            
            tmp(jj) = ck-2;
          % For next element  
          Z = Z - binomial(tmp(jj),k);      
          % k is second element of combination
          k = k - 1;          
         end
         
    JJ(N-ii+1,:) = tmp;
    
end

output = JJ;