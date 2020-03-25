

function dec = BitoDe(bin)

row = size(bin,1);

for j = 1:row
bin1 = bin(j,:);
len = length(bin1);
i = 1:len;
b = 2.^(len-i);

dec(j) = sum(b.*bin1(i));

end