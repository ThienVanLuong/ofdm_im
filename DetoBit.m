
% Input=0:M^K-1   ///   a=K

function output = DetoBit(input,a)
% =========  십진수 --> 이진 행렬 ===================

% binary = zeros(1,)
% dec decimal;
% a : bit length;
for ii = 1:length(input)
    dec = input(ii);
    binary = zeros(1,a);
    i=1;
    while dec > 0 
        e = rem(dec,2);

        if e == 0 
            binary(i)=0;
            dec = dec/2;
        else
            binary(i)=1;
            dec = floor(dec/2);
        end 
        i=i+1;
    end

    binary = binary(end:-1:1);

    if a~=length(binary)
        binary = [zeros(1,a-length(binary)) binary];
    end
    output(ii,:) = binary;
end