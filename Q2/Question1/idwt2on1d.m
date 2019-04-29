function y = idwt2on1d(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[len,~] = size(x);
new_size = sqrt(len/4);
cA = reshape(x(1:len/4),new_size,new_size);
cH = reshape(x(len/4 + 1:len/2),new_size,new_size);
cV = reshape(x(len/2 + 1:3*len/4),new_size,new_size);
cD = reshape(x(3*len/4 + 1:len),new_size,new_size);
y = reshape(idwt2(cA,cH,cV,cD,'db1'),len,1);
end

