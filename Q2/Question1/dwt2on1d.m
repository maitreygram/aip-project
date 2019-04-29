function y = dwt2on1d(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Reshape x to its 2d form
[len,~] = size(x);
X = reshape(x,sqrt(len),sqrt(len));
[cA,cH,cV,cD] = dwt2(X,'db1');
[row,~] = size(cA);
y = [reshape(cA,row^2,1)',reshape(cH,row^2,1)',reshape(cV,row^2,1)',reshape(cD,row^2,1)' ]';

end

