x = zeros(100,1);
x([5,36,76,12,78,45,89,100,22,67]) = 10*randn(10,1);
[size_x,~] = size(x);
h = [1,2,3,4,3,2,1]/16;
[~,size_h] = size(h);
y_true = conv(x,h);
[size_y,~] = size(y_true);
y = y_true + (0.05 * norm(x))*randn(size_y,1); % this is producing RMSE in y and y_true around 1.0833. Hence resulting in High rmse value of reconstruction
% y = y_true;
% y([5,36,76,12,78,45,89,100,22,67]) = y([5,36,76,12,78,45,89,100,22,67]) + (0.05 * norm(x))*randn(10,1);
A_temp = zeros(size_y,size_y+size_h-1);
for i = 1:size_y
    A_temp(i,i:i+size_h-1) = h;
end

A = A_temp(:,size_h:size_h + size_x -1);

reg = 1;
alpha = max(eig(A'*A));
x_recon = ISTA(y,A,reg,0.001,alpha,0);

RMSE = sqrt(sum((x_recon - x).^2)/sum(x.^2));
fprintf("Root Mean squared error in mesurement i.e. y = %d\n",sqrt(sum((y-y_true).^2)/sum(y_true.^2)));
fprintf('Root Mean squared error of reconstruction = %d\n',RMSE);