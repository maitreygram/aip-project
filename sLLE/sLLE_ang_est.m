function [ord_idx] = sLLE_ang_est(projections,epsilon)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[N,P] = size(projections);


%% Fourier Transform of Each Projection
FTprojections = zeros(N,P);
for i = 1:N
    FTprojections(i,:) = fft(projections(i,:));
end
%% Construction of Weight Matrix

W = zeros(N,N);

for i = 1:N
    for j = i:N
        W(i,j) = exp(-norm(FTprojections(i,:) - FTprojections(j,:))/(2*epsilon));
        W(j,i) = W(i,j);
    end
end

%% Dimensionality reduction goes here

temp_mat = eye(N)- W;
M = temp_mat' * temp_mat;
Y = zeros(N,2);
B = eye(N);
while True
%     Y STEP: min tr(Y'MY) s.t. Y'BY = I
    [V,D] = eig(B\M);
    [~,idx] = sort(diag(D),'descend');
    Y = V(:,idx(1:2));   

%     B STEP: sqrt(B)YY'sqrt(B) all diag elem = 1
    B_temp = B;
    for i = 1:N
        B(i,i) = 1/(Y(i,:)*Y(i,:)');
    end
    
    if norm(B_temp -B) < break_err
        break;
    end
    
end

% Final Embeddings
Z = sqrt(B)*Y;

ini_ang = zeros(N,1);
for i = 1:N
    ini_ang(i) = atan(Z(i,2)/Z(i,1));
end

[~,ord_idx] = sort(ini_ang);

end

