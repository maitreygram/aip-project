function [ord_idx] = sLLE_ang_est(projections)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[N,P] = size(projections);


%% Fourier Transform of Each Projection
FTprojections = zeros(N,P);
for i = 1:N
    FTprojections(i,:) = abs(fft(projections(i,:)));
end
%% Construction of Weight Matrix

% weights = load('sLLE_360_rand.mat');
% W = weights.W; 

W = zeros(N,N);

for i = 1:N
    k = 60;
    tempPro = FTprojections;
    tempPro(i,:) = [];
    Idx = knnsearch(tempPro,FTprojections(i,:),'K',k);
%     Idx(1) = [];
    cvx_begin quiet
        variable w(k)
        minimize(norm(FTprojections(i,:) - w'*tempPro(Idx,:)))
        subject to 
            sum(w) == 1
    cvx_end
    
    q = find(Idx >= i);
    Idx(q) = Idx(q) + 1;
    W(i,Idx) = w;
end

% epsilon = 10^6;
% 
% for i = 1:N
%     for j = i+1:N
%         W(i,j) = exp(-norm(FTprojections(i,:) - FTprojections(j,:))^2/(2*epsilon));
%         W(j,i) = W(i,j);
%     end
% end

% % implement equation 8b given in paper  
% W  = W./sum(W,2);
%% Dimensionality reduction goes here

temp_mat = eye(N)- W;
M = temp_mat' * temp_mat;
Y = zeros(N,2);
B = eye(N);

% break_err  = 0.1;
itr = 0;
while true
%     Y STEP: min tr(Y'MY) s.t. Y'BY = I
    Y_temp = Y;
    [V,D] = eig(B\M);
    [~,idx] = sort(diag(real(D)));
    Y = real(V(:,idx(1:2)));
%     Y(:,1) = V(:,1);
%     Y(:,2) = V(:,2);

%     B STEP: sqrt(B)YY'sqrt(B) all diag elem = 1
    B_temp = B;
    for i = 1:N
        
        B(i,i) = 1/(Y(i,:)*Y(i,:)');
    end
%     fprintf('Error in B = %d\n',norm(B_temp - B));
    fprintf('Error in Y = %d   Error in B = %d\n',norm(Y_temp - Y,'fro'),norm(B_temp - B,'fro'));
    itr = itr + 1;
    if itr > 100
        break;
    end
end

% Final Embeddings
Z = sqrt(B)*Y;

ini_ang = zeros(N,1);
for i = 1:N
    ini_ang(i) = atan2d(Z(i,2),Z(i,1));
end

[~,ord_idx] = sort(ini_ang);

end

