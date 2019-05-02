close all;
clear;

% Read image
P = phantom('Modified Shepp-Logan',200);
P = mat2gray(P);

% Get projections
numAngles = 2048;
% angles = rand(numAngles,1) * 360; % random angles
angles = (1:numAngles)'*360/numAngles;

Perm = randperm(length(angles))';
angles = angles(Perm);

R = radon(P, angles);
[Z,~] = size(R);

% weight matrix
W = zeros(numAngles);
D = zeros(numAngles);
for i=1:numAngles
    for j=i:numAngles
        W(i,j) = norm(R(:,i) - R(:,j))^2;
        W(j,i) = W(i,j);
    end
end

epsilon = 50;
Wt = exp(-W/(2*epsilon));

for i=1:numAngles
    D(i,i) = sum(Wt(i,:));
end

% Wt = normc(normr(W));
L = D - Wt;

% eigen vectors
[EigVec,EigVal] = eig(L,D);
EigVal = real(EigVal);
EigVec = real(EigVec);

absEigVal = max(abs(EigVal));

[SortEig,SortInd] = sort(absEigVal);

trivialVec = ones(numAngles,1);
trivialVec = trivialVec/norm(trivialVec);

i = 1;
while norm(abs(EigVec(:,SortInd(i))/norm(EigVec(:,SortInd(i)))) - trivialVec) < 1e-10
    i = i + 1;
end

Phi1 = EigVec(:,SortInd(i));
Phi2 = EigVec(:,SortInd(i+1));

% redDimn = Phi2./Phi1;
Psi = atan2d(Phi2,Phi1);
[~,PsiOrd] = sort(Psi);

anglesEst = (1:numAngles)*360/numAngles;

recons_P = iradon(R(:,PsiOrd),anglesEst);
recons_P = recons_P(2:end-1,2:end-1);

% show result
figure, imshow(mat2gray(recons_P));

% % %%
% % L1 = D - W;
% % [EigVec1,EigVal1] = eig(L1,D);
% % 
% % Y = EigVec1(:,1:2);
% % anglesTan = Y(:,2)./Y(:,1); % check
% % 
% % anglesEst = atand(anglesTan) + 90.0;
% 
% [~,angOrd] = sort(angles);
% [~,PsiOrd] = sort(Psi);
% % [~,estOrd] = sort(anglesEst);
% 
% anglesEst = (1:numAngles)*180/numAngles;
% REst = R(:,PsiOrd);
% 
% recons_P = iradon(REst,anglesEst);
% recons_P = recons_P(2:end-1,2:end-1);
% 
% % show result
% figure, imshow(mat2gray(recons_P));
% l2_norm = norm(recons_P-double(P))