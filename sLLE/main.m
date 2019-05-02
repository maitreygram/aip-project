close all;
clear;

% Read image
P = phantom('Modified Shepp-Logan',200);
imshow(P)

% Get projections
numAngles = 360;
angles = rand(numAngles,1) * 360; % random angles
% angles = randperm(numAngles)*360/numAngles;
% angles = (1:numAngles)*360/numAngles;
% past = load('sLLE_360_rand.mat');
% angles = past.angles;

% angles = randperm(360);
R = radon(P, angles);

% estimate angles

order = sLLE_ang_est(R');

uni_ang = (1:numAngles)*360/numAngles;
recons_P = iradon(R(:,order),uni_ang,'Ram-Lak');
recons_P = recons_P(2:end-1,2:end-1);
% show result
figure, imshow(mat2gray(recons_P));
l2_norm = norm(recons_P-double(P))