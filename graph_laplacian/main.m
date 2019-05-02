close all;
clear;

% Read image
P = phantom('Modified Shepp-Logan',200);
imshow(P)

% Get projections
numAngles = 180;
angles = rand(numAngles,1) * 180; % random angles
% angles = 1:180;

R = radon(P, angles);

% estimate angles
anglesEst = angles;

% reconstruct using fbp
recons_P = iradon(R,anglesEst);
recons_P = recons_P(2:end-1,2:end-1);

% show result
figure, imshow(mat2gray(recons_P));
l2_norm = norm(recons_P-double(P))