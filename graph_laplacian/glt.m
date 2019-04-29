close all;
clear;

% Read image
P = phantom('Modified Shepp-Logan',200);
% imshow(P)

% Get projections
numAngles = 180;
% angles = rand(numAngles,1) * 180; % random angles
angles = 1:180;

R = radon(P, angles);
[Z,~] = size(R);

W = zeros(numAngles);
D = zeros(numAngles);
for i=1:numAngles
    for j=1:numAngles
        W(i,j) = norm(R(:,i) - R(:,j))^2;
    end
    D(i,i) = sum(W(i,:));
end

epsilon = 0.05;
W = exp(-W/epsilon);

