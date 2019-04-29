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
[X,Y] = size(P);
[Z,~] = size(R);
img_center = [ceil(X/2),ceil(Y/2)];
proj_center = ceil(Z/2);

ram_lak = zeros(Z,1);
for i=1:Z/2
    ram_lak(i) = Z/2 - i;
    ram_lak(Z+1-i) = ram_lak(i);
end

F_R = fft(R); % fourier transform of all radon transforms
F_R = fftshift(F_R);
F_R = F_R .* ram_lak;
F_R = ifftshift(F_R);
R_new = ifft(F_R);

recons_P = zeros(X,Y);
for i=1:X
    for j=1:Y
        for k = 1:numAngles
            pixel_loc = round((i-img_center(1))*sind(angles(k)) + (j-img_center(2))*cosd(angles(k))) + proj_center;
            recons_P(i,j) = recons_P(i,j) + R_new(pixel_loc,k);
        end
    end
end

% show result
figure, imshow(mat2gray(recons_P));
l2_norm = norm(recons_P-double(P))