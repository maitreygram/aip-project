close all;
image = imread("slice_55.png"); % read image
[X,Y] = size(image);

angles = rand(18,1) * 180; % random angles
% angles = 1:180;
R_Tr = radon(image,angles); % get radon transforms
[Z,~] = size(R_Tr);
image_center = [ceil(X/2),ceil(Y/2)];
projection_center = ceil(Z/2);

F_R = fft(R_Tr); % fourier transform of all radon transforms

% ram-lak filter
ram_lak = zeros(Z,1);
for i=1:Z/2
    ram_lak(i) = Z/2 - i;
    ram_lak(Z+1-i) = ram_lak(i);
end

F_R = fftshift(F_R);
F_R = F_R .* ram_lak;
F_R = ifftshift(F_R);
R_Tr_new = ifft(F_R);

recons_Img = zeros(X,Y);
for i=1:X
    for j=1:Y
        for k = 1:18
            pixel_loc = round((i-image_center(1))*sind(angles(k)) + (j-image_center(2))*cosd(angles(k))) + projection_center;
            recons_Img(i,j) = recons_Img(i,j) + R_Tr_new(pixel_loc,k);
        end
    end
end

% sum_squared_error = sum(sum((mat2gray(recons_Img)-double(image)).^2))
figure, imshow(mat2gray(recons_Img));