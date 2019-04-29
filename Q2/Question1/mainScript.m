im = double(imread('barbara256.png'));
[row,col] = size(im); 
%% Part(a) Part (b) and Part (c)
im_noisy = double(im) + 2*randn(row,col); % For part (a)

% Phi for part (b) and (c)
p_row = 32;
p_col = 64;
phi = randn(p_row,p_col);
% 2-d DCT matrix for part (a),(b)
U = kron(dctmtx(8),dctmtx(8));
% Haar Wavelet basis for part(c)
% H = 3;
% lamdas for (a),(b),(c)
reg_a = 1;
reg_b = 1;
reg_c = 1;

recon_with_overlap_a = zeros(row,col);
recon_with_overlap_b = zeros(row,col);
recon_with_overlap_c = zeros(row,col);

result_img_a = zeros(row,col);
result_img_b = zeros(row,col);
result_img_c = zeros(row,col);
%%
for i = 1:row-7
    for j = 1:col-7
        yp_a = reshape(im_noisy(i:i+7,j:j+7),64,1);
        xp_bc = reshape(im(i:i+7,j:j+7),64,1);
        yp_bc = phi*xp_bc;
        
        % for part (a)
        A_a = U;
        alpha_a = max(eig(A_a'*A_a));
        theta_a = ISTA(yp_a,A_a,reg_a,0.1,alpha_a,0); % second last argument is extent of convergence parameter
        xexp_a = U*theta_a;
        recon_with_overlap_a(i:i+7,j:j+7) = recon_with_overlap_a(i:i+7,j:j+7) + reshape(xexp_a,8,8);
%         
        % with 2d - DCT basis
        A_b = phi*U;
        alpha_b = max(eig(A_b'*A_b));
        theta_b = ISTA(yp_bc,A_b,reg_b,0.1,alpha_b,0);
        xexp_b = U*theta_b;
        recon_with_overlap_b(i:i+7,j:j+7) = recon_with_overlap_b(i:i+7,j:j+7) + reshape(xexp_b,8,8); 
        
%         % with Haar wavelet basis
        A_c = phi;
        alpha_c = max(eig(A_c'*A_c));
        theta_c = ISTA(yp_bc,A_c,reg_c,10,alpha_c,1);
        xexp_c = idwt2on1d(theta_c);
        recon_with_overlap_c(i:i+7,j:j+7) = recon_with_overlap_c(i:i+7,j:j+7) + reshape(xexp_c,8,8); 
%         fprintf("Patch No.(%d,%d)\n",i,j)
    end
    fprintf("Patch No.(%d,--)\n",i)
end

%%
for i = 1:row
    for j = 1:col
        result_img_a(i,j) = recon_with_overlap_a(i,j)/(min(min(i,row-i+1),8)*min(min(j,1+col-j),8));
        result_img_b(i,j) = recon_with_overlap_b(i,j)/(min(min(i,row-i+1),8)*min(min(j,1+col-j),8));
        result_img_c(i,j) = recon_with_overlap_c(i,j)/(min(min(i,row-i+1),8)*min(min(j,1+col-j),8));
    end
end

%%
figure;
imshow(mat2gray(result_img_a))
title("Part (a)")
figure;
imshow(mat2gray(result_img_b))
title("Part (b)")
figure;
imshow(mat2gray(result_img_c))
title("Part (c)")

%%
RMSE_a = sqrt(immse(result_img_a,double(im)))/sqrt(immse(double(im),zeros(row,col)));
RMSE_b = sqrt(immse(result_img_b,double(im)))/sqrt(immse(double(im),zeros(row,col)));
RMSE_c = sqrt(immse(result_img_c,double(im)))/sqrt(immse(double(im),zeros(row,col)));

fprintf('Root Mean squared error of reconstruction (a)= %d\n',RMSE_a);
fprintf('Root Mean squared error of reconstruction (b)= %d\n',RMSE_b);
fprintf('Root Mean squared error of reconstruction (c)= %d\n',RMSE_c);
