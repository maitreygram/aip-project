close all;
image = imread("slice_50.png"); % read image
[X,Y] = size(image);
Z = sqrt( X^2 + Y^2 );

rho_ind = -Z/2:Z/2;
% angles = [0:179] * pi/180;
angles = rand(18,1) * pi;

L = length(rho_ind);
num_angles = length(angles);

X_ind = -floor(X/2):floor(X/2);
Y_ind = -floor(Y/2):floor(Y/2);

R = false(L*num_angles,X*Y);
Rt = false(X*Y,L*num_angles);

for theta = 1 : num_angles

    if( angles(theta) >= pi/4 && angles(theta) <= 3*pi/4 )

        t =  -Z:1/sqrt(2):+Z;
        T = length(t);

        rho_mat = repmat(rho_ind(:),1,T);
        ty = repmat(t(:)',L,1);
        tx = (rho_mat - ty*cos(angles(theta)))./sin(angles(theta));

        for rho = 1:L
            X_indices = round(ty(rho,:) - min(Y_ind)) + 1;
            Y_indices = round(tx(rho,:) - min(X_ind)) + 1;  
            indices = X_indices >= 1 & X_indices <= Y & Y_indices >= 1 & Y_indices <= X;

            R((theta-1)*L + rho, unique(sub2ind([X Y], Y_indices(indices), X_indices(indices)))) = 1;
            Rt(unique(sub2ind([X Y], Y_indices(indices), X_indices(indices))),(theta-1)*L + rho) = 1;
        end

    else

        t =  -Z:1/sqrt(2):+Z;
        T = length(t);

        rho_mat = repmat(rho_ind(:)',T,1);    
        tx = repmat(t(:),1,L);
        ty = (rho_mat - tx*sin(angles(theta)))./cos(angles(theta));

        for rho = 1:L
            X_indices = round(ty(:,rho) - min(Y_ind)) + 1;
            Y_indices = round(tx(:,rho) - min(X_ind)) + 1;  
            indices = X_indices >= 1 & X_indices <= Y & Y_indices >= 1 & Y_indices <= X;

            R((theta-1)*L + rho, unique(sub2ind([X Y], Y_indices(indices), X_indices(indices)))) = 1;
            Rt(unique(sub2ind([X Y], Y_indices(indices), X_indices(indices))),(theta-1)*L + rho) = 1;
        end
    end
end

R = double(sparse(R));
Rt = double(sparse(Rt));

I = false(X*Y,X*Y);
A = false(L*num_angles,X*Y);
At = false(X*Y,L*num_angles);

for i=1:X*Y
    I(i,i) = 1;
end
I = double(sparse(I));

for i=1:X*Y
    temp = double(sparse(dct2(I(:,i))));
    A(:,i) = R*temp;
    At(i,:) = temp'*Rt;
    if mod(i,1000) == 0
        i
    end
end
% U = double(sparse(U));
% Ut = double(sparse(Ut));

% A = R*U;
% At = Ut*Rt;
% U = dctmtx(X*Y);
%%
y=R*double(image(:));
y = mat2gray(y);

[x,status] = l1_ls(A,At,L*num_angles,X*Y,y,1,0.0001);

%%

recons_img = idct2(x);
recons_img = reshape(recons_img,X,Y);
imshow(mat2gray(recons_img));