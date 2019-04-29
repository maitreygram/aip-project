P = imread("slice_50.png");
P = P(1:10,1:10);
[X,Y] = size(P);
P = padarray(P,[max(0,Y-X) max(0,X-Y)], 'post'); % make it square
imshow(P)
title('Original image')

Q = reshape(P,[],1);

R = radon(P,0:179);

A = radon(eye(19*309,X*Y),0:179);
R2 = A*Q;
R2 = reshape(R2,309,180);

I1 = iradon(R,0:179);
I2 = iradon(R,0:179,'linear','none');

I3 = iradon(R2,0:179);
I4 = iradon(R2,0:179,'linear','none');

figure
subplot(1,2,1)
imshow(I1,[])
title('Filtered Backprojection')
subplot(1,2,2)
imshow(I3,[])
title('My filtered Backprojection')