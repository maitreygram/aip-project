function theta = ISTA(y,A,reg,con_param,alpha,mode)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   Mode == 0 for 2d-DCT basis and Mode == 1 for 2d-Haar wavelets where A =
%   phi
% con_param = 0.1;
[~,col] = size(A);
% alpha = max(eig(A'*A));
theta = zeros(col,1); % initialize theta    
while 1
    if mode == 0
        theta_temp = wthresh(theta + (1/alpha)*A'*(y - A*theta),'s',reg/(2*alpha));
    end
    if mode == 1
        theta_temp = wthresh(theta + (1/alpha)*dwt2on1d(A'*(y - A*idwt2on1d(theta))),'s',reg/(2*alpha));
    end
    res = norm(theta_temp  -  theta);
    theta = theta_temp;
    if ~(res > con_param)
        break;
    end
end

end

