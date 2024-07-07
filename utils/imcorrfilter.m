function rho_XY = imcorrfilter(X,Y,w1,w2)
% % 创建一个3x3的局部区域
% h = ones(w1,w2)/w1/w2;
% 
% % 使用filter2计算图像X和Y的局部均值
% mean_X = filter2(h, X);
% mean_Y = filter2(h, Y);
% 
% % 计算图像X和Y的局部方差
% var_X = filter2(h, (X - mean_X).^2);
% var_Y = filter2(h, (Y - mean_Y).^2);
% 
% % 计算图像X和Y的局部协方差
% cov_XY = filter2(h, (X - mean_X).*(Y - mean_Y));
% 
% % 计算局部相关系数
% rho_XY = cov_XY./(sqrt(var_X).*sqrt(var_Y));

% 使用filter2计算图像X和Y的局部均值
mean_X = boxfilter(X, w1, w2);
mean_Y = boxfilter(Y, w1, w2);

% 计算图像X和Y的局部方差
var_X = boxfilter((X - mean_X).^2, w1, w2);
var_Y = boxfilter((Y - mean_Y).^2, w1, w2);

% 计算图像X和Y的局部协方差
cov_XY = boxfilter((X - mean_X).*(Y - mean_Y), w1, w2);

% 计算局部相关系数
rho_XY = cov_XY./(sqrt(var_X).*sqrt(var_Y));

end




function imDst = boxfilter(imSrc, r1, r2)

%   BOXFILTER   O(1) time box filtering using cumulative sum
%
%   - Definition imDst(x, y)=sum(sum(imSrc(x-r:x+r,y-r:y+r)));
%   - Running time independent of r; 
%   - Equivalent to the function: colfilt(imSrc, [r1, r2], 'sliding', @mean);
%   - But much faster.
% r1 = (w1-1)/2;
% r2 = (w2-1)/2;
[hei, wid, ~] = size(imSrc);
imDst = zeros(size(imSrc));

%cumulative sum over Y axis
imCum = cumsum(imSrc, 1);
%difference over Y axis
imDst(1:r1+1, :, :) = imCum(1+r1:2*r1+1, :, :);
imDst(r1+2:hei-r1, :, :) = imCum(2*r1+2:hei, :, :) - imCum(1:hei-2*r1-1, :, :);
imDst(hei-r1+1:hei, :, :) = repmat(imCum(hei, :, :), [r1, 1]) - imCum(hei-2*r1:hei-r1-1, :, :);

%cumulative sum over X axis
imCum = cumsum(imDst, 2);
%difference over Y axis
imDst(:, 1:r2+1, :) = imCum(:, 1+r2:2*r2+1, :);
imDst(:, r2+2:wid-r2, :) = imCum(:, 2*r2+2:wid, :) - imCum(:, 1:wid-2*r2-1, :);
imDst(:, wid-r2+1:wid, :) = repmat(imCum(:, wid, :), [1, r2]) - imCum(:, wid-2*r2:wid-r2-1, :);
%average
imDst = imDst/(2*r1+1)/(2*r2+1);
end

