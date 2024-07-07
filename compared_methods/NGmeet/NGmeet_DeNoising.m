function [E_Img]= NGmeet_DeNoising( N_Img, O_Img, Par )
% reference paper: Non-local Meets Global: An Integrated Paradigm for Hyperspectral Denoising
% oriData3_noise  N_Img        Input noisy 3-D image
% OriData3        O_Img        Reference image, for the PSNR computing of each step
% k_subspace  Par.k_subspace   The initial spectral rank, can be estimated by HySime
% delta=2                       iteration of k_subspace
% for pavia city and CAVE k_subspace = 5+2*(iter-1);
delta =2;
E_Img            = N_Img;                                                         % Estimated Image
[Height, Width, Band]  = size(E_Img);  
N = Height*Width;
TotalPatNum      = (Height-Par.patsize+1)*(Width-Par.patsize+1);         % Total Patch Number in the image
Average          = mean(N_Img,3);                      % Calculate the average band for fast spatial non-local searching
[Neighbor_arr, Num_arr, Self_arr] =	NeighborIndex(Average, Par);   
% PreCompute all the patch index in the searching window 

for iter = 1 : Par.Iter 
%First step: spectral dimension reduction 
   k_subspace = Par.k_subspace+delta*(iter-1);
   Y = reshape(E_Img, N, Band)';
   E_Img1 = E_Img;
%    [w Rw] = estNoise(Y,'additive');
%    Rw_ori = Rw;
%    Y = sqrt(inv(Rw_ori))*Y;
%    img_ori = reshape(Y', Height, Width, Band);
%    [w Rw] = estNoise(Y,'additive');
%    [~, E]=hysime(Y,w,Rw);
   [E,~,~]= svd(Y,'econ');
    E=E(:,1:k_subspace);

    E_Img = reshape((E'*Y)', Height,Width, k_subspace);
    N_Img1 = reshape((E'*reshape(N_Img, N, Band)')', Height,Width, k_subspace); %%% add change N_Img1 as N_Img 
    Band1=k_subspace;

% %non-local patch grouping and noise estimation
    Average             =   mean(E_Img,3);
    [CurPat, Mat, Sigma_arr]	=	Cub2Patch( E_Img, N_Img1, Average, Par );

    if (mod(iter-1,2)==0)
        Par.patnum = Par.patnum - 10;                                          % Lower Noise level, less NL patches
        NL_mat  =  Block_matching(Mat, Par, Neighbor_arr, Num_arr, Self_arr);  % Caculate Non-local similar patches for each
        if(iter==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr))*sqrt(k_subspace/Band);                      % First Iteration use the input noise parameter
        end
    end
%     time2=toc
% non-local low-rank denoising
    [Spa_EPat, Spa_W]    =  NLPatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par); 
% reconstruct patches to 3-D image
    [Spa_Img, Spa_Wei]   =  Patch2Cub( Spa_EPat, Spa_W, Par.patsize, Height, Width, Band1 );       % Patch to Cubic
    E_Img = Spa_Img./Spa_Wei;
    E_Img = reshape(reshape(E_Img, Height*Width, k_subspace)*E',Height,Width, Band);

% time3 = toc
% estimation, can be ignored to speed up
%     [PSNR,SSIM,~,~] = evaluate(O_Img/255,E_Img/255,Height,Width);PSNR = mean(PSNR);SSIM = mean(SSIM);
%    PSNR = mean(0);SSIM = mean(0);
    PSNR  = csnr( O_Img, E_Img, 0, 0 );
    SSIM  = cal_ssim( O_Img, E_Img, 0, 0 );
    fprintf( 'Iter = %2.3f, PSNR = %2.2f, SSIM = %2.3f, NoiseLevel = %2.3f \n', iter, PSNR, SSIM, sum(Sigma_arr)/TotalPatNum);
    if iter<Par.Iter
    E_Img = 0.1*N_Img+0.9*E_Img;
    else
    end
end
end





function s=csnr(A,B,row,col)

[n,m,ch]=size(A);
summa = 0;
if ch==1
   e=A-B;
   e=e(row+1:n-row,col+1:m-col);
   me=mean(mean(e.^2));
   s=10*log10(255^2/me);
else
    for i=1:ch
        e=A-B;
        e=e(row+1:n-row,col+1:m-col,i);
        mse = mean(mean(e.^2));
        s  = 10*log10(255^2/mse);
        summa = summa + s;
    end
        s = summa/ch;
end
end

function ssim  =  cal_ssim( im1, im2, b_row, b_col )

[h w ch]  =  size( im1 );
ssim  = 0;
if ch==1
    ssim   = ssim_index( im1(b_row+1:h-b_row, b_col+1:w-b_col), im2( b_row+1:h-b_row, b_col+1:w-b_col) );
else
    for i = 1:ch
        ssim   = ssim + ssim_index( im1(b_row+1:h-b_row, b_col+1:w-b_col, i), im2( b_row+1:h-b_row, b_col+1:w-b_col, i) );
    end
    ssim   =  ssim/ch;
end
end




function [mssim, ssim_map] = ssim_index(img1, img2, K, window, L)

%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author was with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University, USA. He is
%currently with Department of Electrical and Computer Engineering,
%University of Waterloo, Canada.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 4, Apr. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim ssim_map] = ssim_index(img1, img2, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================


if (nargin < 2 | nargin > 5)
   mssim = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   mssim = -Inf;
   ssim_map = -Inf;
   return;
end

[M N] = size(img1);

if (nargin == 2)
   if ((M < 11) | (N < 11))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;								      %
   L = 255;                                  %
end

if (nargin == 3)
   if ((M < 11) | (N < 11))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   [H W] = size(window);
   if ((H*W) < 4 | (H > M) | (W > N))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 5)
   [H W] = size(window);
   if ((H*W) < 4 | (H > M) | (W > N))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
	denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

end