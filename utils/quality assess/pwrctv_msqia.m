function [mpsnr,mssim,ergas,sam]=msqia(imagery1, imagery2)
% Evaluates the quality assessment indices for two HSIs.
%
% Syntax:
%   [mpsnr, mssim,ergas ] = MSIQA(imagery1, imagery2)
% Input:
%   imagery1 - the reference HSI data array
%   imagery2 - the target HSI data array
%   NOTE: MSI data array  is a M*N*K array for imagery with M*N spatial
%	pixels, K bands and DYNAMIC RANGE [0,1];
[M,N,p]  = size(imagery1);
psnrvector=zeros(1,p);
for i=1:1:p
    J=255*imagery1(:,:,i);
    I=255*imagery2(:,:,i);
    psnrvector(i)=PSNR_c(J,I,M,N);
end 
mpsnr = mean(psnrvector);

SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*imagery1(:,:,i);
    I=255*imagery2(:,:,i); 
    [ SSIMvector(i),~] = ssim(J,I);
end
mssim=mean(SSIMvector);

ergas = ErrRelGlobAdimSyn(255*imagery1, 255*imagery2);

sam = SAM(255*imagery1, 255*imagery2);
end

function ergas = ErrRelGlobAdimSyn(imagery1, imagery2)

%==========================================================================
% Evaluates Erreur Relative Globale Adimensionnelle de Synthèse (ERGAS, 
%   dimensionless global relative error of synthesis)[1] for two MSIs.
%
% Syntax:
%   ergas = ErrRelGlobAdimSyn(imagery1, imagery2)
%
% Input:
%   imagery1 - the reference MSI data array
%   imagery2 - the target MSI data array
% NOTE: MSI data array  is a M*N*K array for imagery with M*N spatial
%	pixels, K bands and DYNAMIC RANGE [0, 255]. If imagery1 and imagery2
%	have different size, the larger one will be truncated to fit the
%	smaller one.
%
% [1] L. Wald, "Data Fusion: Definitions and Architectures—Fusion of
%     Images of Different Spatial Resolutions", Les Presses, Ecole des
%     Mines de Paris, Paris, France (2002) 200 pp.
%
% See also StructureSIM, FeatureSIM and SpectAngMapper
%
% by Yi Peng
%==========================================================================

[m, n, k] = size(imagery1);
[mm, nn, kk] = size(imagery2);
m = min(m, mm);
n = min(n, nn);
k = min(k, kk);
imagery1 = imagery1(1:m, 1:n, 1:k);
imagery2 = imagery2(1:m, 1:n, 1:k);

ergas = 0;
for i = 1:k
    ergas = ergas + mse(imagery1(:, :, i) - imagery2(:, :, i))/mean2(imagery1(:, :, i));
end
ergas = 100*sqrt(ergas/k);
end


function result = PSNR_c(ReferBuffer,UnReferBuffer,lHeight,lWidth)

    result = 0;
	for j=1:lWidth*lHeight
		temp = ReferBuffer(j)-UnReferBuffer(j);
		result = result + double(temp*temp);
    end;
    
	if (result==0)
		result =100;
	else 
		result = 10*log10(255*255/result*lWidth*lHeight);
    end;
end



function [mssim, ssim_map] = ssim(img1, img2, K, window, L)

% ========================================================================
% SSIM Index with automatic downsampling, Version 1.0
% Copyright(c) 2009 Zhou Wang
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for calculating the
% Structural SIMilarity (SSIM) index between two images
%
% Please refer to the following paper and the website with suggested usage
%
% Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
% quality assessment: From error visibility to structural similarity,"
% IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,
% Apr. 2004.
%
% http://www.ece.uwaterloo.ca/~z70wang/research/ssim/
%
% Note: This program is different from ssim_index.m, where no automatic
% downsampling is performed. (downsampling was done in the above paper
% and was described as suggested usage in the above website.)
%
% Kindly report any suggestions or corrections to zhouwang@ieee.org
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
%            has a smaller size than the input images. The actual size
%            depends on the window size and the downsampling factor.
%
%Basic Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim, ssim_map] = ssim(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim, ssim_map] = ssim(img1, img2, K, window, L);
%
%Visualize the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%========================================================================


if (nargin < 2 || nargin > 5)
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
   if ((M < 11) || (N < 11))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   K(1) = 0.01;					% default settings
   K(2) = 0.03;					%
   L = 255;                                     %
end

if (nargin == 3)
   if ((M < 11) || (N < 11))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
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
   if ((H*W) < 4 || (H > M) || (W > N))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
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
   if ((H*W) < 4 || (H > M) || (W > N))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
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


img1 = double(img1);
img2 = double(img2);

% automatic downsampling
f = max(1,round(min(M,N)/256));
%downsampling by f
%use a simple low-pass filter 
if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 && C2 > 0)
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

% return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Spectral Angle Mapper (SAM).
% 
% Interface:
%           [SAM_index,SAM_map] = SAM(I1,I2)
%
% Inputs:
%           I1:         First multispectral image;
%           I2:         Second multispectral image.
% 
% Outputs:
%           SAM_index:  SAM index;
%           SAM_map:    Image of SAM values.
% 
% References:
%           [Yuhas92]   R. H. Yuhas, A. F. H. Goetz, and J. W. Boardman, "Discrimination among semi-arid landscape endmembers using the Spectral Angle Mapper (SAM) algorithm," 
%                       in Proceeding Summaries 3rd Annual JPL Airborne Geoscience Workshop, 1992, pp. 147?49.
%           [Vivone20]  G. Vivone, M. Dalla Mura, A. Garzelli, R. Restaino, G. Scarpa, M.O. Ulfarsson, L. Alparone, and J. Chanussot, "A New Benchmark Based on Recent Advances in Multispectral Pansharpening: Revisiting pansharpening with classical and emerging pansharpening methods", 
%                       IEEE Geoscience and Remote Sensing Magazine, doi: 10.1109/MGRS.2020.3019315.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SAM_index,SAM_map] = SAM(I1,I2)

[M,N,~] = size(I2);
prod_scal = dot(I1,I2,3); 
norm_orig = dot(I1,I1,3);
norm_fusa = dot(I2,I2,3);
prod_norm = sqrt(norm_orig.*norm_fusa);
prod_map = prod_norm;
prod_map(prod_map==0)=eps;
SAM_map = acos(prod_scal./prod_map);
prod_scal = reshape(prod_scal,M*N,1);
prod_norm = reshape(prod_norm, M*N,1);
z=find(prod_norm==0);
prod_scal(z)=[];prod_norm(z)=[];
angolo = sum(sum(acos(prod_scal./prod_norm)))/(size(prod_norm,1));
SAM_index = real(angolo)*180/pi;

end