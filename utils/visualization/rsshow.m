function A = rsshow(I, scale, ignore_value)
% Remote sensing image enhancement for visualization.
%
% Usage:
% display an image: rsshow(I, 0.05)
% write an image: A = rsshow(I, 0.05); imwrite(A, 'output.jpg')
%
% If the code is used in your scientific research, please cite the paper.
% [1] Shuang Xu, Xiangyong Cao, Jiangjun Peng, Qiao Ke, Cong Ma and Deyu
% Meng. Hyperspectral Image Denoising by Asymmetric Noise Modeling. IEEE
% TGRS, 2023.

if nargin==1
    scale = 0.005;
    ignore_value = NaN;
elseif nargin==2
    ignore_value = NaN;
end
I = double(I);

if size(I,3)>=3
    C = size(I,3);
    band = [C, uint8(C*0.5), uint8(C*0.1)];
    band(band<1) = 1;
    I = I(:,:,band);
end

Iq = I;
if ~isnan(ignore_value)
    Iq(Iq==ignore_value) = nan;
end

if ismatrix(I)
    q = quantile(Iq(:),[scale, 1-scale]);
    [low, high] = deal(q(1),q(2));
    I(I>high) = high;
    I(I<low) = low;
    I = (I-low)/(high-low);
else
    for i=1:size(I,3)
        temp = Iq(:,:,i);
        q = quantile(temp(:),[scale, 1-scale]);
        [low, high] = deal(q(1),q(2));
        temp = I(:,:,i);
        temp(temp>high) = high;
        temp(temp<low) = low;
        temp = (temp-low)/(high-low);
        I(:,:,i) = temp;
    end
end

if nargout==1
    A=I;
else
    imshow(I)
end
end