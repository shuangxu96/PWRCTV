addpath('compared_methods\AWLP\')
addpath('utils\')
addpath('utils\quality assess')
%% load data
clc
clear
load('data\Milan\Milan.mat')

I_GT = im2double(hsi);
I_HSI_LR = imresize(I_GT, 1/ratio);
I_NHSI_LR = I_HSI_LR + randn(size(I_HSI_LR))*10/255;
I_NHSI_LR(I_NHSI_LR>1) = 1;
I_NHSI_LR(I_NHSI_LR<0) = 0;
I_PAN = im2double(pan);

ratio = 6;
metrics = zeros(3,4);

% pan-sharpening with noisy LR-HSI (Fig. 10c)
I_Fus_wN = AWLP(imresize(I_NHSI_LR,ratio), I_PAN,ratio);
[mpsnr,mssim,ergas,sam] = pwrctv_msqia(I_GT, I_Fus_wN);
metrics(1,1:4) = [mpsnr,mssim,ergas,sam];

% denoised LR-HSI + pan-sharpening (Fig. 10d)
beta = 100;
lambda = 1;
tau = 0.4*[1,1];
q = 10;
r = 4;
I_DNHSI_LR = PWRCTV(I_NHSI_LR, imresize(I_PAN,1/ratio), beta, lambda, tau, r, q);
I_Fus_DP = AWLP(imresize(I_DNHSI_LR,ratio), I_PAN,ratio);
[mpsnr,mssim,ergas,sam] = pwrctv_msqia(I_GT, I_Fus_DP);
metrics(2,1:4) = [mpsnr,mssim,ergas,sam];

% pan-sharpening + denoised LR-HSI (Fig. 10e)
beta = 100;
lambda = 1;
tau = 0.4*[1,1];
q = 10;
r = 4;
I_Fus_PD = PWRCTV(I_Fus_wN, I_PAN, beta, lambda, tau, r, q);
[mpsnr,mssim,ergas,sam] = pwrctv_msqia(I_GT, I_Fus_PD);
metrics(3,1:4) = [mpsnr,mssim,ergas,sam];

% save result
savepath = 'result\pansharpening';
mkdir(savepath)
save(fullfile(savepath,'Milan_iidGauss.mat'), 'I_NHSI_LR', 'I_Fus_DP', 'I_Fus_PD','I_Fus_wN')

rgb_index = [58,47,36];
imwrite(rsshow(I_GT(:,:,rgb_index)), fullfile(savepath,'Milan_Pansharpening_GT.jpg'))
imwrite(rsshow(I_NHSI_LR(:,:,rgb_index)), fullfile(savepath,'Milan_Pansharpening_LR_NHSI.jpg'))
imwrite(rsshow(I_Fus_wN(:,:,rgb_index)), fullfile(savepath,'Milan_Pansharpening_Fus_noisy.jpg'))
imwrite(rsshow(I_Fus_DP(:,:,rgb_index)), fullfile(savepath,'Milan_Pansharpening_Fus_D+P.jpg'))
imwrite(rsshow(I_Fus_PD(:,:,rgb_index)), fullfile(savepath,'Milan_Pansharpening_Fus_P+D.jpg'))
imwrite(rsshow(I_PAN), fullfile(savepath,'Milan_Pansharpening_PAN.jpg'))