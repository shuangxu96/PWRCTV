clear all;clc;	
load('Simu_indian.mat')
%load('pure_DCmall.mat')
Ohsi = Normalize(Ori_H);
Nhsi      = Ohsi;
[M,N,p]   = size(Ohsi);

noiselevel = 0.075*rand(p,1);
ratio = 0.15*rand(p,1);
%% Gaussian noise
for i = 1:p
     Nhsi(:,:,i)=Ohsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
for i = 1:p
     Nhsi(:,:,i)=imnoise(Nhsi(:,:,i),'salt & pepper',ratio(i));
end

r=13;
beta = 50;
lambda = 1;% 5,0.5
tau = [0.8,0.8];% need to fine tune
k=4;
tic;
output_image = MoG_RBTV(Nhsi, beta,lambda, tau, r, k);
time = toc;
[mpsnr,mssim,ergas]=msqia(Ohsi, output_image)