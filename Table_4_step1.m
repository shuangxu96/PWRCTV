%% load data
load('data\Urban\Urban_F210.mat')
load('data\Urban\ikonos_spec_resp.mat')
wvl = load('data\Urban\URBAN.wvl');

Nhsi = reshape(Y, [nBand,nRow,nCol]);
Nhsi = permute(Nhsi, [2,3,1]);

Chsi = Nhsi(:,:,wvl(:,4)==1 & wvl(:,2)>=350 & wvl(:,2)<=1035);
working_wvl = wvl(wvl(:,4)==1 & wvl(:,2)>=350 & wvl(:,2)<=1035,2);
srf = zeros(size(Chsi,3),1);
pan = zeros(size(Chsi,1),size(Chsi,2));
for i=1:length(working_wvl)
    temp_wvl = working_wvl(i);
    [~,index] = min(abs(temp_wvl - ikonos_sp(:,1)));
    srf(i) = ikonos_sp(index,2);
    pan = pan + srf(i)*Chsi(:,:,i);
end
pan = pan/sum(srf);
plot(working_wvl,srf)
title('Spectral Response Function')

%% Preprocess
Nhsi = Nhsi/max(Nhsi(:));
pan  = pan/max(pan(:));
[M,N,B] = size(Nhsi);

%% set savepath
savepath = fullfile('result','urban');
mkdir(savepath)

%% save simulations
save([savepath,'\Noisy.mat'],'Nhsi','pan')

%% TCTV
% parameters
opts = [];
opts.rho = 1.25;
opts.directions = [1,2,3];
% run algorithm
disp('--------------- Run TCTV --------------- ')
tic
output = TCTV_TRPCA(Nhsi, opts);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\TCTV.mat'],'output')

%% LMHTV
% parameters
r = 4;
tau = 0.0005;
lambda = 10/sqrt(M*N);
% run algorithm
disp('--------------- Run LMHTV --------------- ')
tic
output = LMHTV(Nhsi, tau, lambda, r, [1,1,0]);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\LMHTV.mat'],'output')

%% LTHTV
% parameters
tau = 0.04;
lambda = 1000/sqrt(M*N);
ten_rank = [ceil(M*0.7),ceil(N*0.7),r];
% run algorithm
disp('--------------- Run LTHTV --------------- ')
tic
output = LTHTV(Nhsi, tau,lambda,ten_rank,[1,1,0]);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\LTHTV.mat'],'output')

%% LRTV
% parameters
if size(Nhsi,3) > 100
    tau = 0.015;
    lambda = 20/sqrt(M*N);
    rank = 10;
else
    tau = 0.01;
    lambda = 10/sqrt(M*N);
    rank = 5;
end
% run algorithm
disp('--------------- Run LRTV --------------- ')
tic
output = LRTV(Nhsi, tau, lambda, rank);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\LRTV.mat'],'output')



%% RCTV
% parameters
beta = 200;
lambda = 1;
tau = 0.8*[1,1];
r = 4;
q = 10;
% run algorithm
disp('--------------- Run RCTV --------------- ')
tic
output = MoG_RBTV(Nhsi, beta, lambda, tau, r);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\RCTV.mat'],'output')

%% BALMF
% parameters
r = 4;
% run algorithm
disp('--------------- Run BALMF --------------- ')
tic
output = BALMF(Nhsi, r);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\BALMF.mat'],'output')



%% CTV
% parameters
opts = [];
opts.rho = 1.5;
% run algorithm
disp('--------------- Run CTV --------------- ')
tic
output = ctv_rpca(Nhsi, opts);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\CTV.mat'],'output')


%% PWRCTV
% parameters
beta = 100;
lambda = 1;
tau = 0.5*[1,1];
r = 3;
q = 2;
% run algorithm
disp('--------------- Run PWRCTV --------------- ')
tic
output = PWRCTV(Nhsi, pan, beta, lambda, tau, r, q);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\PWRCTV.mat'],'output')

%% TDL
Ohsi = im2double(output);
% parameters
noiselevel = std(reshape(Ohsi-Nhsi, [M*N,B]));
vstbmtf_params.peak_value = 1;
vstbmtf_params.nsigma = mean(noiselevel);
% run algorithm
disp('--------------- Run TDL --------------- ')
tic
output = TensorDL(Nhsi, vstbmtf_params);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\TDL.mat'],'output')

%% NGMeet
% parameters
noiselevel = std(reshape(Ohsi-Nhsi, [M*N,B]));
Par   = ParSetH(255*mean(noiselevel),B);
% run algorithm
disp('--------------- Run NGMeet --------------- ')
tic
output = NGmeet_DeNoising( 255*Nhsi, 255*Ohsi, Par);  %NGmeet denoisng function
elapsed_time = toc;
% save data
output = uint8(output);
save([savepath,'\NGMeet.mat'],'output')

%% WNLRATV
% parameters
noise     = reshape(Nhsi - Ohsi, M*N,B);
Sigma_ratio  = std(noise(:));
initial_rank  = 3;
Rank = 6;
ModelPar.alpha = 30;
ModelPar.belta = 1;
ModelPar.gamma = 0.08;
param   = SetParam_NWT(Nhsi, Sigma_ratio);
param.initial_rank = initial_rank;
param.maxiter = 15;
param.patnum        = 200;
param.lambda        = 2e-1;
[prior, model] = InitialPara( param,0,B);
% run algorithm
disp('--------------- Run WNLRATV --------------- ')
tic
output = WNLRATV2(Nhsi,Ohsi, Rank,ModelPar, param, model, prior);
elapsed_time = toc;
% save data
output = uint8(255*output);
save([savepath,'\WNLRATV.mat'],'output')

