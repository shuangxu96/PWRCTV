%% =========================== Frist part notes ===========================
% ADMM algorithm: tensor denosing
%    solve the following problem 
%           (SSTV regularized low rank tucker decomposition problem)
%    
%      The approximate model:  
%                                min tau*||F||_1 + lambda*||E||_1
%                             s.t. D = X+E, X = Z,DZ = F
%                                  X = Core*U1*U2*U3,rank(Core_i)=ri 
%      where,D is SSTV difference operator in the above model
%
%      The lagrange function is :
%       tau*||F||_1 + lambda*||E||_1 + <M1, D-X-E> + <M2,X-Z> +<Gamma,DZ - F>
%           + beta/2*( ||D-X-E||_F^2 + ||X-Z||_F^2 + ||DZ - F||_F^2 )
%
%% =========================== Second part notes===========================

% Reference paper: Hyperspectral Image Restoration via Total Variation 
%                      Regularized Low-rank Tensor Decomposition
% Author: Yao Wang, Jiangjun Peng
% E-mail addresses: andrew.pengjj@gmail.com
% -------------------------------------------------------------------------

%% =========================== Thrid part notes =========================== 
% INPUT:
%   Noi:     noisy 3-D image of size M*N*p normalized to [0,1] band by band
%   tau:     the trade-off parameter (recommended value 1)              
%   lambda:  sparse noise coefficient
%   rank:    rank constraint,[0.8*M,0.8*N,r]
% OUTPUT:
%  clean_iamge:   3-D denoised image
%  S:             The noise term
%  out_value:     MPSNR and MSSIM and ERGAS valuses of each iteration 
%  ========================================================================
function [clean_image,S,out_value,time] = LTNTV(Noi, tau,lambda,rank,beta)
tic
sizeD           = size(Noi);

normD           = norm(Noi(:)); 
n               = prod(sizeD);
maxIter         = 40;
epsilon         = 1e-6;  
mu1            = 0.01;             % The ascending multiplier value
mu2            = 0.01;
mu3 = 0.01;

out_value       = [];
out_value.SSIM  = [];
out_value.PSNR  = [];
out_value.ERGAS = [];

h               = sizeD(1);
w               = sizeD(2);
d               = sizeD(3);
%% 
% Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
% Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
% Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_x   = beta(1)^2*( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = beta(2)^2*( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z   = beta(3)^2*( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
determ  =  Eny_x + Eny_y + Eny_z;

%%  Initialization 
X               = zeros(sizeD);      % X : The clean image
Z               = X;                % Z : auxiliary variable for X
S               = zeros(sizeD);     % S : sparse noise 
F               = zeros(3*n,1);     % F : auxiliary variable for tv
C           = F;                % The multiplier for DZ-F
B              = zeros(size(Noi)); % The multiplier for          
A              = B;

%% main loop

for iter = 1: maxIter
    preX       = X;
    %% - update Core and U_i and X
    temp       = (mu1*(Z-A/mu1)+mu2*(Noi-S+B/mu2))/(mu1+mu2);
    X          = tucker_hooi(temp,rank);
    
    %% - update Z
%     diffT_p  = diffT3( mu3*F - C, sizeD );
    diffT_p  = diffT3_weight( mu3*F - C, sizeD,beta );
    numer1   = reshape( diffT_p + mu1*X(:) + A(:), sizeD);
    z        = real( ifftn( fftn(numer1) ./ (mu3*determ + mu1) ) );
    Z        = reshape(z,sizeD);
    
    %% - update F
%     diff_Z     = diff3(Z(:), sizeD);
    diff_Z     = diff3_weight(Z(:), sizeD,beta); 
    F          = prox_half(diff_Z+ C/mu3, tau/mu3 );  
   
    %% - update S 
    S          = softthre(Noi-X+B/mu2,lambda/mu2);% sparse
    
    %% - update M
    B         = B + mu2*(Noi-X-S);
    A         = A + mu1*(X-Z); 
    C         = C + mu3*(diff_Z-F);            
    mu1       = min(mu1 * 1.5,1e6); 
    mu2       = min(mu2 * 1.5,1e6); 
    mu3       = min(mu3 * 1.5,1e6); 
    
    %% compute the error
    errList    = norm(X(:)-preX(:)) / normD;
    fprintf('LRTDTV: iterations = %d   difference=%f\n', iter, errList);
    if errList < epsilon
        break;  
    end 
    %% output SSIM and PSNR values of each step
%     load simu_indian
%     OriData3 = simu_indian;
%     [out_value.PSNR(iter),out_value.SSIM(iter),out_value.ERGAS(iter)]=msqia(OriData3,X);
end
%% the final clean image
clean_image = X;
fprintf('LRTDTV ends: total iterations = %d,difference=%f\n\n', iter, errList);
toc
time=toc; 
end

function z = prox_half(x,gamma)
z = (2/3)*x.*(abs(x)>(gamma^(2/3).*54^(1/3)/4)).*(1+cos(2*pi/3-2*acos((abs(x)/3).^(-1.5)*gamma/8)/3));
end