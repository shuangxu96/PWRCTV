function [ output_image,output_sparse,output_noise] = LMHTV_Lq(Y_tensor, tau, lambda, r, beta, shrink_mode)
% Solve problem
% solve the following problem (TV regularized LR and MC problem)
%  argmin   ||L||_nuclear + tao * ||X||_TV+lanbda*||S||_1
%                       s.t. X = L D = L + S and rank(L)<=r;
%     via IALM
% -------------------------------------------------------------------------------------------------------------
% Reference paper: W. He, H. Zhang, L. Zhang, and  H. Shen, ¡°Total-Variation-Regularized
% Low-Rank Matrix Factorization for Hyperspectral Image Restoration,¡± IEEE Trans. Geosci. Remote Sens.,
% vol. 54, pp. 178-188, Jan. 2016.
% Author: Wei He (November, 2014)
% E-mail addresses:(weihe1990@whu.edu.cn)
% --------------------------------------------------INPUT-----------------------------------------------------
%  Y_tensor                            noisy 3-D image of size M*N*p normalized to [0,1] band by band
%  tao                                       (recommended value 0.01)
%  lambda
%  G1(omega)                                 all one matrix of size (M*N)*p
%  G0(omega~)                                all zero matrix of size (M*N)*p
%  r                                         rank constraint
% --------------------------------------------------OUTPUT-----------------------------------------------------
%  output_iamge                              3-D denoised image
%  out_value                                 MPSNR and MSSIM valuses of each iteration
% -------------------------------------------------------------------------------------------------------------
% Note: the parameters G0 and G1 do not have any effect. these two
% parameters are used to solve the impainting problem with the location of
% missing pixels to be known.
% -------------------------------------------------------------------------------------------------------------

%% Preprocessing Data
[M,N,p] = size(Y_tensor);
Y = reshape(Y_tensor, [M*N,p]);
d_norm = norm(Y, 'fro');

%% Input variables
if nargin==1
    tau = 0.05;
    lambda = 10/sqrt(M*N);
    r = 10;
    beta = [1,1,0.5];
    shrink_mode = 'L1/2';
elseif nargin==2
    lambda = 10/sqrt(M*N);
    r = 10;
    beta = [1,1,0.5];
    shrink_mode = 'L1/2';
elseif nargin==3
    r = 10;
    beta = [1,1,0.5];
    shrink_mode = 'L1/2';
elseif nargin==4
    beta = [1,1,0.5];
    shrink_mode = 'L1/2';
elseif nargin==5
    shrink_mode = 'L1/2';
end

%% Parameters
sv      = 10;
mu1     = 1e-2; % The ascending multiplier value
mu2     = 1e-2;
mu3     = 1e-2;
max_mu  = 1e6;  % max value for mu
rho     = 1.5;  % ascending factor
tol1    = 1e-6; % tolerance for converge
tol2    = 1e-6;
maxIter = 30;

%% Initializing Variables
Z = zeros(size(Y));  % auxiliary variable for X
S = sparse(Z); % sparse noise

A = zeros(size(Y));  % The multiplier for Z-X
B = zeros(size(Y));  % The multiplier for Y-X-S

F1 = zeros(M,N,p); F2 = zeros(M,N,p); F3 = zeros(M,N,p); % auxiliary variable for DZ
C1 = zeros(M,N,p); C2 = zeros(M,N,p); C3 = zeros(M,N,p); % The multiplier for DZ-F

%% 3D-TV
% beta = [1 1 0.5];
[D,Dt] = defDDt(beta);
eigDtD = abs(beta(1)*fftn([1 -1],  [M N p])).^2 + abs(beta(2)*fftn([1 -1]', [M N p])).^2;
if p>1
    d_tmp(1,1,1)= 1; d_tmp(1,1,2)= -1;
    eigEtE  = abs(beta(3)*fftn(d_tmp, [M N p])).^2;
else
    eigEtE = 0;
end


%% main loop
iter = 0;
tic
while iter<maxIter
    iter = iter + 1;
    %%     Update X
    X = ( mu1*Z + mu2*(Y-S) + (A+B) ) / (mu1+mu2);
    [X,sv] = prox_nuclear(X, mu1, mu2, p, sv, r);
    
    %%       Update Z
    numer1 = reshape(X-A/mu1,M,N,p);
    numer2 = Dt(F1-(1/mu3)*C1,  F2-(1/mu3)*C2, F3-(1/mu3)*C3);
    rhs  = fftn( (mu1/mu3)*numer1 + numer2 );
    lhs  = (mu1/mu3) + eigDtD + eigEtE;
    f    = real(ifftn(rhs./lhs));
    Z    = reshape(f,M*N,p);
    
    %%       update F1 F2 F3
    [DZ1,DZ2,DZ3] = D(f); % D(Z)
    %     F1 = prox_half(DZ1+(1/mu3)*C1, tau/mu3);
    %     F2 = prox_half(DZ2+(1/mu3)*C2, tau/mu3);
    %     if beta(3)==0
    %         F3 = 0;
    %     else
    %         F3 = prox_half(DZ3+(1/mu3)*C3, tau/mu3);
    %     end
    F1 = prox_sparse(DZ1+(1/mu3)*C1, tau/mu3, shrink_mode);
    F2 = prox_sparse(DZ2+(1/mu3)*C2, tau/mu3, shrink_mode);
    if beta(3)==0
        F3 = 0;
    else
        F3 = prox_sparse(DZ3+(1/mu3)*C3, tau/mu3, shrink_mode);
    end
    
    %% update S
    S = prox_L1(Y - X + B/mu2, lambda/mu2);
    
    %% stop criterion
    leq1 = Z - X;
    leq2 = Y -X -S ;
    stopC1 = max(max(abs(leq1)));
    stopC2 = norm(leq2, 'fro') / d_norm;
    disp(['iter ' num2str(iter) ',mu=' num2str(mu1,'%2.1e')  ...
        ',rank = ' num2str(rank(X))  ',stopALM=' num2str(stopC2,'%2.3e')...
        ',stopE=' num2str(stopC1,'%2.3e')]);
    
    if stopC1<tol1  && stopC2<tol2
        break;
    else
        A  = A + mu1*leq1;
        B  = B + mu2*leq2;
        C1 = C1 + mu3*(DZ1 - F1);
        C2 = C2 + mu3*(DZ2 - F2);
        C3 = C3 + mu3*(DZ3 - F3);
        
        mu1 = min(max_mu,mu1*rho);
        mu2 = min(max_mu,mu2*rho);
        mu3 = min(max_mu,mu3*rho);
        
    end
end
toc
output_image = reshape(X,[M,N,p]);
output_sparse = reshape(S,[M,N,p]);
output_noise = Y_tensor-output_image-output_sparse;
end

function [D,Dt] = defDDt(beta)
D  = @(U) ForwardD(U, beta);
Dt = @(X,Y,Z) Dive(X,Y,Z, beta);
end

function [Dux,Duy,Duz] = ForwardD(U, beta)
frames = size(U, 3);
Dux = beta(1)*[diff(U,1,2), U(:,1,:) - U(:,end,:)];
Duy = beta(2)*[diff(U,1,1); U(1,:,:) - U(end,:,:)];
Duz(:,:,1:frames-1) = beta(3)*diff(U,1,3);
Duz(:,:,frames)     = beta(3)*(U(:,:,1) - U(:,:,end));
end

function DtXYZ = Dive(X,Y,Z, beta)
frames = size(X, 3);
DtXYZ = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
DtXYZ = beta(1)*DtXYZ + beta(2)*[Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];
Tmp(:,:,1) = Z(:,:,end) - Z(:,:,1);
Tmp(:,:,2:frames) = -diff(Z,1,3);
DtXYZ = DtXYZ + beta(3)*Tmp;
end

function z = prox_sparse(x, gamma, shrink_mode)
switch shrink_mode
    case 'L0'
        z = prox_L0(x,gamma);
    case 'L1/2'
        z = prox_half(x,gamma);
    case 'L1'
        z = prox_L1(x,gamma);
end
end

function z = prox_half(x,gamma)
z = (2/3)*x.*(abs(x)>(gamma^(2/3).*54^(1/3)/4)).*(1+cos(2*pi/3-2*acos((abs(x)/3).^(-1.5)*gamma/8)/3));
end

function x = prox_L0(x,gamma)
x(abs(x)<gamma)=0;
end

function z = prox_L1(x,gamma)
z = sign(x).*max((abs(x)-gamma),0);
end

function [X,sv] = prox_nuclear(temp, mu1, mu2, p, sv, r)
if  choosvd(p,sv) ==1
    [U, sigma, V] = lansvd(temp, sv, 'L');
else
    [U,sigma,V] = svd(temp,'econ');
end
sigma = diag(sigma);
svp = min(length(find(sigma>1/(mu1+mu2))),r);
if svp<sv
    sv = min(svp + 1, p);
else
    sv = min(svp + round(0.05*p), p);
end
X = U(:, 1:svp) * diag(sigma(1:svp) - 1/(mu1+mu2)) * V(:, 1:svp)';
end