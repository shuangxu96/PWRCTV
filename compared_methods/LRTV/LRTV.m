function [ output_image output_sparse ] = LRTV_accelerate(Y, tau,lambda,r)
% Solve problem
% solve the following problem (TV regularized LR and MC problem)
%  argmin   ||L||_nuclear + tao * ||X||_TV+lanbda*||E||_1
%                       s.t. X = L D = L + E and rank(L)<=r;
%     via IALM
% -------------------------------------------------------------------------------------------------------------
% Reference paper: W. He, H. Zhang, L. Zhang, and  H. Shen, “Total-Variation-Regularized
% Low-Rank Matrix Factorization for Hyperspectral Image Restoration,” IEEE Trans. Geosci. Remote Sens.,
% vol. 54, pp. 178-188, Jan. 2016.
% Author: Wei He (November, 2014)
% E-mail addresses:(weihe1990@whu.edu.cn)
% --------------------------------------------------INPUT-----------------------------------------------------
%  oriData3_noise                            noisy 3-D image of size sizeY(1)*N*p normalized to [0,1] band by band
%  tao                                       (recommended value 0.01)
%  lambda
%  G1(omega)                                 all one matrix of size (sizeY(1)*N)*p
%  G0(omega~)                                all zero matrix of size (sizeY(1)*N)*p
%  r                                         rank constraint
% --------------------------------------------------OUTPUT-----------------------------------------------------
%  output_iamge                              3-D denoised image
%  out_value                                 MPSNR and MSSIM valuses of each iteration
% -------------------------------------------------------------------------------------------------------------
% Note: the parameters G0 and G1 do not have any effect. these two
% parameters are used to solve the impainting problem with the location of
% missing pixels to be known.
% -------------------------------------------------------------------------------------------------------------
sizeY = size(Y);
sizeD = [sizeY(1)*sizeY(2), sizeY(3)];
D = reshape(Y, sizeD); clear Y

tol = 1e-5;
maxIter = 100;
rho = 1.25;

% Initialize mu
normD   = norm(D,'fro');
[~,norm_two,~] = svdsecon(D, 1);
norm_inf = max(abs(D(:)))/lambda;
dual_norm = max(norm_two, norm_inf);
mu = 1.25/dual_norm; % this one can be tuned
max_mu1 = mu * 1e8;
mu1 = mu; mu2 = mu; mu3 = mu;
%%
debug=1;
%% Initializing optimization variables
% intialize
% L = rand(sizeD);
[uu,ss,vv] = svdsecon(D,6);
L = uu*ss*vv';
f = reshape(L, sizeY);
S = zeros(sizeD);

M1 = zeros(sizeD);
M2 = zeros(sizeD);

% F1 = zeros(sizeY); F2 = zeros(sizeY); 
M31 = zeros(sizeY); M32 = zeros(sizeY); 

% for the 3D-TV norm
% define operators
Eny_x   = ( abs(psf2otf([+1; -1], sizeY)) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], sizeY)) ).^2  ;
eigDtD  =  Eny_x + Eny_y;

[diff,diff_t]      = defDDt();


% main loop
iter = 0;
tic
while iter<maxIter
    iter = iter + 1;
    % update sizeY(1) = [u1 u2 u3]
    [Df1, Df2] = diff(f);
    F1 = softthre(Df1+(1/mu3)*M31, tau/mu3);
    F2 = softthre(Df2+(1/mu3)*M32, tau/mu3);

    % Update X
    numer1 = reshape(L-M1/mu1, sizeY);
    numer2 = diff_t(F1-(1/mu3)*M31,  F2-(1/mu3)*M32);
    rhs    = fftn( (mu1/mu3)*numer1 + numer2 );
    lhs    = (mu1/mu3) + eigDtD;
    f      = real(ifftn(rhs./lhs));
    X      = reshape(f,sizeD);

    % update L
    L = prox_nn((mu1*X +mu2* (D-S) + (M1+M2))/(mu1+mu2), 1/(mu1+mu2), r);

    % update E
    S = softthre(D - L + M2/mu2, lambda/mu2);

    % stop criterion
    stopC = norm(D-L-S,'fro')/normD;
    converge_flag = stopC<tol;
    if debug && (mod(iter,5)==0 || iter==1 || iter==maxIter || converge_flag)
        disp(['iter ' num2str(iter) ',mu=' num2str(mu1,'%2.1e')  ...
            ',Y-X-S=' num2str(stopC,'%2.3e')]);
    end

    if  stopC<tol
        break;
    else
        M1 = M1 + mu1*(X-L);
        M2 = M2 + mu2*(D-L-S);
        
        M31   = M31 + mu3*(Df1-F1  );
        M32   = M32 + mu3*(Df2-F2  );

        mu1 = min(max_mu1,mu1*rho);
        mu2 = min(max_mu1,mu2*rho);
        mu3 = min(max_mu1,mu3*rho);

    end


end
toc
output_image = reshape(L,sizeY);
output_sparse = reshape(S,sizeY);
end

function [D,Dt] = defDDt()
D  = @(U) ForwardD(U);
Dt = @(X,Y) Dive(X,Y);
end

function [Dux,Duy] = ForwardD(U)

Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];

end

function DtXYZ = Dive(X,Y)

DtXYZ = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
DtXYZ = DtXYZ + [Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];

end

