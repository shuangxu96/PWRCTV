function [output_image,U,V,MPSNR,MSSIM,ERGAS,SAM] = PWRCTV(Nhsi, Pan, beta, lambda, tau, r, q)
%% Initializing admm variables

if nargin < 2
    beta = 0.5;
end
if nargin < 3
    lambda = 1;
end
if nargin < 4
    tau = [0.01,0.01];
end
if nargin < 5
    r = 6;
end



tol     = 1e-5;
tol_rho = 100;
maxIter = 200;
% rho     = 1.2;
% mu0     = 0.09; 
rho     = 1.5; mu0 = 1;
max_mu  = 1e6;
eps     = 1e-4;

[height,width,nband] = size(Nhsi);
sizeU           = [height,width,r];
D = reshape(Nhsi, [height*width,nband]);

[~,norm_two,~] = svdsecon(D,1);
mu = mu0/norm_two;
normD   = norm(D,'fro');

%% FFT setting
Eny_x   = ( abs(psf2otf([+1; -1], [height,width,r])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [height,width,r])) ).^2  ;
determ  =  Eny_x + Eny_y;

%% Initializing optimization variables
[u,s,v]= svd(D,'econ');
U              = u(:,1:r)*s(1:r,1:r);
V              = v(:,1:r);

S              = zeros(height*width,nband);%Sparse
E              = zeros(height*width,nband);%Gaussian

M1 = zeros(height*width*r,1);  % multiplier for Dx_U-F1
M2 = zeros(height*width*r,1);  % multiplier for Dy_U-F2
M3 = zeros([height*width,nband]);  % multiplier for D-UV^T-E
W3 = 1;

G1 = diff_x(Pan, size(Pan)); G1=reshape(G1,size(Pan));
G2 = diff_y(Pan, size(Pan)); G2=reshape(G2,size(Pan));

W1 = repmat(G1, [1,1,r]);  W1 = abs(W1);
W2 = repmat(G2, [1,1,r]);  W2 = abs(W2);
W1 = (1-abs(W1)).^q;
W2 = (1-abs(W2)).^q;


rho_XY1 = ones(height,width,r);
rho_XY2 = ones(height,width,r);
update_rho = 0;

MPSNR = zeros(maxIter,1);
MSSIM = zeros(maxIter,1);


%% main loop

disp('Stage 1')
iter = 0;
while iter<maxIter
    iter          = iter + 1;
    %% -Update W1,W2
    W1 = W1.*rho_XY1;
    W2 = W2.*rho_XY2;
    %% -Update F1,F2
    F1 = diff_x(U,sizeU)+M1/mu;
    F2 = diff_y(U,sizeU)+M2/mu;
    F1 = softthre_s(F1,tau(1)/mu,W1(:));
    F2 = softthre_s(F2,tau(2)/mu,W2(:));
    %% -Update U
    diffT_p  = diff_xT(F1-M1/mu,sizeU)+diff_yT(F2-M2/mu,sizeU);
    temp     = (D-E-S+M3/mu)*V;
    numer1   = reshape( diffT_p + temp(:), sizeU);
    U_tensor        = real( ifftn( fftn(numer1) ./ (determ + 1+eps) ) );

    if update_rho
        rho_XY1 = zeros(height,width,r);
        rho_XY2 = zeros(height,width,r);
        UG1 = diff_x(U_tensor, size(U_tensor)); UG1=reshape(UG1,size(U_tensor));
        UG2 = diff_y(U_tensor, size(U_tensor)); UG2=reshape(UG2,size(U_tensor));

        for i=1:r
            rho_XY1(:,:,r) = abs(imcorrfilter(UG1(:,:,r),G1,2,2));
            rho_XY2(:,:,r) = abs(imcorrfilter(UG2(:,:,r),G2,2,2));
        end
    end

    U = reshape(U_tensor,[height*width,r]);



    %% -Update V
    [u,~,v]     = svdecon((D-E-S+M3/mu)'*U);
    V           = u*v';
    %% -Update E
    E = mu*(D-U*V'-S+M3/mu)/(2*beta+mu);

    %% -Update S
    if lambda >100
        % 忽略掉高斯噪音影响
        S  = 0;
    else
        S = softthre_s(D-U*V'-E+M3/mu,lambda/mu,W3);
    end

%     output_image = reshape(U*V',[height,width,nband]);
%     [MPSNR(iter), MSSIM(iter),ERGAS(iter),SAM(iter)] = pwrctv_msqia(Ohsi,output_image);


    %% -Update Multiplier
    leq1 = diff_x(U,sizeU)-F1;
    leq2 = diff_y(U,sizeU)-F2;
    leq3 = D-U*V'-E-S;
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = norm(leq2,'fro')/normD;
    stopC3 = norm(leq3,'fro')/normD;
    if mod(iter,10)==0
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
            ',U_x rele = ' num2str(stopC1,'%2.3e') ',U_y rele = ' num2str(stopC2,'%2.3e')...
            ',X-UV = ' num2str(stopC3,'%2.3e')]);
    end
    if stopC1<tol && stopC2<tol && stopC3 <tol && update_rho
        disp(['Converged at iteration ', num2str(iter)])
        break;
    elseif stopC1<tol_rho*tol && stopC2<tol_rho*tol && stopC3 <tol_rho*tol && ~update_rho
        update_rho = 1;
        disp(['Stage 2 at iteration ', num2str(iter)])
    else
        M1 = M1 + mu*leq1;
        M2 = M2 + mu*leq2;
        M3 = M3 + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
end
output_image = reshape(U*V',[height,width,nband]);
end
function out=softthre_s(a,tau,w)
out = sign(a).* max( abs(a) - tau*w, 0);
end