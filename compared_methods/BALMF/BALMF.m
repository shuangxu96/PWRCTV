function [X_new,U,V,model] = BALMF(Y, Rank, varargin)
% This code carries out the BALMF (band-wise asymmetric Laplacian noise
% based matrix factorization) method proposed by [1]. 
% 
% Input
%   Y: Noisy HSI in tensor format.
%   Rank: Rank.
%   varargin: Other arguments, including
%        'MaxIter': the maximum number of the outer iteration. (Default: 20)
%        'InMaxIter': the maximum number of the inner iteration. (Default: 30)
%        'Tol': the convergence tolerance of the outer iteration. (Default: 1e-4)
%        'InTol': the convergence tolerance of the inner iteration. (Default: 1e-5)
%        'gamma': scaling parameter. (See details in line 10 of Algorithm1 of [1], Default: 1.5)
%        'rho1': penalty parameter. (Default: 10)
%        'GT': ground truth. If provided, the code will print PSNR in each step; otherwise not (Default: [])
%
% Output
%   X_new: Denoised HSI.
%   U: Factor matrix.
%   V: Factor matrix.
%   model: AL model, a structure with following fields:
%       model.kappa: asymmetric parameter in AL distribution.
%       model.lambda: scale parameter in AL distribution.
% 
% Copyright (c) 2023 Shuang Xu
% Email: xu.s@outlook.com; xs@nwpu.edu.cn
%
% If the code is used in your scientific research, please cite the paper.
% [1] Shuang Xu, Xiangyong Cao, Jiangjun Peng, Qiao Ke, Cong Ma and Deyu
% Meng. Hyperspectral Image Denoising by Asymmetric Noise Modeling. IEEE
% TGRS, 2023.


% =========================
% parse the input arguments
% =========================
[MaxIter, InMaxIter, Tol, InTol, gamma, rho1, GT] = ...
    process_options(varargin, ...
    'MaxIter',     20, ...
    'InMaxIter',   30, ...
    'Tol',         1e-4, ...
    'InTol',       1e-5, ...
    'gamma',       1.5, ...
    'rho1',        10,...
    'GT',          []);
nway = ndims(Y);
if nway==3
    [h,w,B] = size(Y);
    Y = reshape(Y, [], size(Y,3)); % convert into matrix format
end

N = h*w;
if ~isempty(GT) % make sure GT is in tensor format
    if ismatrix(GT)
        GT = reshape(GT, h,w,B);
    end
end

% =========================
% initialization
% =========================
start_time = tic;
X = Init_BALMF(Y, Rank);
[V, S, U] = svdsecon(X', Rank);
init_time = toc(start_time);
U = bsxfun(@times, U, sqrt(diag(S))');
V = bsxfun(@times, V, sqrt(diag(S))');

X_old = U*V';
kappa = 0.5*ones(1,B);

fprintf('\n Model has been initialized (%6.3fs Time Elapsed) \n', init_time)

% =========================
% main loop
% =========================
RelErrorList = zeros(1,MaxIter);
M1 = []; 
disp('-----------------------------------------------------------------------------------------------------');
disp(['    OutIter    ','InIter       ','InRelErr    ', 'InTimeElapsed      ', 'OutRelErr        ', 'OutTimeElapsed  ']);
disp('-----------------------------------------------------------------------------------------------------');


for iter = 1:MaxIter
    start_time = tic;
    Noise = Y-X_old;

    % Update lambda
    rho = bsxfun(@times, Noise>0, kappa)+bsxfun(@times, Noise<=0, 1-kappa);
    lambda = N./sum( abs(Noise).*rho  ,1);

    % Update kappa
    eta = sum(Noise,1).*lambda;
    kappa = (2*N+eta-sqrt(4*N^2+eta.^2))./(2*eta);

    % Update U,V
    W = bsxfun(@times, rho, lambda);
    if sum(isnan(W(:)))==1
        W_max = max(max(W(~isnan(W))));
        W(isnan(W)) = W_max;
    end
    [U, V] = WL1LRMF_HQS(Y, W, Rank, InMaxIter, U, V, InTol, iter, gamma, rho1);
    X_new = U*V';

    % convergence?
    con_time = toc(start_time);
    RelError = norm(X_new(:)-X_old(:))/norm(X_old(:));
    RelErrorList(iter) = RelError;

    disp([sprintf('     %3d     ',iter), '    END   ',...
        '    n/a' , '             n/a', ...
        sprintf('              %11.7f',RelError), sprintf('       %6.3f',con_time)  ]);
    if ~isempty(GT)
        X = reshape(X_new, h,w,B);
        M1 = [M1,PSNR3D(255*X,255*GT)]; 
        fprintf('   ------> MPSNR = %6.3f <----- \n', m1);
    end
    disp('-----------------------------------------------------------------------------------------------------');
    disp(['    OutIter    ','InIter       ','InRelErr    ', 'InTimeElapsed      ', 'OutRelErr        ', 'OutTimeElapsed  ']);
    disp('-----------------------------------------------------------------------------------------------------');

    if RelError<Tol
        break
    end

    X_old = X_new;
end

model.lambda = lambda;
model.kappa = kappa;
if ~isempty(GT)
    model.MPSNR = M1(1:iter);
end
model.RelError = RelErrorList(1:iter);
if nway==3
    X_new = reshape(X_new, h,w,B);
end

end % end BALMF


function [U, V] = WL1LRMF_HQS(Y, W, k, max_iter, U, V, tol, out_iter, gamma, rho1)

% rho1 = 10;
rho2 = rho1*mean(W(:))/2;
X_old = U*V';

for iter=1:max_iter
    start_time = tic;
    Z = Shrink(W.*(Y-X_old), 1/rho1);
    numer = W.^2.*Y - W.*Z + rho2*U*V'/rho1;
    denom = W.^2+rho2/rho1;
    X_new = numer./denom; % X_new = clamp(X_new, 0, 1);
    [V, S, U] = svdsecon(X_new', k);
    U = bsxfun(@times, U, sqrt(diag(S))');
    V = bsxfun(@times, V, sqrt(diag(S))');
    con_time = toc(start_time);
    RelError = norm(X_new(:)-X_old(:))/norm(X_old(:));

    if iter == 1 || mod(iter, 10) == 0
        disp([sprintf('     %3d     ',out_iter), sprintf('   %3d    ',iter),...
            sprintf('  %11.7f',RelError), sprintf('      %6.3f',con_time), ...
            '              n/a' , '              n/a'  ]);
    end

    if RelError<tol
        break
    end
    X_old = X_new;
    rho1 = gamma*rho1;
    rho2 = gamma*rho2;
end

end % end WL1LRMF_HQS


function [U,S,V] = svdsecon(X,k)
% Input:
% X : m x n matrix
% k : extracts the first k singular values
%
% Output:
% X = U*S*V' approximately (up to k)
%
% Description:
% Does equivalent to svds(X,k) but faster
% Requires that k < min(m,n) where [m,n] = size(X)
% This function is useful if k is much smaller than m and n
% or if X is sparse (see doc eigs)
%
% Vipin Vijayan (2014)

%X = bsxfun(@minus,X,mean(X,2));
[m,n] = size(X);
assert(k <= m && k <= n, 'k needs to be smaller than size(X,1) and size(X,2)');

if  m <= n
    C = X*X';
    [U,D] = eigs(C,k);
    clear C;
    if nargout > 2
        V = X'*U;
        s = sqrt(abs(diag(D)));
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
    C = X'*X;
    [V,D] = eigs(C,k);
    clear C;
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(abs(diag(D)));
    U = bsxfun(@(x,c)x./c, U, s');
    S = diag(s);
end
end % end - svdsecon

function [varargout] = process_options(args, varargin)

args = prepareArgs(args); % added to support structured arguments
% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
    error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
    error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
    warn = 1;
    nout = n / 2;
else
    warn = 0;
    nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
    varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
for i=1:2:length(args)
    found = 0;
    for j=1:2:n
        if strcmpi(args{i}, varargin{j}) || strcmpi(args{i}(2:end),varargin{j})
            varargout{(j + 1)/2} = args{i + 1};
            found = 1;
            break;
        end
    end
    if (~found)
        if (warn)
            warning(sprintf('Option ''%s'' not used.', args{i}));
            args{i}
        else
            nunused = nunused + 1;
            unused{2 * nunused - 1} = args{i};
            unused{2 * nunused} = args{i + 1};
        end
    end
end

% Assign the unused arguments
if (~warn)
    if (nunused)
        varargout{nout} = unused;
    else
        varargout{nout} = cell(0);
    end
end

end % end process_options

function out = prepareArgs(args)
% Convert a struct into a name/value cell array for use by process_options
%
% Prepare varargin args for process_options by converting a struct in args{1}
% into a name/value pair cell array. If args{1} is not a struct, args
% is left unchanged.
% Example:
% opts.maxIter = 100;
% opts.verbose = true;
% foo(opts)
%
% This is equivalent to calling
% foo('maxiter', 100, 'verbose', true)

% This file is from pmtk3.googlecode.com


if isstruct(args)
    out = interweave(fieldnames(args), struct2cell(args));
elseif ~isempty(args) && isstruct(args{1})
    out = interweave(fieldnames(args{1}), struct2cell(args{1}));
else
    out = args;
end

end % end prepareArgs


function [L, S] = Init_BALMF(X, Rank)

[M, N] = size(X);

lambda = 1 / sqrt(max(M,N));
mu = 10*lambda;
max_iter = 7;
% initial solution
L = zeros(M, N);
S = zeros(M, N);
Y = zeros(M, N);

for iter = (1:max_iter)
    L = MatShrink(X - S + (1/mu)*Y, 1/mu, Rank);
    S = Shrink(X - L + (1/mu)*Y, lambda/mu);
    Z = X - L - S;
    Y = Y + mu*Z;
end
end % end Init_BALMF

function y = Shrink(x, gamma)
y=sign(x).*max(abs(x)-gamma,0);
end % end Shrink

function y = MatShrink(x, gamma, Rank)
[u, s, v] = svdsecon(x, Rank);
y = u*Shrink(s, gamma)*v';
end % end MatShrink

function y = clamp(x,a,b)

% clamp - clamp a value
%
%   y = clamp(x,a,b);
%
% Default is [a,b]=[0,1].
%
%   Copyright (c) 2004 Gabriel Peyr�

if nargin<2
    a = 0;
end
if nargin<3
    b = 1;
end

if iscell(x)
    for i=1:length(x)
        y{i} = clamp(x{i},a,b);
    end
    return;
end

y = max(x,a);
y = min(y,b);
end % end clamp