function x = prox_nn(input, gamma, rank)
if nargin<=2
    rank=inf;
end
[uu,sigma,vv] = svdecon(input);
sigma = diag(sigma);
svp = min(length(find(sigma>gamma)),rank);
if svp>=1
    sigma = sigma(1:svp)-gamma;
else
    svp = 1;
    sigma = 0;
end
x = uu(:,1:svp)*diag(sigma)*vv(:,1:svp)';
end