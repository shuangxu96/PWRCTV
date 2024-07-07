function myshow(X, scale, band)

if nargin==1
    band = NaN;
    scale = 0.01;
end
if nargin==2
    band = NaN;
end

C = size(X,3);

% figure()
if C == 3 || C==1
    imshow(X)
else
    if isnan(band)
        band = [C, uint8(C*0.5), uint8(C*0.1)];
        band(band<1) = 1;
    end
    rsshow(X(:,:,band), scale)
end
