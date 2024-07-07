function norm_tensor = Normalize(orginal_tensor)
[m,n,p] = size(orginal_tensor);
norm_tensor = zeros([m,n,p]);
for band =1:p
    tmp = orginal_tensor(:,:,band);
    tmp = (tmp-min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
    norm_tensor(:,:,band) = tmp;
end

% max_value = max(orginal_tensor(:));
% min_value = min(orginal_tensor(:));
% norm_tensor = (max_value-orginal_tensor)/(max_value-min_value);
end