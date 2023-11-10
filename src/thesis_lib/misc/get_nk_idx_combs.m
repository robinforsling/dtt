function idx_list = get_nk_idx_combs(n,k)
% --- get_nk_idx_combs() --------------------------------------------------
% Returns a list with all combinations of indices according to n and k.
%
% 2023-10-30 Robin Forsling

if k > n; error('k > n'); end

nc = factorial(n) / (factorial(k)*factorial(n-k));

idx_list = zeros(nc,k);
prev_idx = 1:k;
prev_idx(end) = prev_idx(end) - 1;

for i = 1:nc
    prev_idx(k) = prev_idx(k) + 1;
    for j = k:-1:2
        if prev_idx(j) > n + j - k
            pj = prev_idx(j-1) + 1;
            prev_idx(j-1:k) = pj:pj+(k-j+1);
        end
    end
    idx_list(i, :) = prev_idx;
end