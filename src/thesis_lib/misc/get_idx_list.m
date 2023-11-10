function [idx_list,N] = get_idx_list(n,m)
% --- get_idx_list() ------------------------------------------------------
% Returns list with indices.
%
% 2023-10-30 Robin Forsling

N = get_binom_coeff(n,m);
idx_list = zeros(N,m);
prev_idx = 1:m; 
prev_idx(end) = prev_idx(end)-1;

for i = 1:N
    prev_idx(m) = prev_idx(m)+1;
    for j = m:-1:2
        if prev_idx(j) > n+j-m; pj = prev_idx(j-1) + 1; prev_idx(j-1:m) = pj:pj+(m-j+1); end
    end
    idx_list(i,:) = prev_idx;
end

end