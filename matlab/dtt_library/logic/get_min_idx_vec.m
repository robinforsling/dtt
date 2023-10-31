function idx_vec = get_min_idx_vec(val,m)
% --- get_min_idx_vec() ---------------------------------------------------
% Return indices corresponding to the m smallest elements of the vector or
% matrix val. If val is a matrix only its diagonal entries are considered.
%
% 2023-10-30 Robin Forsling

if ~isvector(val); val = diag(val); end

idx_vec = zeros(1,m);
for k = 1:m
    [~,idx] = min(val);
    idx_vec(k) = idx;
    val(idx) = Inf;
end