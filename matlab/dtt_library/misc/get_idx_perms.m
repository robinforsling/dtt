function idx_list = get_idx_perms(n,m)
% --- get_idx_perms() -----------------------------------------------------
% Returns a list with all permutations of indices according to n and m.
%
% Robin Forsling 2023-11-02

if m > n; error('m > n'); end

basis_list = get_nk_idx_combs(n,m);
perm_list = flipud(perms(1:m));
idx_list = [];

for k = 1:size(perm_list,1)
    idx_list = [idx_list ; basis_list(:,perm_list(k,:))];
end



