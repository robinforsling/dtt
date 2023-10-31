function C = generate_all_mat_combs(A)
% --- all_mat_combs(A) ----------------------------------------------------
% Input is a multidimensional cell array A of size n x m. Each element in A
% is a vector. The function returns a N x 1 cell array containing all N
% combinations of of n x m matrices given the corresponding elements of A.
%
% 2023-10-30 Robin Forsling

[n,m] = size(A); N = n*m;
nc_vec = get_number_of_comb_vec(A);
nc_tot = prod(nc_vec);

C = cell(nc_tot,1);
for icomb = 1:nc_tot; M = zeros(n,m); C{icomb} = M; end

nc_for = 1;
nc_aft = nc_tot;
for k = 1:N
    icomb = 1;
    nci = nc_vec(k);
    nc_aft = nc_aft / nci;
    for f = 1:nc_for
        for i = 1:nci
            for a = 1:nc_aft 
                C{icomb}(k) = A{k}(i);
                icomb = icomb + 1;
            end
        end
    end
    nc_for = nc_for*nci;
end

end


% Get number of combinations as vector
function nc_vec = get_number_of_comb_vec(A)
    
    [n,m] = size(A); na = n*m;
    nc_vec = zeros(1, na);
    for k = 1:na; nc_vec(k) = length(A{k}); end
end