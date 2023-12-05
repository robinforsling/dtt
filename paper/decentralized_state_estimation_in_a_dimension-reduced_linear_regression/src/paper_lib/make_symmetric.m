function M = make_symmetric(M)
% --- make_symmetric() ----------------------------------------------------
% Force matrix to be symmetric.
%
% 2023-10-30 Robin Forsling

M = (M+M.')/2;