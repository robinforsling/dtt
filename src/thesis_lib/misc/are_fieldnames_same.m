function same = are_fieldnames_same(a,b)
% --- are_fieldnames_same() -----------------------------------------------
% Check if the structs a and b have the same number of fields and that the
% fields are the same.
%
% 2023-10-30 Robin Forsling

same = 0;

if ~isstruct(a) || ~isstruct(b); return; end

fa = fieldnames(orderfields(a));
fb = fieldnames(orderfields(b));

if length(fa) == length(fb)
    same = 1;
    for k = 1:length(fa)
        if ~strcmp(fa{k},fb{k}); same = 0; end
    end
end

