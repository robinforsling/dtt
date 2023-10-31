function my_struct = set_struct_fields_empty(my_struct)
% --- set_struct_fields_empty() -------------------------------------------
% Sets all fields of my_struct to empty fields.
%
% 2023-10-30 Robin Forsling

if isstruct(my_struct)
    fn = fieldnames(my_struct);
    for k = 1:length(fn)
        my_struct.(fn{k}) = [];
    end
else
    error('set_struct_fields_empty: input is not a struct...')
end