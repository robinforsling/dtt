function is_soc = is_string_or_char(my_var)
% --- is_string_or_char() -------------------------------------------------
% Checks whether input is string or char array.
%
% 2023-10-30 Robin Forsling

is_soc = isstring(my_var) || ischar(my_var);