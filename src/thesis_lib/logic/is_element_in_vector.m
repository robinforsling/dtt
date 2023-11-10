function [in_vector,idx] = is_element_in_vector(element,vector)
% --- is_element_in_vector() ----------------------------------------------
% Checks if vector contains specified element.
%
% 2023-10-30 Robin Forsling

idx = [];
in_vector = 0;
n = length(vector);
for k = 1:n
    if element == vector(k)
        in_vector = 1;
        idx = [idx k];
    end
end

