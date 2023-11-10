function [in_vector,idx,abs_diff] = is_element_in_vector_approx(element,vector,varargin)
% --- is_element_in_vector_approx() ---------------------------------------
% Checks if element is approximately equal to entry in vector.
%
% 2023-10-30 Robin Forsling

tolerance = 1e-6;

if nargin > 2; tolerance = varargin{1}; end

idx = [];
abs_diff = [];
in_vector = 0;
n = length(vector);
for k = 1:n
    if abs(element - vector(k)) <= tolerance
        in_vector = 1;
        idx = [idx k];
        abs_diff = [abs_diff abs(element - vector(k))];
    end
end