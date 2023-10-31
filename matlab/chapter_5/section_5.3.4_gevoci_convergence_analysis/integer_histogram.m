function [xhist,yhist] = integer_histogram(data)
% --- integer_histogram() -------------------------------------------------
% Returns frequency of each integer in the matrix data. Elements in M not 
% an integer will be rounded to nearest integer.
%
% 2023-10-30 Robin Forsling

data = round(data(:));
n = length(data);

int_min = min(data);
int_max = max(data);
xhist = int_min:int_max;
yhist = zeros(1,length(xhist));

for k = 1:length(xhist)
    yhist(k) = n - nnz(data-xhist(k));
end