function [typ,ave,stdev] = convergence_rate_stats(v)
% --- convergence_rate_stats() --------------------------------------------
% Compute statistics for the convergence rate analysis.
%
% 2023-10-30 Robin Forsling

% TYPICAL VALUE:
[x,h] = integer_histogram(v);
[~,idx] = max(h);
typ = x(idx);

% AVERAGE:
ave = mean(v);

% STANDARD DEVIATION:
stdev = std(v);