function example_gevo_extra_parameters
% --- example_gevo_extra_parameters() -------------------------------------
% Example related to extra parameters to be transmitted.
%
% 2023-10-30 Robin Forsling

PM = [4 2 ; 6 2 ; 6 3 ; 9 3 ; 9 5 ; 9 7 ; 12 3 ; 12 6 ; 15 9];
fex = @(m,p) (m-1) / (8*(2*p-m+3));

fprintf('\n--- EXTRA PARAMETERS ---\n')
for i = 1:size(PM,1)
    m = PM(i,2);
    p = PM(i,1);
    ex = 100*fex(m,p);
    fprintf('%d & %d & %1.2f%% \n',m,p,ex)
end
fprintf('\n')