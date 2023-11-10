function par = motion_primitive_parameters
% --- motion_primitive_parameters() ---------------------------------------
% Returns struct of motion primitive parameters.
%
% 2023-10-30 Robin Forsling

par.N = [];                                                                 % Number of points
par.X = [];                                                                 % Trajectory/sequence of points
par.head = [];                                                              % Heading corresponding to each point
par.L = [];                                                                 % Length of continuous curve
par.Ld = [];                                                                % Length of discretized curve