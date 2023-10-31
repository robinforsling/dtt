function [h,J] = transform_cart_to_pol()
% --- transform_cart_to_pol() ---------------------------------------------
% Transform from Cartesian (2D) to polar coordinates. Includes Jacobian.
%
% 2023-10-30 Robin Forsling

h = @(x) [  sqrt(x(1)^2+x(2)^2) ; ...
            atan2(x(2),x(1))];

J = @(x) [  x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) ; ...
            -x(2)/(x(1)^2+x(2)^2) x(1)/(x(1)^2+x(2)^2)];
 