function xcart = convert_pol_to_cart(x)
% --- convert_pol_to_cart() -----------------------------------------------
% Converts from polar to Cartesian coordinates (2D).
%
% 2023-10-30 Robin Forsling

xcart = [   x(1)*cos(x(2)) ; ...
            x(1)*sin(x(2)) ];