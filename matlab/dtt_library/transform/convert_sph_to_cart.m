function xcart = convert_sph_to_cart(x)
% --- convert_sph_to_cart() -----------------------------------------------
% Converts from spherical to Cartesian coordinates (3D).
%
% 2023-10-30 Robin Forsling

xcart =  [  x(1)*cos(x(2))*cos(x(3)) ; ...
            x(1)*sin(x(2))*cos(x(3)) ; ...
            x(1)*sin(x(3)) ];