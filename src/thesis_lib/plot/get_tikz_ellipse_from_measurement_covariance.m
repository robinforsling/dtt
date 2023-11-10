function str = get_tikz_ellipse_from_measurement_covariance(xs,xt,R,length_scale,varargin)
% --- get_tikz_ellipse_from_measurement_covariance() ----------------------
% Returns TikZ code for measurement ellipse.
%
% 2023-10-30 Robin Forsling

if nargin > 4; gamma = varargin{1}; else; gamma = 1; end

xr = (xt-xs)*length_scale;
r = sqrt(xr(1:2)'*xr(1:2));
az = r2d*atan2(xr(2),xr(1));
sr = sqrt(R(1,1))*length_scale;
stheta = sqrt(R(2,2));
g = sqrt(gamma);

str = sprintf('ellipse[x radius = %4.3f, y radius = %4.3f, rotate = %4.3f];',g*sr,g*r*stheta,az);