function a = d2r(varargin)
% --- d2r() ---------------------------------------------------------------
% Converts from degrees to radians.
%
% 2023-10-30 Robin Forsling

a = pi/180;
if nargin >= 1
    a = varargin{1}*a;
end
