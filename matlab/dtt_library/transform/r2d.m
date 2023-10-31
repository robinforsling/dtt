function a = r2d(varargin)
% --- r2d() ---------------------------------------------------------------
% Converts from radians to degrees.
%
% 2023-10-30 Robin Forsling

a = 180/pi;
if nargin >= 1
    a = varargin{1}*a;
end
