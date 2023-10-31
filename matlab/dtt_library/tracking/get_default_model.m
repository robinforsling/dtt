function [F,Q,h,J,R] =  get_default_model(varargin)
% --- get_default_model() -------------------------------------------------
% Get default sensor and process models.
%
% 2023-10-30 Robin Forsling

if nargin > 0; n = varargin{1}; else; n = 4; end
if nargin > 1; Ts = varargin{2}; else; Ts = 1; end

d2r = pi/180;

F = @(T) [];

m = n/2;

F = @(T) [eye(m) T*eye(m) ; zeros(m) eye(m)];
Q = @(T) [T^3/3*eye(2) T^2/2*eye(m) ; T^2/2*eye(m) T*eye(m)];

if m == 2
    [h,J] = transform_cart_to_pol;
    R = diag([50 1*d2r].^2);
elseif m == 3
    [h,J] = transform_cart_to_sph;
    R = diag([50 1*d2r 1*d2r ].^2);
end




