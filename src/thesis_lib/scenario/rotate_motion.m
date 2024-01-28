function Xnew = rotate_motion(X,alpha,varargin)
% --- rotate_motion() -----------------------------------------------------
% Rotates a 2D trajectory X by angle alpha. For rotation around a point xc
% other than the origin, set vargargin{1} = xc.
%
% 2023-10-30 Robin Forsling

if nargin > 2; xc = varargin{1}; else; xc = [0;0]; end

T = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)];
pos = T*(X(1:2,:)-xc)+xc;
vel = T*X(3:4,:);
acc = T*X(5:6,:);
Xnew = [pos;vel;acc];