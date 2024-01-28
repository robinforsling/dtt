function par = circle_arc_fixed_length_right(varargin)
%  --- circle_arc_fixed_length_right() ------------------------------------
% Creates 2D circle arc trajectory of length L and with radius r. The turn
% is to the right. pos0 defines the first point and head0 defines the 
% initial heading.
%
% Possible function calls:
% 1. circle_arc_fixed_length_right(in_par,r,L,N,Ts) or
% 2. circle_arc_fixed_length_right(pos0,head0,r,L,N,Ts)
% where N and Ts are optional
%
% 2023-10-30 Robin Forsling

par = motion_primitive_parameters;
N = 2;
nx = 6;
Ts = 1;

% --- HANDLE INPUTS -------------------------------------------------------
if nargin > 0 
    if are_fieldnames_same(varargin{1},par)
        in_par = varargin{1};
        if nargin >= 3 
            pos0 = in_par.X(1:2,end); head0 = in_par.head(end); r = varargin{2}; L = varargin{3};
            if nargin >= 4; N = varargin{4}; end
        else; error('illegal input...')
        end
    else
        if nargin >= 4
            pos0 = varargin{1}; head0 = varargin{2}; r = varargin{3}; L = varargin{4};
            if nargin >= 5; N = varargin{5}; end
        else; error('illegal input...')
        end
    end
end


X = zeros(2,N); 
turn_ang = L/r;
ang_vec = turn_ang*linspace(0,N,N)/N;
Ld = 0;

v = L/(Ts*(N-1));
a = v^2/r;

for k = 2:N
    ang = ang_vec(k);
    X(1:2,k) = r*[sin(ang) ; (cos(ang)-1)];
    X(3:4,k) = v*[cos(ang) ; -sin(ang)];
    X(5:6,k) = a*[-sin(ang) ; -cos(ang)];
    Ld = Ld + norm(X(1:2,k)-X(1:2,k-1));
end
ang_vec = -ang_vec;

par.N = N;
par.X = rotate_motion(X,head0) + [pos0 ; zeros(nx-2,1)];
par.head = ang_vec + head0;
par.L = L;
par.Ld = Ld;

