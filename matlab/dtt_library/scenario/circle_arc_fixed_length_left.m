function par = circle_arc_fixed_length_left(varargin)
% --- circle_arc_fixed_length_left() --------------------------------------
% Creates 2D circle arc trajectory of length L and with radius r. The turn
% is to the left. x0 defines the first point and head0 defines the initial
% heading.
%
% Possible function calls:
% 1. circle_arc_fixed_length_left(in_par,r,L,N) or
% 2. circle_arc_fixed_length_left(x0,head0,r,L,N)
% where N is optional
%
% 2023-10-30 Robin Forsling

par = motion_primitive_parameters;
N = 2;

% --- HANDLE INPUTS -------------------------------------------------------
if nargin > 0 
    if are_fieldnames_same(varargin{1},par)
        in_par = varargin{1};
        if nargin >= 3 
            x0 = in_par.X(:,end); head0 = in_par.head(end); r = varargin{2}; L = varargin{3};  
            if nargin >= 4; N = varargin{4}; end
        else; error('illegal input...')
        end
    else
        if nargin >= 4
            x0 = varargin{1}; head0 = varargin{2}; r = varargin{3}; L = varargin{4}; 
            if nargin >= 5; N = varargin{5}; end
        else; error('illegal input...')
        end
    end
end


X = zeros(2,N); 
turn_ang = L/r;
ang_vec = turn_ang*linspace(0,N,N)/N;
Ld = 0;

for k = 2:N
    ang = ang_vec(k);
    X(:,k) = r*[sin(ang) ; (1-cos(ang))];
    Ld = Ld+norm(X(:,k)-X(:,k-1));
end

par.N = N;
par.X = rotate_motion(X,head0) + x0(1:2);
par.head = ang_vec + head0;
par.L = L;
par.Ld = Ld;

