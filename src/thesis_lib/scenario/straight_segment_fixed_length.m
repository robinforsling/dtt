function par = straight_segment_fixed_length(varargin)
% --- straight_segment_fixed_length() -------------------------------------
% Creates a 2D straight line segment of length L.
%
% Possible function calls:
% 1. straight_segment_fixed_length(in_par,L,N,Ts) or
% 2. straight_segment_fixed_length(pos0,head0,L,N,Ts)
% where N and Ts are optional
%
% 2024-01-16 Robin Forsling

par = motion_primitive_parameters;
N = 2;
nx = 6;
Ts = 1;

% --- HANDLE INPUTS -------------------------------------------------------
if nargin > 0 
    if are_fieldnames_same(varargin{1},par)
        in_par = varargin{1};
        if nargin >= 2 
            pos0 = in_par.X(1:2,end); head0 = in_par.head(end); L = varargin{2};
            if nargin >= 3; N = varargin{3}; end
            if nargin >= 4; Ts = varargin{4}; end
        else; error('circle_arc_fixed_length: illegal input...')
        end
    else
        if nargin >= 3
            pos0 = varargin{1}; head0 = varargin{2}; L = varargin{3};
            if nargin >= 4; N = varargin{4}; end
            if nargin >= 5; Ts = varargin{5}; end
        else; error('circle_arc_fixed_length: illegal input...')
        end
    end
end

v = L/(Ts*(N-1));
pos = [linspace(0,L,N) ; zeros(1,N)];
vel = [v*ones(1,N) ; zeros(1,N)];
acc = zeros(2,N);    
X = [pos ; vel ; acc]; 

par.N = N;
par.X = rotate_motion(X,head0) + [pos0 ; zeros(nx-2,1)];
par.head = head0*ones(1,N);
par.L = L;
par.Ld = L;

