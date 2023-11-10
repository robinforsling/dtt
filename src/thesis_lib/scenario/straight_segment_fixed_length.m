function par = straight_segment_fixed_length(varargin)
% --- straight_segment_fixed_length() -------------------------------------
% Creates a 2D straight line segment of length L.
%
% Possible function calls:
% 1. straight_segment_fixed_length(in_par,L,N) or
% 2. straight_segment_fixed_length(x0,head0,L,N)
% where N is optional
%
% 2023-10-30 Robin Forsling

par = motion_primitive_parameters;
N = 2;

% --- HANDLE INPUTS -------------------------------------------------------
if nargin > 0 
    if are_fieldnames_same(varargin{1},par)
        in_par = varargin{1};
        if nargin >= 2 
            x0 = in_par.X(:,end); head0 = in_par.head(end); L = varargin{2};
            if nargin >= 3; N = varargin{3}; end
        else; error('circle_arc_fixed_length: illegal input...')
        end
    else
        if nargin >= 3
            x0 = varargin{1}; head0 = varargin{2}; L = varargin{3};
            if nargin >= 4; N = varargin{4}; end
        else; error('circle_arc_fixed_length: illegal input...')
        end
    end
end


X = [linspace(0,L,N) ; zeros(1,N)]; 

par.N = N;
par.X = rotate_motion(X,head0) + x0(1:2);
par.head = head0*ones(1,N);
par.L = L;
par.Ld = L;

