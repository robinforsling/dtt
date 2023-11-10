function varargout = plot_ellipse(x0,S,varargin)
% --- plot_ellipse() ------------------------------------------------------
% Plot ellipse with shape matrix S around x0.
%
% 2023-10-30 Robin Forsling

L = chol(S,'lower');
N = 1000;
a = linspace(0,2*pi,N);
x = L*[cos(a) ; sin(a)];
h = plot(x0(1)+x(1,:),x0(2)+x(2,:),varargin{:});

if nargout > 0; varargout{1} = h; end