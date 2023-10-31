function varargout = plot_confidence_ellipse(x0,S,gamma,varargin)
% --- plot_confidence_ellipse() -------------------------------------------
% Plot confidence ellipse around x0 with gamma derived from chi2
% distribution.
%
% 2023-10-30 Robin Forsling

L = chol(S,'lower');
N = 1000;
a = linspace(0,2*pi,N);
x = L*sqrt(gamma)*[cos(a) ; sin(a)];
h = plot(x0(1)+x(1,:),x0(2)+x(2,:),varargin{:});
if nargout > 0; varargout{1} = h; end