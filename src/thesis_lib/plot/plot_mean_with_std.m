function varargout = plot_mean_with_std(x,mu,sigma,varargin)
% --- plot_mean_with_std() ------------------------------------------------
% Plots mean mu and confidence area between [mu-sigma,mu+sigma], evaluated 
% at x, where std is the standard deviation.
%
% 2024-06-12 Robin Forsling

[x,mu,sigma] = reshape_inputs(x,mu,sigma);

% DEFAULTS:
clr = 0.5*[1 1 1];
opacity = 0.25;

% HANDLE INPUTS:
if nargin > 3; clr = varargin{1}; end
if nargin > 4; opacity = varargin{2}; end

% CONFIDENCE AREA:
ylow = mu-sigma;
yhigh = mu+sigma;

xr = [x fliplr(x)];
yr = [ylow fliplr(yhigh)];

% PLOT:
hold on

harea = fill(xr,yr,clr);
harea.FaceAlpha = opacity;
harea.EdgeColor = clr;
harea.EdgeAlpha = opacity;

hmean = plot(x,mu,'-');
hmean.LineWidth = 1.5;
hmean.Color = clr;

if nargout > 0; varargout{1} = hmean; end
if nargout > 1; varargout{2} = harea; end

end


function varargout = reshape_inputs(varargin)
    for i = 1:nargin
        [m,n] = size(varargin{i});
        varargout{i} = reshape(varargin{i},1,m*n);
    end
end