function h = plot_confidence_interval(x,ylow,yhigh,varargin)
% --- plot_confidence_interval() ------------------------------------------
% Plots region/area between lower and upper confidence bounds.
%
% 2024-01-18 Robin Forsling

% DEFAULTS:
clr = 0.5*[1 1 1];
opacity = 0.25;

% HANDLE INPUTS:
if nargin > 3; clr = varargin{1}; end
if nargin > 4; opacity = varargin{2}; end

% CONFIDENCE AREA:
xr = [x fliplr(x)];
yr = [ylow fliplr(yhigh)];

% PLOT:
h = fill(xr,yr,clr);
h.FaceAlpha = opacity;
h.EdgeColor = clr;
h.EdgeAlpha = opacity;
