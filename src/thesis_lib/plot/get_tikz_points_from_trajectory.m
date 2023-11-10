function str = get_tikz_points_from_trajectory(varargin)
% --- get_tikz_points_from_trajectory() -----------------------------------
% Returns str of tikz curve points from input vectors. 
%
% Output format:
% '(x0,y0) -- (x1,y1) -- ... -- (xN,yN);'
%
% 2023-10-30 Robin Forsling

if nargin > 1; X = varargin{1}; Y = varargin{2}; 
else; X = varargin{1}(1,:); Y = varargin{1}(2,:);
end

N = length(X);
if N ~= length(Y); error('input vectors must be of equal length'); end

% INITIALIZE
x0 = X(1); y0 = Y(1);
if mod(x0*10,1) ~= 0 || mod(y0*10,1) ~= 0
    str = sprintf('(%4.3f,%4.3f)',x0,y0);
elseif mod(x0,1) ~= 0 || mod(y0,1) ~= 0
    str = sprintf('(%4.1f,%4.1f)',x0,y0);
else 
    str = sprintf('(%d,%d)',x0,y0);
end

% LOOP OVER POINTS
for k = 2:N
    x = X(k); y = Y(k);
    if mod(x*10,1) ~= 0 || mod(y*10,1) ~= 0
        str = strcat(str,'--',sprintf('(%4.3f,%4.3f)',x,y));
    elseif mod(x0,1) ~= 0 || mod(y0,1) ~= 0
        str = strcat(str,'--',sprintf('(%4.1f,%4.1f)',x,y));
    else 
        str = strcat(str,'--',sprintf('(%d,%d)',x,y));
    end
end

str = strcat(str,';');
