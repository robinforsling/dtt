function str = get_tikz_plot_coordinates(varargin)
% --- get_tikz_plot_coordinates() -----------------------------------------
% Generate tikz code adapted for "\draw plot coordinates" commands.
%
% 2023-10-30 Robin Forsling

if nargin == 1; x = varargin{1}(1,:); y = varargin{1}(2,:); 
elseif nargin == 2; x = varargin{1}; y = varargin{2}; 
end
N = length(x);
str = '{';
for k = 1:N
    if is_numeric_integer(x(k)) && is_numeric_integer(y(k))
        str = strcat(str,sprintf('(%d,%d)',round(x(k)),round(y(k))));
    elseif is_numeric_integer(x(k))
        str = strcat(str,sprintf('(%d,%3.3f)',round(x(k)),y(k)));
    elseif is_numeric_integer(y(k))
        str = strcat(str,sprintf('(%3.3f,%d)',x(k),round(y(k))));
    else
        str = strcat(str,sprintf('(%3.3f,%3.3f)',x(k),y(k)));
    end
end
str = strcat(str,'};')                          ;

end

function b = is_numeric_integer(v)
    g = 1e-5;
    if abs(mod(v,1)) <= g; b = true; else; b = false; end
end