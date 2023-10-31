function s = get_datetime(varargin)
% --- get_datetime() ------------------------------------------------------
% Returns current date and time.
%
% 2023-10-30 Robin Forsling

if nargin > 0; fmt = varargin{1}; else; fmt = 1; end

switch fmt
    case 1; s = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    otherwise; error('get_datetime: unknown format')
end