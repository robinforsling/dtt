function remove_ticks_and_ticklabels(varargin)
% --- remove_ticks_and_ticklabels() ---------------------------------------
% Remove all ticks and tick labels.
%
% 2023-10-30 Robin Forsling

if nargin > 0; axis = varargin{1}; else; axis = 'both'; end

% X
if strcmpi(axis,'x') || strcmpi(axis,'both')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
end

% Y
if strcmpi(axis,'y') || strcmpi(axis,'both')
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end