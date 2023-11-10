function set_fontsize_all(varargin)
% --- set_fontsize_all() --------------------------------------------------
% Sets all text fonts to specified size.
%
% 2023-10-30 Robin Forsling

if nargin > 0; font_size = varargin{1}; else; font_size = 12; end

set(findall(gcf,'-property','FontSize'),'FontSize',font_size);