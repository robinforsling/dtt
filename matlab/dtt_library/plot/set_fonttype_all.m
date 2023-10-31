function set_fonttype_all(font_str)
% --- set_fonttype_all() --------------------------------------------------
% Changes all fonts according to font_str
%
% 2023-10-30 Robin Forsling

if isempty(font_str) || isnumeric(font_str)
    warning('set_fonttype_all: input must be characters of valid font type')
    disp('valid fonts:')
    listfonts
    return
end

fh = findall(0,'Type','Figure');
txt_obj = findall(fh,'Type','text');
set(txt_obj,'FontName',font_str);