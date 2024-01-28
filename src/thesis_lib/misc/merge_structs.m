function s = merge_structs(varargin)
% --- merge_structs() -----------------------------------------------------
% Merges input structs into one struct.
%
% 2024-01-16 Robin Forsling

if nargin == 0; s = []; return; end

if isstruct(varargin{1}); s = varargin{1};
else; error('input #1 is not a valid struct...'); 
end

for i = 2:nargin
    v = varargin{i};
    if ~isstruct(v); error(sprintf('input #%d is not a valid struct...',i)); end % Check input is a struct
    fv = fieldnames(v); nv = length(fv);
    fs = fieldnames(s); ns = length(fs);
    for j = 1:nv
        for k = 1:ns
            if strcmp(fv(j),fs(k)); error('overlapping struct fields...'); end % Check if fields overlap
        end
        s.(fv{j}) = v.(fv{j});
    end
end