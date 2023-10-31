function varargout = ellipse_shape_rotation(varargin)
% --- ellipse_shape_rotation() --------------------------------------------
% Calculates parametrization of an ellipse by its major axis, minor axis
% and its rotation.
% 
% x radius is the major axis
% y radius is the minor axis
% rotation in the normal sense, around the "z axis" where no (zero) 
% rotation aligns the ellipse major axis with the "x axis"
%
% 2023-10-30 Robin Forsling

r2d = 180/pi;
S = cell(nargin,3);

if nargin == 1 && iscell(varargin{1})
    VA = varargin{1};
else
    VA = varargin;
end

NA = length(VA);

for k = 1:NA
    R = VA{k};
    [U,D] = eig(R);
    if size(R,1) ~= 2 || size(R,2) ~= 2 || ~all(diag(D) >= 0); error('ellipse_shape_rotation: dimensionality error...'); end
    rmajor = sqrt(D(2,2));
    rminor = sqrt(D(1,1));
    angle = atan2(U(2,2),U(1,2))*r2d;
    str = sprintf('%3.2f/%3.2f/%3.2f',rmajor,rminor,angle);
    gen_str = sprintf('(%3.2f,%3.2f,%3.2f)',rmajor,rminor,angle);
    tikz_str = sprintf('[x radius = %3.2f, y radius = %3.2f, rotate = %3.2f]',rmajor,rminor,angle);
    S{k,1} = str; S{k,2} = gen_str; S{k,3} = tikz_str;
    %if NA < 2; display(str); display(tikz_str); end
end

loop_str = ''; % Suitable for tikz foreach loops
for k = 1:NA
    loop_str = strcat(loop_str, S{k,1});
    if k < NA; loop_str = strcat(loop_str,','); end
end


if nargout == 0
    fprintf('\nParameters:')
    for k = 1:NA; fprintf('\n%s',S{k,1}); end  
    fprintf('\n\nGeneric TikZ ellipse parameters:')
    for k = 1:NA; fprintf('\n%s',S{k,2}); end  
    fprintf('\n\nTikZ ellipse options:')
    for k = 1:NA; fprintf('\n%s;',S{k,3}); end
    fprintf('\n\n')
elseif nargout == 1
    varargout{1} = S;
else
    varargout{1} = S;
    varargout{2} = loop_str;
end

