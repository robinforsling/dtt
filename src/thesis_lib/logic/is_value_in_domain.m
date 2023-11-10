function in_domain = is_value_in_domain(val,dom)
% --- is_value_in_domain() ------------------------------------------------
% Check if value is in domain.
%
% 2023-10-30 Robin Forsling

if length(dom) ~= 2; error('domain must be an array of two element corresponding to the boundary of the domain...'); end

dom = sort(dom);
in_domain = 0;
if val >= dom(1) && val <= dom(2)
    in_domain = 1; 
end
