function [xhat,P,varargout] = extract_estimate(est)
% --- extract_estimate() --------------------------------------------------
% Extract xhat, P and H (optional) from estimate struct est.
%
% 2023-10-30 Robin Forsling

xhat = est.xhat;
P = est.P;

if nargout > 2
    if isfield(est,'H'); varargout{1} = est.H;
    else; varargout{1} = eye(length(xhat));
    end
end