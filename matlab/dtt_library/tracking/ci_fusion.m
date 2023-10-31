function [xhat,P,varargout] = ci_fusion(varargin)
% --- ci_fusion() ---------------------------------------------------------
% Fusion of y1 and y2 using covariance intersection.
%
% Function call:
%   Standard CI:
%       1. ci_fusion(y2,R2,y2,R2)
%       2. ci_fusion(y1,R1,H1,y2,R2,H2)
%   Precomputed weight parameter w:
%       3. ci_fusion(y2,R2,y2,R2,w)
%       4. ci_fusion(y1,R1,H1,y2,R2,H2,w)
%
% 2023-10-30 Robin Forsling

optimize_omega = 1;

if nargin > 3 && nargin < 6
    y1 = varargin{1}; R1 = varargin{2}; y2 = varargin{3}; R2 = varargin{4}; H1 = eye(length(y1)); H2 = eye(length(y2)); 
    if nargin == 5; w = varargin{5}; optimize_omega = 0; end
elseif nargin < 8 
    y1 = varargin{1}; R1 = varargin{2}; H1 = varargin{3}; y2 = varargin{4}; R2 = varargin{5}; H2 = varargin{6}; 
    if nargin == 7; w = varargin{7}; optimize_omega = 0; end
else 
    error('invalid inputs')
end

I1 = inv(R1); I2 = inv(R2); 

if optimize_omega
    f = @(w) trace(inv(w*H1'*I1*H1+(1-w)*H2'*I2*H2) ); 
    w = fminbnd(f,0,1,optimset('Display','off'));
end

P = inv(w*H1'*I1*H1+(1-w)*H2'*I2*H2);
xhat = P*(w*H1'*I1*y1+(1-w)*H2'*I2*y2);

if nargout > 2; varargout{1} = w; end
