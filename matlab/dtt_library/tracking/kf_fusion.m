function [xhat,P] = kf_fusion(varargin)
% --- kf_fusion() ---------------------------------------------------------
% Fusion of y1 and y2 under the naive assumption that y1 and y2 are
% mutually uncorrelated.
%
% Function call:
%   1. kf_fusion(y2,R2,y2,R2)
%   2. kf_fusion(y1,R1,H1,y2,R2,H2)
%
% 2023-10-30 Robin Forsling

if nargin == 4; y1 = varargin{1}; R1 = varargin{2}; y2 = varargin{3}; R2 = varargin{4}; H1 = eye(length(y1)); H2 = eye(length(y2)); 
elseif nargin == 6; y1 = varargin{1}; R1 = varargin{2}; H1 = varargin{3}; y2 = varargin{4}; R2 = varargin{5}; H2 = varargin{6}; 
else; error('invalid inputs')
end

I1 = inv(R1); I2 = inv(R2); 
P = inv(H1'*I1*H1+H2'*I2*H2);
xhat = P*(H1'*I1*y1+H2'*I2*y2);