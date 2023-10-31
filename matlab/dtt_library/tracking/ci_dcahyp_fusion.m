function [xhat,P,varargout] = ci_dcahyp_fusion(y1,R1,y2,D2,H2,varargin)
% --- ci_dcahyp_fusion() --------------------------------------------------
% CI fusion of (y1,R1) and (y2,D2) where D2 is an diagonal covariance
% approximation (DCA) of R2. The method is called DCA hyperrectangle
% enclosing.
%
% 2023-10-30 Robin Forsling

[n2,nx] = size(H2);
I1 = inv(R1); I2 = inv(D2);

if nargin > 5
    w = varargin{1};
else
    w = sdpvar(n2+1,1);
    Y = sdpvar(nx);
     
    options = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);   % SDP options  
    
    U = w(n2+1)*I1 + H2'*(diag(w(1:n2)).*I2)*H2;
    F = [sum(w) == 1, 0 <= w(:) <= 1, [Y eye(nx) ; eye(nx) U] >= 0]; 
    J = trace(Y);
    res = optimize(F,J,options);
    w = value(w); 
    
    if sum(w) <= 0.99; w(:) = 1/(n+1); end
end

P = inv(w(n2+1)*I1 + H2'*(diag(w(1:n2)).*I2)*H2);
xhat = P*(w(n2+1)*I1*y1 + H2'*(diag(w(1:n2)).*I2)*y2); 

if nargout > 2; varargout{1} = w; end