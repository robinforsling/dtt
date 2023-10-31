function [xhat,P,varargout] = le_fusion(varargin)
% --- le_fusion() ---------------------------------------------------------
% Fusion of y1 and y2 using the largest ellipsoid method
%
% Function call:
%   Standard LE:
%       1. le_fusion(y1,R1,y2,R2)
%       2. le_fusion(y1,R1,y2,R2,H)
%   Precomputed gains:
%       3. le_fusion(y1,R1,y2,R2,K1,K2,R12)
%       4. le_fusion(y1,R1,y2,R2,H,K1,K2,R12)
%
% 2023-10-30 Robin Forsling

y1 = varargin{1}; R1 = varargin{2}; y2 = varargin{3}; R2 = varargin{4};
n = length(y1);

precomp_gains = 0;

switch nargin 
    case 4; H2 = eye(n);
    case 5; H2 = varargin{5}; 
    case 7; H2 = eye(n); K1 = varargin{5}; K2 = varargin{6}; R12 = varargin{7}; precomp_gains = 1;
    case 8; H2 = varargin{5}; K1 = varargin{6}; K2 = varargin{7}; R12 = varargin{8}; precomp_gains = 1;
    otherwise; error('invalid inputs')
end

R1 = real(make_symmetric(R1));
R2 = real(make_symmetric(R2));

if ~precomp_gains % Standard version

    I1 = inv(R1); I2 = make_symmetric(H2'/R2*H2);
    [U1,D1] = eig(make_symmetric(I1));
    T1 = inv(sqrtm(D1))*U1.';
    [U2,D2] = eig(make_symmetric(T1*I2*T1'));
    T = U2.'*T1;
    Ti = inv(T);
    
    c_regularize = 0.99;
    Gammai = Ti*diag(min([ones(1,n) ; diag(D2)']))*Ti';
    R12 = R1*Gammai*H2'*R2; R12 = c_regularize*real(R12);
    
    y = [y1;y2]; R = [R1 R12 ; R12' R2]; H = [eye(n);H2];
    P = make_symmetric(inv(H'/R*H));
    K = P*H'/R; K1 = K(:,1:n); K2 = K(:,n+1:end);
    xhat = K*y;

else % Precomputed gains

    A = H2*R1-R12';
    B = H2*R1*H2'+R2-H2*R12-R12'*H2';
    xhat = K1*y1 + K2*y2;
    P = make_symmetric(R1 - A'/B*A);
end

if nargout == 5
    varargout{1} = K1;
    varargout{2} = K2;
    varargout{3} = R12;
end

