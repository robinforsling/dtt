function [Psi,varargout] = dr_gevo(R1,R2,m,varargin)
% --- dr_gevo() -----------------------------------------------------------
% Dimension-reduction (DR) based on generalized eigenvalue optimization 
% (GEVO).
%
% VARARGIN
%   fusion_model = varargin{1}: fusion method as string or number
%   R12 = varargin{2}; xcov, only valid for some fusion methods
%
% 2023-10-30 Robin Forsling

fusion_method = 1;

% HANDLE INPUTS
if nargin > 3
    var = varargin{1};
    if is_string_or_char(var)
        switch lower(var)
            case {'kf','naive','sf-formula'}; fusion_method = 1;
            case 'ci'; fusion_method = 2;
            case {'bsc','known-xcorr'} 
                fusion_method = 3;
                if nargin < 5; error('too few inputs'); else; R12 = varargin{2}; end
            case 'le'; fusion_method = 4;         
            otherwise; error('unknown input string')
        end
    elseif isnumeric(var) && isscalar(var); fusion_method = var;
    else; error('unknown input')
    end
end

R1 = make_symmetric(R1);
R2 = make_symmetric(R2);


% COMPUTE Psi
switch fusion_method 
    case 1; [Psi,Jmin] = gevo_kf(R1,R2,m); 
    case 2; [Psi,Jmin] = gevo_ci(R1,R2,m);
    case 3; [Psi,Jmin] = gevo_bsc(R1,R2,R12,m);
    case 4; [Psi,Jmin] = gevo_le(R1,R2,m);
    otherwise; error('unknown fusion model')
end

if nargout > 1; varargout{1} = Jmin; end

end


% --- KF ------------------------------------------------------------------
function [Psi,Jmin] = gevo_kf(R1,R2,m)
    Q = R1*R1;
    S = R1+R2; 
    [X,D] = eig(Q,S,'qz');
    idx_vec = get_max_idx_vec(D,m);
    Psi = gram_schmidt_process(X(:,idx_vec))';
    Jmin = kf_loss_function(R1,R2,Psi);
end


% --- CI ------------------------------------------------------------------
function [Psi,Jmin] = gevo_ci(R1,R2,m)
    max_iter = 50; J_thres = 0.0001; wmin = 0.001; wmax = 1;                   % Algorithm settings
    I1 = inv(R1);   

    % INITIALIZE
    w = 0.5;  
    J = Inf;

    % ITERATE
    for k = 1:max_iter
        J_prev = J;
        Q = R1*R1/(w*w);
        S = R1/w + R2/(1-w);

        % SOLVE FOR Psi
        [X,D] = eig(Q,S,'qz');
        idx_vec = get_max_idx_vec(D,m);
        Psi = X(:,idx_vec)';

        % SOLVE FOR w
        f = @(w) trace(inv(w*I1 + (1-w)*Psi'/(Psi*R2*Psi')*Psi));
        w = fminbnd(f,wmin,wmax);
        J = ci_loss_function(R1,R2,Psi,w);
    
        if (J_prev-J)/J < J_thres; break; end
    end

    Psi = gram_schmidt_process(Psi);

    if k == max_iter; warning('gevo_ci: reached maximum number of iterations'); end
    Jmin = J;
end


% --- BSC -----------------------------------------------------------------
function [Psi,Jmin] = gevo_bsc(R1,R2,R12,m)
    Q = make_symmetric((R1-R12)'*(R1-R12)); 
    S = make_symmetric(R1+R2-R12-R12');
    [X,D] = eig(Q,S,'qz');
    idx_vec = get_max_idx_vec(D,m);
    Psi = gram_schmidt_process(X(:,idx_vec))';
    Jmin = bsc_loss_function(R1,R2,R12,Psi);
end


% --- LE ------------------------------------------------------------------
function [Psi,Jmin] = gevo_le(R1,R2,m)   
    n = size(R1,1); 
    if m == 1; c = 1; else; c = 0.999; end

%     [U1,D1] = eig(make_symmetric(inv(R1))); 
%     T1 = inv(sqrtm(D1))*U1';
%     [U2,D2] = eig(make_symmetric(T1/R2*T1'));
%     T = U2'*T1; Ti = inv(T);
%     Gammai = Ti*diag(min([ones(1,n) ; diag(D2)']))*Ti';
%     R12 = R1*Gammai*R2; R12 = c*R12;
%     [Psi,Jmin] = gevo_bsc(R1,R2,R12,m); 

    [U1,D1] = eig(make_symmetric(R1)); 
    T1 = inv(sqrtm(D1))*U1';
    [U2,D2] = eig(make_symmetric(T1*R2*T1'));
    T = U2'*T1; Ti = inv(T);
    R12 = c*Ti*diag(min([ones(1,n) ; diag(D2)']))*Ti';
    [Psi,Jmin] = gevo_bsc(R1,R2,R12,m); 
end

