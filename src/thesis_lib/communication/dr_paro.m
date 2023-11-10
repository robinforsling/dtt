function [H,varargout] = dr_paro(R1,R2,m,varargin)
% --- dr_paro() -----------------------------------------------------------
% Dimension-reduction (DR) based on principle axis restricted optimization 
% (PARO).
%
% VARARGIN
%   fusion_model = varargin{1}: fusion method as string or number
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
            otherwise; error('unknown input string')
        end
    elseif isnumeric(var) && isscalar(var); fusion_method = var;
    else; error('unknown input')
    end
end

% SELECT Psi
switch fusion_method 
    case 1; [H,Jmin] = paro_kf(R1,R2,m); 
    case 2; [H,Jmin] = paro_ci(R1,R2,m);
    otherwise; error('unknown fusion model')
end

if nargout > 1; varargout{1} = Jmin; end

end


% --- KF ------------------------------------------------------------------
function [H,Jmin] = paro_kf(R1,R2,m)
    n = size(R1,1); 
    idx_list = generate_nk_idx_combs(n,m); nidx = length(idx_list);
    J_vec = zeros(nidx,1);
    [U,D] = eig(R2); 
    Di = inv(D); I1 = inv(R1);
    for i = 1:nidx
        idx = idx_list(i,:);
        H = U(:,idx)';
        IH = Di(idx,idx);
        J_vec(i) = loss_kf(I1,IH,H);
    end
    [Jmin,imin] = min(J_vec);
    H = U(:,idx_list(imin,:))';
end

function J = loss_kf(I1,IH,H)
    J = trace(inv(I1+H'*IH*H));
end


% --- CI ------------------------------------------------------------------
function [H,Jmin] = paro_ci(R1,R2,m)
    n = size(R1,1); 
    idx_list = generate_nk_idx_combs(n,m); nidx = length(idx_list);
    J_vec = zeros(nidx,1);
    [U,D] = eig(R2); 
    Di = inv(D); I1 = inv(R1);
    for i = 1:nidx
        idx = idx_list(i,:);
        H = U(:,idx)';
        IH = Di(idx,idx);
        J_vec(i) = loss_ci(I1,IH,H);
    end
    [Jmin,imin] = min(J_vec);
    H = U(:,idx_list(imin,:))';
end

function J = loss_ci(I1,IH,H)
    f = @(w) trace(inv(w*I1+(1-w)*H'*IH*H));
    w = fminbnd(f,0,1,optimset('Display','off'));
    J = trace(inv(w*I1+(1-w)*H'*IH*H));
end

