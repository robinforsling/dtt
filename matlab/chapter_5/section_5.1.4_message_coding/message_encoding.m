function msg = message_encoding(yPsi,RPsi,Psi)
% --- message_encoding() --------------------------------------------------
% Section 5.1.4 Meessage Coding - Encoding Algorithm
%
% 2023-10-30 Robin Forsling

[m,p] = size(Psi);
excl_comp = zeros(1,m*(m-1)/2);

if m == 1
    
    Psi = Psi/norm(Psi);
    Y = Psi*RPsi;
    
elseif m > 1
    
    orth_threshold = 1e-6;
    det_threshold = 1e-5; 
    
    % Check if Phi has orthonormal rows:
    C = Psi*Psi';
    if sum(sum(abs(C-eye(m))))/m > orth_threshold 
        warning('Phi does not have orthonormal rows'); 
        for i = 1:m; Psi(i,:) = Psi(i,:)/norm(Psi(i,:)); end
    end

    U = Psi;
    Y = U(1,:)*RPsi(1,1);
    excl_comp = [];
    for i = 1:m-1
        idx_vec = get_idx_perms(p,i);
        for j = 1:size(idx_vec,1)
            idx = idx_vec(j,:);
            A = [];
            for k = 1:i; A = [A ; U(k,idx)*RPsi(k,k)]; end
            if abs(det(A)) >= det_threshold
                v = U(i+1,:)*RPsi(i+1,i+1); 
                v(idx) = [];
                Y = [Y v];
                excl_comp = [excl_comp idx];
                break;
            end
        end
    end 
end

msg.yPsi = yPsi;
msg.PhiK = Y;
msg.J = excl_comp;

end


