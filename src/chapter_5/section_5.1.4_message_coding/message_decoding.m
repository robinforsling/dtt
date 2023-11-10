function [yPsi,RPsi,Psi] = message_decoding(msg)
% --- message_decoding() --------------------------------------------------
% Section 5.1.4 Meessage Coding - Decoding Algorithm
%
% 2023-10-30 Robin Forsling

[p,m] = get_dim(msg);
excl_comp = get_excl_comp(msg,m);

Y = msg.PhiK;
R = zeros(m);
U = zeros(m,p);
C = cell(1,m);

for i = 1:m
    C{i} = Y(1:p-i+1);
    Y(1:p-i+1) = [];
end

R(1,1) = norm(C{1});
U(1,:) = C{1}/R(1,1);
for i = 1:m-1
    vidx_excl = excl_comp{i};
    vidx_incl = 1:p; vidx_incl(vidx_excl) = [];
    vincl = C{i+1};
    A = []; b = zeros(i,1);
    for j = 1:i 
        u = C{j}; 
        uincl = u(vidx_incl); uexcl = u(vidx_excl);
        A = [A ; uexcl];
        b(j) = -uincl*vincl';
    end
    vexcl = A\b;
    vfull = vec_conc(vincl,vexcl,vidx_excl);
    C{i+1} = vfull;
    R(i+1,i+1) = norm(vfull);
    U(i+1,:) = vfull/R(i+1,i+1);
end

yPsi = msg.yPsi;
RPsi = R;
Psi = U;

end

function [p,m] = get_dim(msg)
    NJ = length(msg.J);
    NP = length(msg.PhiK);
    m = 1/2 + sqrt(1+8*NJ)/2;
    p = (2*NP+m^2-m)/(2*m);
end

function excl_comp = get_excl_comp(msg,m)
    excl_comp = cell(1,m-1);
    idx_last = 0;
    for i = 1:m-1
        idx = idx_last + (1:i);
        excl_comp{i} = msg.J(idx);
        idx_last = idx(end);
    end
end

function w = vec_conc(v,u,uidx)
    p = length(v)+length(u);
    vidx = 1:p; vidx(uidx) = [];
    nv = length(vidx); nu = length(uidx);
    w = zeros(1,p);
    for i = 1:nv; w(vidx(i)) = v(i); end
    for i = 1:nu; w(uidx(i)) = u(i); end
end

        
