function U = gram_schmidt_process(X)
% --- gram_schmidt_process() ----------------------------------------------
% Creates an orthonormal basis U from input basis X. If X is a tall matrix
% the basis is assumed to be given as column vectors, otherwise the basis
% is assumed to be given as row vectors.
% 
% 2023-10-30 Robin Forsling

[n,m] = size(X);
is_tall = n > m;

if ~is_tall; X = X'; [n,m] = size(X); end % Flip matrix (transpose)

U = zeros(n,m);
U(:,1) = X(:,1)/norm(X(:,1));

for i = 2:m
    x = X(:,i); u = x;
    for j = 1:i-1
        uprev = U(:,j);
        u = u - po(uprev,x);
    end
    U(:,i) = u/norm(u);
end

if ~is_tall; U = U'; end  % Flip back

end


% --- PROJECTION OPERATOR -------------------------------------------------
function vu = po(u,v)
    vu = (u'*v/(u'*u))*u;
end

