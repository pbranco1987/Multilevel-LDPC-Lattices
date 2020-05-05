function D = get_distances(H,j)
[m,n] = size(H);

Jold = false(1,n);
Iold = false(m,1);
Jnew = false(1,n);
Jnew(j) = 1;
Inew = logical(H(:,j));
D = ones(m,1)*Inf;
D(Inew) = 1;

for d = 2:m
    Jold = Jold | Jnew;
    Jnew = Inew'*H > 0;
    Jnew = (Jnew-Jold) > 0;
    Iold = Iold | Inew;
    Inew = H*Jnew' > 0;
    Inew = (Inew-Iold) > 0;
    D(Inew) = d;
    if ~any(Inew)
        break
    end
end
