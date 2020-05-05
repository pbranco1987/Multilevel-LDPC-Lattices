function [HH,p,A] = PEG_CS_tri(H,mm)
% check splitting function to create a matrix HH out of matrix H using the
% technique detailed in "Multilevel LDPC Lattices with Efficient Encoding and 
% Decoding and a Generalization of Construction D'" by Paulo Branco and 
% Danilo Silva

[m,n] = size(H);
HH = spalloc(mm,n,sum(sum(H)));

[p,g] = get_parents_tri(H,mm);
fprintf('Matriz PEG CS Triangular\n');

for j = 1:n
    if mod(j,100)==0
        fprintf('(%d/%d)\n',j,n)
    end
    I = find(H(:,j))';
    if j <= mm-g
        k = g+j;
        i0 = p(k);
        if ~ismember(i0,I)
            error('Parent vector is incorrect')
        end
        HH(k,j) = 1;
        I = setdiff(I,i0);
    end
    for i = I
        dcnew = sum(HH,2);
        D = get_distances(HH,j);

        S = (1:min(j+g-1,mm));
        S = intersect(S, find(p==i));
        
      
        if ~isempty(S)
            S = intersect(S, find(D == max(D(S))));
            S = intersect(S, find(dcnew == min(dcnew(S))));
            k = S(randi(length(S)));
            HH(k,j) = 1;
        end
    end
    if has_girth4(HH,j)
        disp(''); disp(['Girth 4 at n = ' num2str(j)])
        HH = HH(:,1:j);
        return
    end
    display_progress(j,n);
end

A = double(ones(m,1)*p(:)' == (1:m)'*ones(1,mm));
if ~all(all(A*HH==H))
    error('A*HH ~= H')
end

function [p,g] = get_parents_tri(H,mm)
m = size(H,1);
g = bandwidth(H,'lower');
if ~all(diag(H,-g))
    error('Input matrix must be approximate upper triangular with ones along the lowest nonzero diagonal.')
end
p = [1:m NaN*ones(1,mm-m)];
for k = m+1:mm
    dc = sum(H,2);
    count = histcounts(p,'BinMethod','integers')';
    dcnorm = dc./count;
    j = k-g;
    S = find(H(:,j));
    S = intersect(S, find(dcnorm == max(dcnorm(S))));
    i = S(randi(length(S)));
    p(k) = i;
end

