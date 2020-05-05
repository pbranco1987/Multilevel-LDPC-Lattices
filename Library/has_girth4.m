function g4 = has_girth4(H,j)
persistent O
if isempty(O) || j == 1
    O = sparse(size(H,1),size(H,1));
end
Ij = find(H(:,j))';
g4 = false;
for i = 1:length(Ij)
    for k = 1:i-1
        if O(Ij(i),Ij(k)) > 0
            g4 = true;
            return
        else
            O(Ij(i),Ij(k)) = 1;
        end
    end
end