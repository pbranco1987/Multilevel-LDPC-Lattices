function [HH,p,A] = PEG_CS(H,mm)
% check splitting construction using progressive edge-growth criteria
% (maximum depth). It takes an H matrix and applies the check splitting 
% procedure to achieve matrix HH. The check-splitting construction using 
% PEG (progressive-edge growth) is explained "Multilevel LDPC Lattices with
% Efficient Encoding and Decoding and a Generalization of Construction D'" by 
% Paulo Branco and Danilo Silva

[m,n] = size(H);
HH = spalloc(mm,n,sum(sum(H)));

p = get_parents(H,mm);

for j = 1:n
%     if mod(j,100)==0
%         fprintf('(%d/%d)\n',j,n)
%     end
    I = find(H(:,j))';
    for i = I
        dcnew = sum(HH,2);
        D = get_distances(HH,j);
        S = find(p==i);
        S = intersect(S, find(D == max(D(S))));
        S = intersect(S, find(dcnew == min(dcnew(S))));
        k = S(randi(length(S)));
        HH(k,j) = 1;
    end
    display_progress(j,n);
%     if has_girth4(HH,j)
%         disp(''); disp(['Girth 4 at n = ' num2str(j)])
%         HH = HH(:,1:j);
%         return
%     end
end

A = double(ones(m,1)*p(:)' == (1:m)'*ones(1,mm));
if ~all(all(A*HH==H))
    error('A*HH ~= H')
end

function p = get_parents(H,mm)
m = size(H,1);
p = [1:m NaN*ones(1,mm-m)];
for k = m+1:mm
    dc = sum(H,2); % acho que poderia estar fora do loop
    count = histcounts(p,'BinMethod','integers')';
    dcnorm = dc./count;
    S = find(dcnorm == max(dcnorm));
    i = S(randi(length(S)));
    p(k) = i;
end
%p = repmat(1:m,1,ceil(mm/m)); p = p(1:mm);
