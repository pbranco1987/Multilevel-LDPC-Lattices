function parity_positions = make_invertible(A,varargin)
A = logical(A);
[m,n] = size(A);
if nargin==1
    g = m;
elseif nargin==2
    if varargin{1}<=m
        g = varargin{1};
    else
        error('Gap must be less than or equal to m');
    end
end
Ainv = [A(g+1:end,:); A(1:g,:)];
for jdx = 1:m-g
    ind = find(Ainv(m-g+1:end,jdx)) + m-g;
    if ~isempty(ind)
        Ainv(ind,:) = xor(Ainv(ind,:),repmat(Ainv(jdx,:),length(ind),1));  
    end
end
pc = m-g+1;
parity_positions = [1:m-g zeros(1,g)];
for jdx = m-g+1:n
    indx = find(Ainv(pc:end,jdx))+pc-1;
    indx = indx(indx>=pc);
    if ~isempty(indx)
        ind = indx(2:end);
        pv = indx(1);
        Ainv(ind,:) = xor(Ainv(ind,:),repmat(Ainv(pv,:),length(ind),1));  
        Ainv([pc pv],:) = Ainv([pv pc],:);
        parity_positions(pc) = jdx;
        pc = pc + 1;
    end
    if pc>m
        break
    end
end






% end
% c = 1;
% pivot_counter = g;
% for jdx = g:-1:1
%     indices = find(Ainv(:,m-g+c));
%     indx = indices(indices>=pivot_counter);
%     if ~isempty(indx)
%         ind = indx(2:end);
%         pv = indx(1);
%         Ainv(ind,m-g+c:end) = xor(Ainv(ind,m-g+c:end),repmat(Ainv(pv,m-g+c:end),length(ind),1));
%         Ainv([jdx pv],:) = Ainv([pv jdx],:);
%     end
%     c = c + 1;
%     pivot_counter = pivot_counter - 1;
% end