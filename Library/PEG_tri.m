function H = PEG_tri(dv,m,g)
% this function does the triangular progressive edge growth construction
% in "Progressive edge-growth Tanner graphs" by Hu, Arnold, and Eleftheriou
% forcing the parity-check matrix to be triangular

% progressive-edge growth
if ~exist('g','var') || isempty(g)
  g = m;
end

n = length(dv);
H = spalloc(m,n,sum(dv));
fprintf('Matriz PEG Triangular\n');
for j = 1:n
    if mod(j,100)==0
        fprintf('(%d/%d)\n',j,n)
    end
    for k = 1:dv(j)
        %Se está na diagonal, aloca o primeiro bit na diagonal
        if (j <= m-g) && (k == 1)
            H(g+j,j) = 1;
            continue
        end
        
        dc = sum(H,2);
        D = get_distances(H,j);
        
        %As possíveis linhas (de 1 até j+h-1 (se na diagonal)
        %                    (de 1 até m (depois da diagonal)
        S = (1:min(j+g-1,m)); 

        if ~isempty(S)    
            %----
            S = intersect(S, find(D==max(D(S))));
            S = intersect(S, find(dc==min(dc(S))));
            i = S(randi(length(S)));            
            %----          
            H(i,j) = 1;
        end
    end
    if has_girth4(H,j)
        disp(''); disp(['Girth 4 at n = ' num2str(j)])
        %H = H(:,1:j);
        %return
    end
end
