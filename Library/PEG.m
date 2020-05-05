function H = PEG(dv,m)
% progressive edge-growth. This is in line with the algorithm presented
% in "Progressive edge-growth Tanner Graphs" by Hu, Arnold, and
% Eleftheriou. The idea is to maximize the neighborhood, that is, the depth
% between nodes.
    dv = sort(dv);
    if max(dv) > m
        error('Number of rows smaller than maximum degree.')
    end
    n = length(dv);
    H = spalloc(m,n,sum(dv));
    dc = zeros(m,1);
    for j = 1:n
        for k = 1:dv(j)
            I = (1:m);
            % restrict to checks of maximum distance to symbol j
            D = get_distances(H,j);
            I = I(D(I) == max(D(I)));
            % restrict to checks of lowest degree
            I = I(dc(I) == min(dc(I)));
            % if multiple remain, choose one at random
            i = I(randi(length(I)));
            H(i,j) = 1;
            dc(i) = dc(i) + 1;
        end
        display_progress(j,n);
    %     if has_girth4(H,j)
    %         disp(''); disp(['Girth 4 at n = ' num2str(j)])
    %         H = H(:,1:j);
    %         return
    %     end
        %spy(H); figure(gcf); pause(0.001)
    end
    %% Check whether the last m=n-k columns of H form a full rank matrix.
    H_sub = full(H);
    H_sub = gf2rref(H_sub(:,n-m+1:end)); % Harvests the last m columns of H and applies reduced echelon form in gf(2).
    rank_sub = rank(full(real(H_sub))); % Calculates the rank of this matrix
    is_fullrank = rank_sub == size(H_sub,2); % Full rank or not
    if ~is_fullrank
        fset = 1:n;
        pp = make_invertible(H);
        dp = setdiff(fset,pp);
        H = [H(:,dp) H(:,pp)];
    end
end