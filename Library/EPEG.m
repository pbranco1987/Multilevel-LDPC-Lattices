function H = EPEG(Dv,M)
% this function builds a nested H parity-check matrix. The Dv input
% tells how nesting is supposed to be done and M indicates the number of 
% rows each sub parity-check matrix must have, as it is a matrix of degrees,
% where each row tells the degree of the respective submatrix. E-PEG's 
% algorithm is described in detail in "Low-Density Parity-Check Lattices: 
% Construction and Decoding Analysis" by Sadeghi, Banihashemi, and Panario

% This algorithm is random, that is, for a given set of parameters
% Dv and M, the output matrix H may be different as the algorithm is run
% multiple times

% Dv is a matrix with with L rows
% M is a vector with the number of parity-checks of the L levels
% M = [m0 m1 ... mL]

[L,n] = size(Dv);
H = spalloc(M(1),n,sum(Dv(1,:)));
dc = zeros(M(1),1);
M = M(:);
M = flipud(M);
if ~issorted(M)
    error('Vector M must be in decreasing order.')
end
Dv = flipud(Dv);
Dv = sort(Dv,2);
if any(max(Dv,[],2) > M)
    error('Number of rows smaller than maximum degree.')
end
for j = 1:n
    for ell = 1:L
        if ell == 1
            d = 0;
            m = 0;
        else
            d = Dv(ell-1,j);
            m = M(ell-1);
        end
        for k = d+1:Dv(ell,j)
            I = m+1:M(ell);
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
    end
    display_progress(j,n)
    %spy(H); figure(gcf); pause(0.001)
end
