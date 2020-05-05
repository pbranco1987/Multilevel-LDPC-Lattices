function H = EPEG_deterministic(Dv,M)

% this function builds a nested H parity-check matrix. The Dv input
% tells how nesting is supposed to be done and M indicates the number of 
% rows each sub parity-check matrix must have, as it is a matrix of degrees,
% where each row tells the degree of the respective submatrix. E-PEG's 
% algorithm is described in detail in "Low-Density Parity-Check Lattices: 
% Construction and Decoding Analysis" by Sadeghi, Banihashemi, and Panario

% This algorithm is deterministic, that is, for a given set of parameters
% Dv and M, the output matrix H will always be the same.

% Dv is a matrix with with L rows
% M is a vector with the number of parity-checks of the L levels
% M = [m0 m1 ... mL]

[L,n] = size(Dv);
H = spalloc(M(1),n,sum(Dv(1,:)));
Dv = flipud(Dv);
M = fliplr(M);

for j = 1:n
    if mod(j,100)==0 
        fprintf('(%d/%d)\n',j,n)
    end
    for ell = 1:L
        if ell == 1
            d = 0;
            m = 0;
        else
            d = Dv(ell-1,j);
            m = M(ell-1);
        end
        for k = d+1:Dv(ell,j)
            dc = sum(H,2);
            D = get_distances(H,j);
            S = m+1:M(ell);
            S = intersect(S, find(D == max(D(S))));
            S = intersect(S, find(dc == min(dc(S))));
            i = S(1);
            H(i,j) = 1;
        end
    end
end
