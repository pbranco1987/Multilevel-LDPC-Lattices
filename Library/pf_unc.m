function pf = pf_unc(sigma,n)
% Frame error probability of ZZ^n constellation

pf = n*2*qfunc(1/(2*sigma));    % low probability approximation
if pf > 1e-6
    pf = 1-(1-2*qfunc(1/(2*sigma)))^n;  % exact value
end
