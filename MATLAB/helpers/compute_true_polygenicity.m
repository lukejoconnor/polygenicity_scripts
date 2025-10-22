function pi = compute_true_polygenicity(perSNPh2,marginal_variance,f,finv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h2 = sum(perSNPh2);
nz = perSNPh2 > 0;
y = 1/h2 * sum(perSNPh2(nz) .* f(marginal_variance(nz)));
if nargin(finv) == 1
    pi = h2 ./ finv(y);
else
    pi = h2 ./ finv(y, max(marginal_variance));
end
end