% Deterministic fixture for MATLAB/helpers/compute_polygenicity.m.

addpath('../MATLAB/helpers')

sigma2 = [0.05, 0.10];
omega = [0.25, 0.75];
h2 = 0.4;

entropy_expected = h2 * exp(sum(omega .* log(1 ./ sigma2)));
effective_expected = h2 / sum(omega .* sigma2);
softmax_expected = -h2 * log(sum(omega .* exp(-1 ./ sigma2)));

assert(abs(compute_polygenicity(sigma2, omega, h2, 'entropy') ...
    - entropy_expected) < 1e-12)
assert(abs(compute_polygenicity(sigma2, omega, h2, 'effective') ...
    - effective_expected) < 1e-12)
assert(abs(compute_polygenicity(sigma2, omega, h2, 'softmax') ...
    - softmax_expected) < 1e-12)

% The previous generic API remains available.
x = sigma2 / h2;
legacy_entropy = compute_polygenicity(x, omega, @log, @exp);
assert(abs(legacy_entropy - entropy_expected) < 1e-12)

disp('compute_polygenicity fixture passed')
