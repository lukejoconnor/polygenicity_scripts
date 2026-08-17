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

% One-component inputs retain one row per replicate.
sigma2_rows = [0.1; 0.2];
omega_rows = [1; 1];
h2_rows = [0.3; 0.4];
assert(all(abs(compute_polygenicity(...
    sigma2_rows, omega_rows, h2_rows, 'effective') - [3; 2]) < 1e-12))

% Reciprocal-first formulas must not overflow on valid subnormal inputs.
tiny = eps(0);
for measure = {'entropy', 'effective', 'softmax'}
    assert(compute_polygenicity(tiny, 1, tiny, measure{1}) == 1)
end

% Zero-weight components must not set the numerical scaling.
sigma2_zero_weight = [tiny, realmax];
omega_zero_weight = [1, 0];
for measure = {'entropy', 'effective', 'softmax'}
    assert(compute_polygenicity(sigma2_zero_weight, omega_zero_weight, ...
        tiny, measure{1}) == 1)
end

disp('compute_polygenicity fixture passed')
