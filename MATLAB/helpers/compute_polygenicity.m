function Pi = compute_polygenicity(sigma2, omega, h2_or_f, measure_or_finv)
%COMPUTE_POLYGENICITY Evaluate Equation 15 of O'Connor & Sella.
%
%   Pi = compute_polygenicity(sigma2, omega, h2, measure)
%
%   Pi = h2 * f^-1(sum_k omega_k * f(1/sigma2_k)),
%
% where:
%   sigma2(k) is the FMR component variance sigma_k^2 in phenotypic-
%       variance units (not divided by h2);
%   omega(:,k) is the normalized fraction of heritability assigned to
%       component k, with sum_k omega(:,k) = 1; and
%   h2 is total SNP heritability.  Thus the absolute heritability weights
%       in Equation 15 are w(:,k) = h2 .* omega(:,k).
%
% measure selects the generator f that acts on component polygenicity
% 1/sigma_k^2:
%   'entropy':   f(u) = log(u),       f^-1(y) = exp(y)
%   'effective': f(u) = 1/u,         f^-1(y) = 1/y
%   'softmax':   f(u) = exp(-u),     f^-1(y) = -log(y)
%
% omega may contain one row per jackknife replicate. sigma2 may be a row
% vector shared across replicates or a matrix of the same size as omega.
%
% Backward compatibility: calls using the previous
% compute_polygenicity(x,w,f,finv) interface retain their original behavior.

if isa(h2_or_f, 'function_handle') && ...
        isa(measure_or_finv, 'function_handle')
    f = h2_or_f;
    finv = measure_or_finv;
    weight_sum = sum(omega, 2);
    normalized_weights = omega ./ weight_sum;
    mean_f = sum(normalized_weights .* f(sigma2), 2);
    if nargin(finv) == 1
        Pi = weight_sum ./ finv(mean_f);
    else
        Pi = weight_sum ./ finv(mean_f, max(sigma2(:)));
    end
    return
end

h2 = h2_or_f;
measure = measure_or_finv;

if isvector(sigma2)
    sigma2 = reshape(sigma2, 1, []);
end
if isvector(omega)
    omega = reshape(omega, 1, []);
end
if size(sigma2, 2) ~= size(omega, 2)
    error('sigma2 and omega must have the same number of columns.');
end
if size(sigma2, 1) == 1 && size(omega, 1) > 1
    sigma2 = repmat(sigma2, size(omega, 1), 1);
elseif size(sigma2, 1) ~= size(omega, 1)
    error(['sigma2 must have one row or the same number of rows as ' ...
        'omega.']);
end
if any(~isfinite(sigma2(:)) | sigma2(:) <= 0)
    error('sigma2 must contain positive finite values.');
end
if any(~isfinite(omega(:)) | omega(:) < 0)
    error('omega must contain nonnegative finite values.');
end

omega_sum = sum(omega, 2);
if any(abs(omega_sum - 1) > 1e-10)
    error('Each row of omega must sum to one.');
end
omega = omega ./ omega_sum;
h2 = h2(:);
if isscalar(h2)
    h2 = repmat(h2, size(omega, 1), 1);
elseif length(h2) ~= size(omega, 1)
    error('h2 must be scalar or have one value per row of omega.');
end
if any(~isfinite(h2) | h2 <= 0)
    error('h2 must contain positive finite values.');
end
if ~(ischar(measure) || (isstring(measure) && isscalar(measure)))
    error('measure must be entropy, effective, or softmax.');
end

switch lower(char(measure))
    case 'entropy'
        mean_f = sum(omega .* log(1 ./ sigma2), 2);
        Pi = h2 .* exp(mean_f);
    case 'effective'
        mean_f = sum(omega .* sigma2, 2);
        Pi = h2 ./ mean_f;
    case 'softmax'
        % Stable evaluation of -h2*log(sum omega*exp(-1/sigma2)).
        % Direct exponentiation underflows for the FMR component scales.
        log_terms = log(omega) - 1 ./ sigma2;
        row_max = max(log_terms, [], 2);
        log_mean_f = -inf(size(row_max));
        finite_rows = isfinite(row_max);
        log_mean_f(finite_rows) = row_max(finite_rows) + log(sum(exp(...
            log_terms(finite_rows,:) - row_max(finite_rows)), 2));
        Pi = -h2 .* log_mean_f;
    otherwise
        error('Unknown measure: %s', measure);
end
end
