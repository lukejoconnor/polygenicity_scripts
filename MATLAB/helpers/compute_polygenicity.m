function Pi = compute_polygenicity(x,w,f,finv)
%compute_polygenicity evaluates the function
% Pi_f(x,w) = h^2 / f^-1 (1/h^2 w_1*f(x_1)+...+w_k*f(x_k))
% where h^2 = w_1+...+w_k and f:(0,infty)->(0,infty) is a continuous
% function with inverse finv.
%   For example:
%   \Pi_{entropy}: f = @log and finv=@exp
%   \Pi_{effective}: f=@(x)x and finv=f
%   \Pi_{softmax}: f=@(x)exp(-1./x) and finv=@(x)-1./log(x)
%   \Pi_{softmax} alternative implementation that avoids overflow issues:
%       f=@(x)exp(1/max(x(:)) - 1./x) and finv=@(y,xmax)-1./(-1/xmax + log(y))
%   \Pi_0 [not recommended]: f=@inv and finv=@inv
% 
% Can be applied to matrix-valued x and w (e.g., the jackknife output of FMR),
% in which case it computes polygenicity for each row of the matrix.

h2 = sum(w,2);
w = w./h2;
y = sum(w .* f(x),2);
if nargin(finv) == 1
    Pi = h2 ./ finv(y);
else
    % To avoid overflow issues
    Pi = h2 ./ finv(y, max(x(:)));
end
end