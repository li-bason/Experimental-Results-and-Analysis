
function v = stch_val(f, lambda, z_star, mu)
% Smooth Tchebycheff via log-sum-exp (numerically stable)
%   STCH_mu(x; lambda, z*) = mu * log( sum_i exp( lambda_i*(f_i-z_i*)/mu ) )
%   Converges to Tchebycheff = max_i lambda_i*(f_i-z_i*) as mu -> 0
    y    = lambda(:) .* (f(:) - z_star(:));
    y_mx = max(y);
    v    = y_mx + mu * log( sum( exp((y - y_mx) / mu) ) );
end