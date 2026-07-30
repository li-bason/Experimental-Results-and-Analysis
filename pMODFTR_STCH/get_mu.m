
%% ==================================================================
%  LOCAL HELPER FUNCTIONS  (pure; no shared workspace)
%% ==================================================================

function mu = get_mu(opts, Delta, q)
% mu*log(q) <= c_mu*Delta^2 matches the function-value error order in
% Assumption 3; hence mu vanishes automatically with the trust region.
    if nargin < 3, q = 2; end
    if strcmp(opts.mu_mode, 'dynamic')
        mu = min(opts.mu_max, opts.mu_factor * opts.c_mu * Delta^2 / max(log(q),1));
        mu = max(mu, 1e-8);
    else
        mu = max(opts.mu_fixed, 1e-8);
    end
end
