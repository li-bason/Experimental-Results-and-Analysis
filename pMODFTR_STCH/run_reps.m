
%% =========================================================================
%  LOCAL FUNCTION  (requires MATLAB R2016b or later)
%% =========================================================================
function hv_mean = run_reps(prob, opts, n_rep, fev_grid)
% Run n_rep independent repetitions of pMODFTR_STCH and return
% the mean hypervolume interpolated onto fev_grid.
    buf = zeros(n_rep, numel(fev_grid));
    for r = 1:n_rep
        rng(r * 1000); % paired random numbers across configurations
        res  = pMODFTR_STCH(prob, opts);
        fh   = res.fev_hist;
        hh   = res.hv_hist;
        % Deduplicate fev values (safe guard for interp1)
        [fh_u, ia] = unique(fh, 'last');
        hh_u       = hh(ia);
        % Clamp grid to observed range, then interpolate
        fg = min(fev_grid, max(fh_u));
        buf(r,:) = interp1(fh_u, hh_u, fg, 'previous', hh_u(1));
    end
    hv_mean = mean(buf, 1);
end
