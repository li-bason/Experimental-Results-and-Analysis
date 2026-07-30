function opts = fill_opts(opts, prob)
% Fill missing fields with sensible defaults
    defs = { 'N',10;  'Delta0',0.3;  'Dmin',1e-7;  'Dmax',1.0; ...
             'eta1',0.1;  'eta2',0.75;  'g1',0.5;  'g2',2.0; ...
             'p',0.9;  'mu_mode','dynamic';  'mu_fixed',0.1; ...
             'c_mu',1.0;  'mu_max',5.0;  'max_fevals',500; ...
             'gamma',2.0; 'model_max_points',30; 'bad_model_fraction',0.35; ...
             'mu_factor',10.0 };
    for k = 1:size(defs,1)
        if ~isfield(opts, defs{k,1})
            opts.(defs{k,1}) = defs{k,2};
        end
    end
    if ~isfield(opts,'ref_pt')
        opts.ref_pt = 1.1 * ones(1, prob.q);
    end
end
