function metrics = evaluate_front(F, problem)
% Distribution and convergence diagnostics for a bi-objective ZDT front.
    F = F(nd_idx(F),:);
    F = sortrows(F, 1);
    if size(F,1) < 2
        metrics = struct('igd',inf,'spacing_cv',inf,'max_gap',inf, ...
                         'coverage',0,'count',size(F,1));
        return
    end

    t = linspace(0,1,1001)';
    switch upper(problem)
        case 'ZDT1'
            P = [t, 1-sqrt(t)];
        case 'ZDT2'
            P = [t, 1-t.^2];
        otherwise
            error('evaluate_front supports ZDT1 and ZDT2.');
    end
    D = zeros(size(P,1),1);
    for i = 1:size(P,1)
        D(i) = min(sqrt(sum((F-P(i,:)).^2,2)));
    end
    gaps = sqrt(sum(diff(F,1,1).^2,2));
    metrics.igd        = mean(D);
    metrics.spacing_cv = std(gaps) / max(mean(gaps),eps);
    metrics.max_gap    = max(gaps);
    metrics.coverage   = max(F(:,1)) - min(F(:,1));
    metrics.count      = size(F,1);
end
