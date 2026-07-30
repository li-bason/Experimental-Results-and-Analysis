function [grd, g0] = lin_model(x0, pts, gvals, Delta, p, n)
% Linear interpolation model m(x) = g0 + grd'*(x - x0)
% x0 is always included in the fitting subset [Fix F5].
    k = size(pts, 1);

    % Locate x0 in pts (it is always present after geometry repair)
    d0        = sqrt(sum((pts - x0').^2, 2));
    [~, id0]  = min(d0);

    % Probabilistic subset selection
    if rand() < p  ||  k <= n+2
        use = (1:k)';
    else
        nu  = max(n+2, round(k*p));
        use = randperm(k, min(nu,k))';
        if ~ismember(id0, use)
            use(1) = id0;           % force x0 into subset [F5]
        end
    end

    Xu = pts(use,:);    % ku x n
    gu = gvals(use);    % ku x 1
    ku = size(Xu, 1);

    % Scale displacements by Delta for numerical conditioning
    S = (Xu - x0') / Delta;    % ku x n
    A = [ones(ku,1), S];       % ku x (n+1)

    % Min-norm least-squares (handles underdetermined case)
    c   = pinv(A) * gu;        % n+1 x 1
    g0  = c(1);
    grd = c(2:end) / Delta;    % rescale gradient to original coordinates
end