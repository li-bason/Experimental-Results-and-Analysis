function hv = hv2d(F, ref)
% Exact 2-D hypervolume indicator via sweep-line (O(k log k))
    nd_i = nd_idx(F);
    Fnd  = F(nd_i,:);
    ok   = Fnd(:,1) < ref(1)  &  Fnd(:,2) < ref(2);
    Fnd  = Fnd(ok,:);
    if isempty(Fnd), hv = 0; return; end
    Fnd = sortrows(Fnd, 1);    % ascending f1 <=> descending f2 (non-dom.)
    k   = size(Fnd, 1);
    hv  = 0;
    for i = 1:k
        if i < k
            dx = Fnd(i+1,1) - Fnd(i,1);
        else
            dx = ref(1) - Fnd(k,1);
        end        
        dy = ref(2) - Fnd(i,2);
        if dx > 0 && dy > 0
            hv = hv + dx*dy;
        end
    end
end