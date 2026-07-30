function idx = nd_idx(F)
% Indices of non-dominated rows in F  (m x q), O(m^2)
    m   = size(F, 1);
    dom = false(m, 1);
    for i = 1:m
        for j = 1:m
            if i ~= j  &&  all(F(j,:) <= F(i,:))  &&  any(F(j,:) < F(i,:))
                dom(i) = true;
                break
            end
        end
    end
    idx = find(~dom);
    if isempty(idx), idx = (1:m)'; end
end