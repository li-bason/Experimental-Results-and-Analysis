function Lambda = make_weights(N, q)
% Uniform weight vectors on the (q-1)-simplex
    if q == 2
        t      = linspace(0,1,N)';
        Lambda = [t, 1-t];
    else
        Lambda = zeros(N,q);
        for i = 1:N
            w = -log(rand(q,1) + 1e-12);
            Lambda(i,:) = (w/sum(w))';
        end
    end
    Lambda = max(Lambda, 1e-4);
    Lambda = Lambda ./ sum(Lambda, 2);
end