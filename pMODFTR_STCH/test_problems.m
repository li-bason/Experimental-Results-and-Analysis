function F = test_problems(name, x)
% ZDT benchmark suite (bi-objective, x in [0,1]^n, n >= 2)
% Usage: F = test_problems('ZDT1', x)
    x = x(:);
    n = length(x);
    switch upper(name)

        case 'ZDT1'
            % Convex Pareto front: f2 = 1 - sqrt(f1), f1 in [0,1]
            f1 = x(1);
            g  = 1 + 9/(n-1) * sum(x(2:end));
            f2 = g * max(0, 1 - sqrt(f1/g));
            F  = [f1; f2];

        case 'ZDT2'
            % Non-convex Pareto front: f2 = 1 - f1^2
            f1 = x(1);
            g  = 1 + 9/(n-1) * sum(x(2:end));
            f2 = g * max(0, 1 - (f1/g)^2);
            F  = [f1; f2];

        case 'ZDT3'
            % Disconnected Pareto front (5 segments)
            f1 = x(1);
            g  = 1 + 9/(n-1) * sum(x(2:end));
            f2 = g * (1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1));
            F  = [f1; f2];

        otherwise
            error('test_problems: unknown problem "%s".\n  Available: ZDT1, ZDT2, ZDT3', name);
    end
end
