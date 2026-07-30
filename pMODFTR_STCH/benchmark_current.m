clc; clear;
prob.n=5; prob.q=2; prob.lb=zeros(5,1); prob.ub=ones(5,1);
opts.N=15; opts.max_fevals=800; opts.Delta0=0.3; opts.Dmin=1e-7;
opts.Dmax=1; opts.eta1=0.1; opts.eta2=0.75; opts.g1=0.5; opts.g2=2;
opts.mu_max=0.5; opts.c_mu=1; opts.p=0.9; opts.mu_mode='dynamic';
opts.ref_pt=[1.1 1.1];
for pi=1:2
    name=sprintf('ZDT%d',pi); prob.f=@(x)test_problems(name,x);
    for r=1:5
        rng(10000+100*pi+r); res=pMODFTR_STCH(prob,opts);
        m=evaluate_front(res.F,name);
        fprintf('%s r%d IGD %.8f CV %.8f GAP %.8f COV %.8f K %d HV %.8f\n', ...
            name,r,m.igd,m.spacing_cv,m.max_gap,m.coverage,m.count,hv2d(res.F,opts.ref_pt));
    end
end
